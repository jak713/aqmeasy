import os
from xml.parsers.expat import model
from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QLabel,
    QListWidget, 
    QListWidgetItem,
    QGroupBox,
    QLineEdit,
    QApplication,
    QSizePolicy,
    QStyle
)
from PySide6.QtCore import  Qt
from PySide6.QtGui import  QIcon, QBrush, QColor    
from aqmeasy.ui.stylesheets import stylesheets
from aqmeasy.ui.icons import Icons

from aqmeasy.controllers.QCORR_controllers.QCORR_controller import FileController
from aqmeasy.controllers.QCORR_controllers.QCORR_worker import QCORRWorker



class FilePanel(QWidget):
    """File browser, batch file selection, drag and drop, file status"""

    def __init__(self, parent, model):
        super().__init__()
        self.parent = parent
        self.model = model
        self.controller = FileController(model, self)
        self.init_ui()
        self.setAcceptDrops(True)

        self.model.filesChanged.connect(self._refresh_file_list)
        self.model.fileStatusChanged.connect(self._update_all_colors)


    def init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)
        input_group = QGroupBox("Drop in files or folders below")
        layout.addWidget(input_group)
        layout = QVBoxLayout()
        input_group.setLayout(layout)

        selecting_layout = QHBoxLayout()
        layout.addLayout(selecting_layout)

        select_files = QPushButton()
        select_files.setIcon(QIcon(Icons.file_open))
        select_files.clicked.connect(self.controller.open_file_dialog)
        selecting_layout.addWidget(select_files)

        run_button = QPushButton("Run QCORR")
        run_button.setStyleSheet(stylesheets.RunButton)

        run_button.clicked.connect(self.parent.run_qcorr)

        run_button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        run_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay))
        selecting_layout.addWidget(run_button)

        clear_files = QPushButton()
        clear_files.setIcon(QIcon(Icons.trash))
        clear_files.clicked.connect(self.controller.clear_file_list)
        selecting_layout.addWidget(clear_files)

        # file list view
        self.file_view = QListWidget()
        # Selection updates model's currently selected file
        self.file_view.itemSelectionChanged.connect(self._on_file_selection_changed)


        layout.addWidget(self.file_view)

        # output dir selection
        output_layout = QHBoxLayout()
        layout.addLayout(output_layout)
        
        output_label = QLabel("Output Directory:")
        output_layout.addWidget(output_label)

        self.output_dir_label = QLineEdit()
        self.output_dir_label.setPlaceholderText("Select output directory...")
        self.output_dir_label.setReadOnly(True)
        output_layout.addWidget(self.output_dir_label)

        select_output_dir = QPushButton()
        select_output_dir.setIcon(QIcon(Icons.folder_open))
        select_output_dir.clicked.connect(self.controller.select_output_directory)
        output_layout.addWidget(select_output_dir)


    def _display_selected_files(self, filenames):
        """Display files in the list widget."""
        for file in filenames:
            # check if file is already in the list
            if not any(item.toolTip() == file for item in self.file_view.findItems("*", Qt.MatchFlag.MatchWildcard)):
                item = QListWidgetItem(self._format_display_name(file))
                item.setToolTip(file)
                self.file_view.addItem(item)
                self._update_item_color(item)

    def _format_display_name(self, filepath):
        """
        Format display name to show file location structure.
        
        Examples:
            /path/QCORR_1/success/calc.log -> ✓ calc.log
            /path/QCORR_1/failed/imag_freq/calc.log -> ✗ [imag_freq] calc.log
        """
        filename = os.path.basename(filepath)
        dirname = os.path.dirname(filepath)
        
        # Check if it's in a success folder
        if dirname.endswith('success'):
            return f"✓ {filename}"
        
        # Check if it's in a failed subfolder
        if 'failed' in dirname:
            # Extract error type (last directory name before the file)
            parts = dirname.split(os.sep)
            if 'failed' in parts:
                failed_idx = parts.index('failed')
                if failed_idx < len(parts) - 1:
                    error_type = parts[failed_idx + 1]
                    return f"✗ [{error_type}] {filename}"
            return f"✗ {filename}"
        
        # Original file (before QCORR)
        return filename
    
    def _refresh_file_list(self, filenames):
        """
        Completely refresh the file list (called after QCORR reorganizes files).
        """
        # Clear existing items
        self.file_view.clear()
        
        # Add all files from model
        for filepath in filenames:
            item = QListWidgetItem(self._format_display_name(filepath))
            item.setToolTip(filepath)  # Store full path in tooltip
            self.file_view.addItem(item)
            self._update_item_color(item)
    
    def _update_item_color(self, item):
        """Update list item color based on file status."""
        filepath = item.toolTip()
        status = self.model.get_file_status(filepath)
        
        if status == 'success':
            item.setForeground(QBrush(QColor(0, 180, 0)))  # Green
        elif status.startswith('failed'):
            item.setForeground(QBrush(QColor(220, 50, 50)))  # Red
        else:  
            item.setForeground(QBrush(QColor(150, 150, 150)))  
    
    def _update_all_colors(self, statuses):
        """Update colors for all items when statuses change."""
        for i in range(self.file_view.count()):
            item = self.file_view.item(i)
            self._update_item_color(item)

    def _on_file_selection_changed(self):
        """Handle file selection change in the list widget."""
        if self.model.isSelectable == False:
            return
        
        current_item = self.file_view.currentItem()
        if current_item is None:
            return
        
        # Get full path from tooltip
        full_path = current_item.toolTip()
        
        self.model.currently_selected_file = full_path
        self.model.currentlySelectedFileChanged.emit(full_path)
    
    def _on_status_changed(self, statuses):
        """Update all item colors when statuses change."""
        for i in range(self.file_view.count()):
            item = self.file_view.item(i)
            self._update_item_color(item)

    def dragEnterEvent(self, event):
        urls = event.mimeData().urls()
        if not urls:
            event.ignore()
            return
        else:
            for url in urls:
                if not url.isLocalFile():
                    event.ignore()
                    return
        event.acceptProposedAction()

    def dropEvent(self, event):
        urls = event.mimeData().urls()
        dirs = [url.toLocalFile() for url in urls if url.isLocalFile() and os.path.isdir(url.toLocalFile())]
        if len(dirs) == 0 and len(urls) == 1 and os.path.isfile(urls[0].toLocalFile()):
            self.controller.open_drop_file(urls[0].toLocalFile())
            return
        else:
            self.controller.open_drop_folder(dirs)
        event.acceptProposedAction()