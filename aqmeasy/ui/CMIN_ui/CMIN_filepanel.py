import os
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
    QStyle,
)
from PySide6.QtCore import  Qt
from PySide6.QtGui import  QIcon, QColor, QBrush
from aqmeasy.ui.icons import Icons
from aqmeasy.ui.stylesheets import stylesheets

from aqmeasy.controllers.CMIN_controller import FileController



class FilePanel(QWidget):
    """File browser, batch file selection, drag and drop, file status"""

    def __init__(self, parent, model):
        super().__init__()
        self.cmin_parent = parent
        self.model = model
        self.controller = FileController(model, self)
        self.init_ui()
        self.setAcceptDrops(True)


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

        self.run_button = QPushButton("Run CMIN")
        self.run_button.setStyleSheet(stylesheets.RunButton)
        self.run_button.clicked.connect(self._on_run_clicked)
        selecting_layout.addWidget(self.run_button)
        self.set_running(False)

        clear_files = QPushButton()
        clear_files.setIcon(QIcon(Icons.trash))
        clear_files.clicked.connect(self.controller.clear_file_list)
        selecting_layout.addWidget(clear_files)



        # file list view
        self.file_view = QListWidget()
        # Selection updates model's currently selected file
        self.file_view.itemSelectionChanged.connect(lambda: self._on_file_selection_changed(self.file_view.currentItem().toolTip() if self.file_view.currentItem() else ""))

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

    def _on_run_clicked(self):
        """Start or stop CMIN depending on current state."""
        if self.cmin_parent is None:
            return

        if self.run_button.text() == "Stop CMIN":
            if hasattr(self.cmin_parent, 'stop_cmin'):
                self.cmin_parent.stop_cmin()
        else:
            if hasattr(self.cmin_parent, 'run_cmin'):
                self.cmin_parent.run_cmin()

    def set_running(self, running):
        """Update run button label and icon for running state."""
        if running:
            self.run_button.setText("Stop CMIN")
            self.run_button.setIcon(
                QApplication.style().standardIcon(QStyle.StandardPixmap.SP_MediaStop)
            )
        else:
            self.run_button.setText("Run CMIN")
            self.run_button.setIcon(
                QApplication.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay)
            )


    def display_selected_files(self, filenames):
        for file in filenames:
            # check if file is already in the list
            if not any(item.toolTip() == file for item in self.file_view.findItems("*", Qt.MatchFlag.MatchWildcard)):
                item = QListWidgetItem(os.path.basename(file))
                item.setToolTip(file)
                # Default color: white/default (input files)
                self.file_view.addItem(item)

    def add_output_file(self, filepath, status='success'):
        """Add output file with color coding
        
        Args:
            filepath: Path to output file
            status: 'success' (green), 'failed' (red), 'processing' (grey),
                    'partial' (orange), 'eliminated' (red)
        """
        # Check if already exists
        if any(item.toolTip() == filepath for item in self.file_view.findItems("*", Qt.MatchFlag.MatchWildcard)):
            return
        
        item = QListWidgetItem(os.path.basename(filepath))
        item.setToolTip(filepath)
        
        # Set background color based on status
        self._apply_item_status(item, status)
        
        self.file_view.addItem(item)

    def set_file_status(self, filepath, status):
        """Update the color status of an existing file item in the list."""
        for item in self.file_view.findItems("*", Qt.MatchFlag.MatchWildcard):
            if item.toolTip() == filepath:
                self._apply_item_status(item, status)
                break

    def _apply_item_status(self, item, status):
        """Colorize one list item according to processing status."""
        if status == 'success':
            item.setBackground(QBrush(QColor(46, 204, 113, 80)))
        elif status in ('failed', 'eliminated'):
            item.setBackground(QBrush(QColor(231, 76, 60, 80)))
        elif status == 'partial':
            item.setBackground(QBrush(QColor(243, 156, 18, 90)))
        elif status == 'processing':
            item.setBackground(QBrush(QColor(149, 165, 166, 80)))
        else:
            item.setBackground(QBrush())

    def _on_file_selection_changed(self, selected_item):
        if not selected_item or self.model.isSelectable == False:
            return
        self.model.currently_selected_file = selected_item
        self.model.currentlySelectedFileChanged.emit(selected_item)

        # Use load_molecules_from_files to load the selected file into the molecule viewer
        self.cmin_parent.molecule_viewer.load_molecules_from_files([selected_item])

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