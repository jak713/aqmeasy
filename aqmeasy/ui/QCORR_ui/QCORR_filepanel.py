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
from PySide6.QtCore import QMimeData, Qt, Signal
from PySide6.QtGui import QDrag, QIcon
from aqmeasy.ui.stylesheets import stylesheets
from aqmeasy.ui.icons import Icons

from aqmeasy.controllers.QCORR_controller import FileController


class FilePanel(QWidget):
    """File browser, batch file selection, drag and drop, file status"""
    
    # fileSelected = Signal(str)

    def __init__(self, model):
        super().__init__()
        self.model = model
        self.controller = FileController(model, self)
        # self.view_panel = view_panel
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

        run_button = QPushButton("Run QCORR")
        run_button.setStyleSheet(stylesheets.RunButton)

        run_button.clicked.connect(self.controller.run_qcorr)
        # run_button.clicked.connect(lambda: self.view_panel.results_view.setRootIndex(root for root,dirs in self.model.__get__w_dir_main__()))

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
        self.file_view.itemSelectionChanged.connect(lambda: self._on_file_selection_changed(self.file_view.currentItem().toolTip()))

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
        for file in filenames:
            # check if file is already in the list
            if not any(item.toolTip() == file for item in self.file_view.findItems("*", Qt.MatchFlag.MatchWildcard)):
                item = QListWidgetItem(os.path.basename(file))
                item.setToolTip(file)
                self.file_view.addItem(item)

    # def _clear_file_list(self):
    #     self.file_view.clear()
    #     self.view_panel.file_viewer.clear()

    def _on_file_selection_changed(self, selected_item):
        self.model.currently_selected_file = selected_item
        self.model.currentlySelectedFileChanged.emit(selected_item)
        print(self.model.__get__currently_selected_file__())
        

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