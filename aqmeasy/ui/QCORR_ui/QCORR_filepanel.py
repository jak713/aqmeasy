import os
from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QFileDialog,
    QLabel,
    QListWidget, 
    QListWidgetItem,
    QGroupBox,
    QLineEdit,
)
from PySide6.QtCore import QMimeData, Qt, Signal
from PySide6.QtGui import QDrag
from aqmeasy.ui.stylesheets import stylesheets


FILE_FILTERS = [
    "Gaussian Output Files (*.log)",
    "ORCA Output Files (*.out)",
]

class FilePanel(QWidget):
    """File browser, batch file selection, drag and drop, file status"""

    # Space for signals when the model is ready

    def __init__(self):
        super().__init__()
        self.init_ui()
        self.setAcceptDrops(True)

    def init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)
        input_group = QGroupBox("Drop in files or folders below")
        input_group.setStyleSheet(stylesheets.QGroupBox)
        layout.addWidget(input_group)
        layout = QVBoxLayout()
        input_group.setLayout(layout)

        selecting_layout = QHBoxLayout()
        layout.addLayout(selecting_layout)

        select_files = QPushButton("Browse Files")
        select_files.setStyleSheet(stylesheets.QPushButton)
        select_files.clicked.connect(self.get_filenames)
        selecting_layout.addWidget(select_files)

        clear_files = QPushButton("Clear Files")
        clear_files.setStyleSheet(stylesheets.QPushButton)
        clear_files.clicked.connect(self.clear_file_list)
        selecting_layout.addWidget(clear_files)

        # file list view
        self.file_view = QListWidget()
        self.file_view.setStyleSheet(stylesheets.QListWidget)
        layout.addWidget(self.file_view)

        # output dir selection
        output_layout = QHBoxLayout()
        layout.addLayout(output_layout)
        output_label = QLabel("Output Directory:")
        output_label.setStyleSheet(stylesheets.QLabel)
        output_layout.addWidget(output_label)
        self.output_dir_label = QLineEdit("Not selected")
        self.output_dir_label.setStyleSheet(stylesheets.QLineEdit)
        self.output_dir_label.setReadOnly(True)
        output_layout.addWidget(self.output_dir_label)
        select_output_dir = QPushButton("Browse")
        select_output_dir.setStyleSheet(stylesheets.QPushButton)
        select_output_dir.clicked.connect(self.get_output_directory)
        output_layout.addWidget(select_output_dir)

    def get_filenames(self):
        filters = ';;'.join(FILE_FILTERS)
        filenames,selected_filter = QFileDialog.getOpenFileNames(self, filter=filters)
        ######
        # need to store filenames on model
        ######

    def clear_file_list(self):
        self.file_view.clear()

    def get_output_directory(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if directory:
            self.output_dir_label.setText(directory)
            ######
            # need to store output directory on model
            ######

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
        if not dirs:
            event.ignore()
            return
        folder_path = dirs[0] if len(dirs) == 1 else dirs

        if folder_path:
            self.display_folder_contents(folder_path)
        else:
            self.display_folder_contents(dirs)
        
        event.acceptProposedAction()