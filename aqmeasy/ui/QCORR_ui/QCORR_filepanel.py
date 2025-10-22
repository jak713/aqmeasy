import os
from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QFileDialog,
    QLabel,
    QProgressBar,
    QFileSystemModel,
    QTreeView,
    QListWidget, 
    QListWidgetItem
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
        selecting_layout = QHBoxLayout()
        layout.addLayout(selecting_layout)

        select_files = QPushButton("Select Files")
        select_files.setStyleSheet(stylesheets.QPushButton)
        select_files.clicked.connect(self.get_filenames)
        selecting_layout.addWidget(select_files)

        select_directory = QPushButton("Select Folder")
        select_directory.setStyleSheet(stylesheets.QPushButton)
        select_directory.clicked.connect(self.get_directory)
        selecting_layout.addWidget(select_directory)

        # file list view

        # self.files = QFileSystemModel(self)
        # self.files.setNameFilters(['*.log', '*.out'])
        # self.files.setNameFilterDisables(False)
        # self.files.setRootPath('')

        # self.file_view = QTreeView()
        # self.file_view.setStyleSheet(stylesheets.QTreeView)
        # # hide kind column
        # self.file_view.setModel(self.files)
        # self.file_view.setColumnWidth(0,200)
        # self.file_view.setColumnHidden(2, True)
        # layout.addWidget(self.file_view)

    def get_filenames(self):
        filters = ';;'.join(FILE_FILTERS)
        filenames,selected_filter = QFileDialog.getOpenFileNames(self, filter=filters)
        ######
        # need to store filenames on model
        ######

    def get_directory(self):
        folder_path = QFileDialog.getExistingDirectory(caption="",dir="",options=QFileDialog.ShowDirsOnly)
        if folder_path:
            self.display_folder_contents(folder_path)

    def display_folder_contents(self, folder_path):
        """Display files in the selected folder in selected_paths QLabel, that match the file filters."""
        matched_files = []
        folders = folder_path if isinstance(folder_path, list) else [folder_path]

        for single_folder in folders:
            for root, dirs, files in os.walk(single_folder):
                for file in files:
                    if file.lower().endswith(('.log', '.out')):
                        matched_files.append(os.path.join(file))
                        self.file_view.setRootIndex(self.files.index(root))


    def dragEnterEvent(self, event):
        urls = event.mimeData().urls()
        if not urls:
            event.ignore()
            return
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