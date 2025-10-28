import os
from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QGroupBox, 
    QFileSystemModel,
    QTreeView,
    QTextBrowser,
)
# from PySide6.QtCore import Qt, QRunnable, Slot
from PySide6.QtGui import QTextOption
from aqmeasy.ui.icons import Icons
from aqmeasy.controllers.QCORR_ViewController import ViewController

class ViewPanel(QWidget):

    def __init__(self,model):
        super().__init__()
        self.model = model
        self.controller = ViewController(model, self)
        self.init_ui()
        self.update_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)

        view_group = QGroupBox("File Contents")
        view_layout = QVBoxLayout()
        view_group.setLayout(view_layout)

        self.file_viewer = QTextBrowser()
        self.file_viewer.setFontFamily("Menlo")
        self.file_viewer.setWordWrapMode(QTextOption.WrapMode.NoWrap)
        view_layout.addWidget(self.file_viewer)

        results_group = QGroupBox("QCORR Results")
        results_layout = QVBoxLayout()
        results_group.setLayout(results_layout)

        self.processed_files = QFileSystemModel()
        self.processed_files.setRootPath('')
        self.results_view = QTreeView()
        self.results_view.setModel(self.processed_files)
        self.results_view.setRootIsDecorated(False)
        results_layout.addWidget(self.results_view)

        layout.addWidget(view_group)
        layout.addWidget(results_group)

    # def display_folder_contents(self, folder_path):
    #     """Display files in the selected folder in selected_paths QLabel, that match the file filters."""
    #     matched_files = []
    #     folders = folder_path if isinstance(folder_path, list) else [folder_path]

    #     for single_folder in folders:
    #         for root, dirs, files in os.walk(single_folder):
    #             for file in files:
    #                 if file.lower().endswith(('.log', '.out')):
    #                     matched_files.append(os.path.join(file))
    #                     self.file_view.setRootIndex(self.files.index(root))


    def update_ui(self):
        """Update the UI elements based on the model's state."""
        w_dir = self.model.__get__w_dir_main__()
        if w_dir:
            self.results_view.setRootIndex(self.processed_files.index(w_dir))

    #TODO:
    # Implement search functionality within the file viewer.
    