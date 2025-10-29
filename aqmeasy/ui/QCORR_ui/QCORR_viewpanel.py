import os
from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QGroupBox, 
    QFileSystemModel,
    QTreeView,
    QTextBrowser,
    QLineEdit,
    QPushButton,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QTextOption, QTextCharFormat, QColor

from aqmeasy.ui.icons import Icons
from aqmeasy.controllers.QCORR_controllers.QCORR_ViewController import ViewController

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
        self.results_view.setRootIsDecorated(True)
        self.results_view.setColumnWidth(0, 200)

        self.results_view.doubleClicked.connect(
            lambda index: self.display_file_content(
                self.processed_files.filePath(index)
            )
        )
        results_layout.addWidget(self.results_view)
        layout.addWidget(view_group)
        layout.addWidget(results_group)

    def update_ui(self):
        """Update the UI elements based on the model's state."""
        w_dir = self.model.__get__w_dir_main__()
        print(f"Updating UI with working directory: {w_dir}")
        if w_dir:
            # self.processed_files.setRootPath(w_dir)
            self.results_view.setRootIndex(self.processed_files.index(w_dir))

    def display_file_content(self, file_path):
        """Display the content of the selected file in the text viewer."""
        try:
            with open(file_path, 'r') as file:
                content = file.read()
                self.file_viewer.setPlainText(content)
        except Exception as e:
            self.file_viewer.setPlainText(f"Error reading file: {e}")

    def _display_selected_files(self, filenames):
        """Display the list of selected files in the file viewer."""
        display_text = "\n".join(filenames)
        self.file_viewer.setPlainText(display_text)

    def clear_file_viewer(self):
        """Clear the file viewer content."""
        self.file_viewer.clear()

# TODO
# class searchDialog