import os
from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QLabel,
    QPlainTextEdit,
    QTextEdit,
    QGroupBox, 
    QFileSystemModel,
    QTreeView,
    QTextBrowser,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QTextOption
from aqmeasy.ui.icons import Icons

class ViewPanel(QWidget):

    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)

        view_group = QGroupBox("File Contents")
        view_layout = QVBoxLayout()
        view_group.setLayout(view_layout)

        self.file_viewer = QTextBrowser()
        self.file_viewer.setFontFamily("Courier")
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
        self.results_view.setAlternatingRowColors(True)
        results_layout.addWidget(self.results_view)

        layout.addWidget(view_group)
        layout.addWidget(results_group)

    def display_file_contents(self, file_path):
        """Display the contents of the selected file in the text viewer."""
        print(f"Displaying contents of file: {file_path}")
        try:
            with open(file_path, 'r') as file:
                content = file.read()
                character_count = len(content)
                
                if character_count > 10000000:
                    self.file_viewer.setText(f"File too large to display efficiently ({character_count} characters).\nThis will be dealt with in due course.")
                else:
                    self.file_viewer.setText(content)
                    
        except Exception as e:
            self.file_viewer.setText(f"Error reading {file_path}: {e}")

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