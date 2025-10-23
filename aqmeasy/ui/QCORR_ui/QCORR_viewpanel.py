import os
from PySide6.QtWidgets import (
    QWidget
)

class ViewPanel(QWidget):


    def __init__(self):
        super().__init__()

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