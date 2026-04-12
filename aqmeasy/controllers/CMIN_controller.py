from PySide6.QtWidgets import QFileDialog
from PySide6.QtCore import Signal, QObject, Slot, QRunnable
from aqme.cmin import cmin
from aqmeasy.utils import discover_aqme_result_files

import os

FILE_FILTERS = [
    "Structure-Data Files (*.sdf)",
    "XYZ Files (*.xyz)",
]

class FileController(QObject):
    """Controller for the CMIN module, connecting the model and the view."""
    def __init__(self, model, view):
        super().__init__()
        self.model = model
        self.view = view
        self.model.wdirChanged.connect(self.set_wdir)

    def open_file_dialog(self):
        """Opens a file dialog to select files and updates the model."""
        filters = ';;'.join(FILE_FILTERS)
        filenames, _ = QFileDialog.getOpenFileNames(self.view, filter=filters)
        if filenames:
            self.model.files = filenames
            self.model.filesChanged.emit(filenames)
            self.view.display_selected_files(filenames)
        
    def clear_file_list(self):
        """Clears the file list in the model and view."""
        self.model.isSelectable = False
        self.view.file_view.clear()
        self.model.files = []
        self.model.filesChanged.emit([])
        self.model.isSelectable = True

    def select_output_directory(self):
        """Opens a dialog to select the output directory and updates the model."""
        directory = QFileDialog.getExistingDirectory(self.view, "Select Output Directory")
        if directory:
            self.model.w_dir_main = directory
            self.model.wdirChanged.emit(directory)
            self.set_wdir(directory)

    def set_wdir(self, directory=None):
        if directory:
            self.view.output_dir_label.setText(directory)
        else:
            directory = self.model.w_dir_main
            self.view.output_dir_label.setText(directory)

    def open_drop_folder(self, dirs):
        """Processes dropped folders to find relevant files and updates the model."""
        all_files = []
        for dir_path in dirs:
            for root, _, files in os.walk(dir_path):
                for file in files:
                    if file.endswith(('.sdf', '.xyz')):
                        all_files.append(os.path.join(root, file))
        if all_files:
            self.model.files = all_files
            self.model.filesChanged.emit(all_files)
            self.view.display_selected_files(all_files)
            self.model.w_dir_main = os.path.dirname(dirs[0])
            self.model.wdirChanged.emit(self.model.w_dir_main)
            print(f"Set working directory to: {self.model.w_dir_main}")

    def open_drop_file(self, file):
        """Processes dropped files and updates the model."""
        if file.endswith(('.sdf', '.xyz')):
            self.model.files.append(file)
            self.view.display_selected_files([file])

    def load_results_from_source(self, base_dir, source="auto"):
        """Load compatible CMIN inputs discovered from AQME result folders."""
        discovered = discover_aqme_result_files(
            base_dir,
            source=source,
            extensions=(".sdf", ".xyz"),
            recursive=True,
        )
        if discovered:
            self.model.files = discovered
            self.model.filesChanged.emit(discovered)
            self.view.display_selected_files(discovered)
            self.model.w_dir_main = os.path.abspath(base_dir)
            self.model.wdirChanged.emit(self.model.w_dir_main)
        return discovered