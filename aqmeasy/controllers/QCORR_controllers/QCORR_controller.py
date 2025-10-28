from PySide6.QtWidgets import QFileDialog
from PySide6.QtCore import Signal, QObject, Slot, QRunnable
from aqme.qcorr import qcorr

import os

FILE_FILTERS = [
    "Gaussian Output Files (*.log)",
    "ORCA Output Files (*.out)",
]

class FileController(QObject):
    """Controller for the QCORR module, connecting the model and the view."""
    def __init__(self, model, view):
        super().__init__()
        self.model = model
        self.view = view

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
            self.view.output_dir_label.setText(directory)
        print(self.model.__get__w_dir_main__())

    def open_drop_folder(self, dirs):
        """Processes dropped folders to find relevant files and updates the model."""
        all_files = []
        for dir_path in dirs:
            for root, _, files in os.walk(dir_path):
                for file in files:
                    if file.endswith(('.log', '.out')):
                        all_files.append(os.path.join(root, file))
        if all_files:
            self.model.files = all_files
            self.model.filesChanged.emit(all_files)
            self.view._display_selected_files(all_files)

    def open_drop_file(self, file):
        """Processes dropped files and updates the model."""
        if file.endswith(('.log', '.out')):
            self.model.files.append(file)
            self.view._display_selected_files([file])