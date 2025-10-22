from PySide6.QtCore import Signal, QObject, Slot
from aqme.qcorr import qcorr
import os
import shutil

class qcorrController(QObject):
    """Controller for the QCORR module, connecting the model and the view."""
    def __init__(self, model, view):
        super().__init__()
        self.model = model
        self.view = view

        # Connect view signals to controller slots
        self.view.file_panel.get_filenames_signal.connect(self.update_files)
        self.view.file_panel.get_directory_signal.connect(self.update_wdir)

        # Connect model signals to view update methods
        self.model.filesChanged.connect(self.view.update_file_list)
        self.model.wdirChanged.connect(self.view.update_working_directory)