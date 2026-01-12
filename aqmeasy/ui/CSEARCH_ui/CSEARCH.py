from PySide6.QtWidgets import QWidget, QVBoxLayout
from PySide6.QtCore import Qt
from aqmeasy.ui.stylesheets import stylesheets

from aqmeasy.ui.CSEARCH_ui.CSEARCH_widget import CSEARCHWidget
from aqmeasy.controllers.CSEARCH_controller import CSEARCHThread
from aqmeasy.models.CSEARCH_model.CSEARCH_model import csv_dictionary
from aqmeasy.models.CSEARCH_model.CSEARCH_command import general_command_model

class CSEARCH(QWidget):
    def __init__(self, parent=None):
        super().__init__()
        if parent:
            self.parent = parent
        self.model = csv_dictionary
        self.worker = CSEARCHThread(self, general_command_model)
        self.main_widget = CSEARCHWidget(self, self.model, general_command_model)
        self.resize(1000, 900)

        self.setStyleSheet(stylesheets.QWidget)
        self.setWindowTitle("CSEARCH")
        layout = QVBoxLayout()
        layout.addWidget(self.main_widget)
        self.setLayout(layout)

    def open_qprep_after_csearch(self, destination_folder: str):
        """Open QPREP from parent (main_window) with the generated SDF files after CSEARCH run."""
        QPREP = self.parent.new_qprep_widget() # type: ignore # 
        QPREP.file_panel.get_files_from_csearch([f"{destination_folder}/{name}_{general_command_model['program']}.sdf" for name in self.model["code_name"] if name])