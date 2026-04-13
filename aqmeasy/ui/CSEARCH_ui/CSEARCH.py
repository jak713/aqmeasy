import os

from PySide6.QtGui import QCloseEvent
from PySide6.QtWidgets import QWidget, QVBoxLayout
from aqmeasy.ui.stylesheets import stylesheets
from aqmeasy.models.QDESCP_model.aqmetab_model import extract_qdescp_prefill_from_sdf

from aqmeasy.ui.CSEARCH_ui.CSEARCH_widget import CSEARCHWidget
from aqmeasy.controllers.CSEARCH_controller import CSEARCHThread
from aqmeasy.models.CSEARCH_model.CSEARCH_model import csv_dictionary
from aqmeasy.models.CSEARCH_model.CSEARCH_command import general_command_model
from aqmeasy.utils import discover_aqme_result_files

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

    @staticmethod
    def _extract_qdescp_prefill_from_sdf(file_paths):
        return extract_qdescp_prefill_from_sdf(file_paths)

    def get_generated_sdf_paths(self):
        destination = str(general_command_model.get("destination", "") or "").strip()
        program = str(general_command_model.get("program", "") or "").strip()
        if not destination or not program:
            return []

        paths = []
        for name in self.model.get("code_name", []):
            if not name:
                continue
            sdf_path = f"{destination}/{name}_{program}.sdf"
            if os.path.exists(sdf_path):
                paths.append(sdf_path)

        deduped = []
        seen = set()
        for path in paths:
            if path in seen:
                continue
            seen.add(path)
            deduped.append(path)
        return deduped

    