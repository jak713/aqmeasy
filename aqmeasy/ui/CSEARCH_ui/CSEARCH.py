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

    def _discover_csearch_sdf_paths(self, destination_folder: str):
        discovered = discover_aqme_result_files(
            destination_folder,
            source="csearch",
            extensions=(".sdf",),
            recursive=True,
        )
        if discovered:
            return discovered

        # Fallback to deterministic naming used by CSEARCH outputs.
        program = str(general_command_model.get("program", "") or "").strip()
        fallback = []
        if program:
            for name in self.model.get("code_name", []):
                if not name:
                    continue
                sdf_path = f"{destination_folder}/{name}_{program}.sdf"
                if os.path.exists(sdf_path):
                    fallback.append(sdf_path)
        return fallback

    def open_cmin_after_csearch(self, destination_folder: str):
        """Open CMIN with generated CSEARCH structures as input."""
        parent_window = getattr(self, "parent", None)
        if parent_window is None or not hasattr(parent_window, "new_cmin_widget"):
            return
        cmin_widget = parent_window.new_cmin_widget()  # type: ignore
        cmin_widget.file_panel.controller.load_results_from_source(destination_folder, source="csearch")

    def open_qprep_after_csearch(self, destination_folder: str):
        """Open QPREP from parent (main_window) with generated CSEARCH SDF files."""
        parent_window = getattr(self, "parent", None)
        if parent_window is None or not hasattr(parent_window, "new_qprep_widget"):
            return
        qprep_widget = parent_window.new_qprep_widget()  # type: ignore
        qprep_widget.file_panel.get_files_from_csearch(self._discover_csearch_sdf_paths(destination_folder))

    def open_qdescp_after_csearch(self, destination_folder: str):
        """Open QDESCP from parent (main_window) with generated CSEARCH SDF files."""
        parent_window = getattr(self, "parent", None)
        if parent_window is None or not hasattr(parent_window, "new_qdescp_widget"):
            return
        qdescp_widget = parent_window.new_qdescp_widget()  # type: ignore
        sdf_paths = self._discover_csearch_sdf_paths(destination_folder)
        payload = {"files": sdf_paths}
        payload.update(self._extract_qdescp_prefill_from_sdf(sdf_paths))
        qdescp_widget.set_input_payload(payload)

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

    