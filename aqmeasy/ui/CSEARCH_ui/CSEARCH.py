from PySide6.QtWidgets import QWidget, QVBoxLayout
from aqmeasy.ui.stylesheets import stylesheets

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
        """Open QPREP with discovered CSEARCH SDF files."""
        qprep_widget = self.parent.new_qprep_widget()  # type: ignore[attr-defined]
        files = discover_aqme_result_files(
            destination_folder,
            source="csearch",
            extensions=(".sdf",),
            recursive=True,
        )
        qprep_widget.file_panel.get_files_from_csearch(files)

    def open_cmin_after_csearch(self, destination_folder: str):
        """Open CMIN with files discovered from CSEARCH output directories."""
        cmin_widget = self.parent.new_cmin_widget()  # type: ignore[attr-defined]
        cmin_widget.file_panel.controller.load_results_from_source(destination_folder, source="csearch")

    def open_qdescp_after_csearch(self, destination_folder: str):
        """Open QDESCP with SDF files discovered from CSEARCH output directories."""
        qdescp_widget = self.parent.new_qdescp_widget()  # type: ignore[attr-defined]
        files = discover_aqme_result_files(destination_folder, source="csearch", extensions=(".sdf",), recursive=True)
        qdescp_widget.set_input_files(files)
