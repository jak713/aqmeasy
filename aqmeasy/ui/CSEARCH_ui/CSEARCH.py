from PySide6.QtWidgets import QWidget, QVBoxLayout
from PySide6.QtCore import Qt
from aqmeasy.ui.stylesheets import stylesheets

from aqmeasy.ui.CSEARCH_ui.CSEARCH_widget import CSEARCHWidget
from aqmeasy.models.CSEARCH_model.CSEARCH_model import csv_dictionary

class CSEARCH(QWidget):
    def __init__(self):
        super().__init__()
        self.model = csv_dictionary
        self.main_widget = CSEARCHWidget(csv_model=self.model)
        self.resize(1000, 900)

        self.setStyleSheet(stylesheets.QWidget)
        self.setWindowTitle("CSEARCH")
        layout = QVBoxLayout()
        layout.addWidget(self.main_widget)
        self.setLayout(layout)