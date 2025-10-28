from PySide6.QtWidgets import QDialog, QTableWidget, QTableWidgetItem, QVBoxLayout, QHBoxLayout, QPushButton, QGroupBox
from PySide6.QtCore import Qt

from aqmeasy.ui.stylesheets import stylesheets
from aqmeasy.utils import smiles2pixmap 

class csv_table(QDialog):
    def __init__(self, csv_model):
        super().__init__()
        self.model = csv_model
        self.setStyleSheet(stylesheets.QWidget)
        self.setWindowTitle("CSV Data")
        self.resize(1000, 500)

        main_layout = QVBoxLayout()
        group = QGroupBox("Double click to edit cells")
        layout = QVBoxLayout()

        table = QTableWidget()

        top_layout = QHBoxLayout()
        intermediate_button = QPushButton("Add Intermediate")
        ts_button = QPushButton("Add Transition State")

        top_layout.addWidget(intermediate_button)
        top_layout.addWidget(ts_button)

        self.table = table
        self.refresh_view()

        layout.addLayout(top_layout)
        layout.addWidget(table)
        group.setLayout(layout)
        main_layout.addWidget(group)
        self.setLayout(main_layout)
        self.table = table
        self.intermediate_button = intermediate_button
        self.ts_button = ts_button

    def refresh_view(self):
        """Check the model to refresh the view after changes"""
        table = self.table
        table.setRowCount(len(self.model["SMILES"]))
        table.setColumnCount(len(self.model.keys()))
        table.setHorizontalHeaderLabels(self.model.keys())
        table.verticalHeader().setDefaultSectionSize(120)
        table.horizontalHeader().setDefaultSectionSize(120)
        for row in range(len(self.model["SMILES"])):
            for col, key in enumerate(self.model.keys()):
                if key == "SMILES":
                    pixmap = smiles2pixmap(self.model[key][row])
                    item = QTableWidgetItem()
                    item.setData(Qt.DecorationRole, pixmap)
                    table.setItem(row, col, item)
                else:
                    item = QTableWidgetItem(str(self.model[key][row]))
                    table.setItem(row, col, item)

        for row in range(table.rowCount()):
            item = table.item(row, 0)
            if item:
                item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)