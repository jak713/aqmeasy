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
        self.resize(1280, 500)

        main_layout = QVBoxLayout()
        group = QGroupBox("Double click to edit cells. Hold ⌘ or Ctrl to select multiple cells.")
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

    def refresh_view(self) -> None:
        """Check the model to refresh the view after changes"""
        table = self.table
        # There was a feedback loop when add_intermediate was triggered so the signals are blocked for the refreshing and then unblocked. Seems to work as intended now.
        table.blockSignals(True)
        table.setRowCount(len(self.model.__getitem__("SMILES")))
        table.setColumnCount(self.model.__len__())
        table.setHorizontalHeaderLabels(self.model.keys())
        table.verticalHeader().setDefaultSectionSize(120)
        table.horizontalHeader().setDefaultSectionSize(120)

        for row in range(len(self.model.__getitem__("SMILES"))):
            for col, key in enumerate(self.model.keys()):
                if key == "SMILES":
                    pixmap = smiles2pixmap(self.model.__getitem__(key)[row])
                    item = QTableWidgetItem()
                    item.setData(Qt.ItemDataRole.DecorationRole, pixmap) # Data in the form of an icon
                    table.setItem(row, col, item)
                else:
                    item = QTableWidgetItem(str(self.model.__getitem__(key)[row]))
                    table.setItem(row, col, item)

        for row in range(table.rowCount()):
            item = table.item(row, 0)
            if item:
                item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
        table.blockSignals(False)

    def get_row_count(self) -> int:
        return self.table.rowCount()
    
    def get_column_count(self) -> int:
        return self.table.columnCount()
    
    def get_item(self, row: int, column: int):
        return self.table.item(row, column)