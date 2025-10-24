from PySide6.QtWidgets import QDialog, QTableWidget, QTableWidgetItem, QVBoxLayout, QHBoxLayout, QPushButton, QGroupBox
from PySide6.QtCore import Qt
from aqmeasy.ui.stylesheets import stylesheets

class csv_table(QDialog):
    def __init__(self, csv_model, smiles2pixmap, parent=None):
        super().__init__(parent)
        self.setStyleSheet(stylesheets.QDialog)
        self.setWindowTitle("CSV Data")
        self.resize(1000, 500)
        self.setModal(False)

        main_layout = QVBoxLayout()
        group = QGroupBox("Double click to edit cells")
        group.setStyleSheet(stylesheets.QGroupBox)
        layout = QVBoxLayout()

        table = QTableWidget()
        table.setStyleSheet(stylesheets.QTableWidget)

        top_layout = QHBoxLayout()
        intermediate_button = QPushButton("Add Intermediate")
        intermediate_button.setStyleSheet(stylesheets.QPushButton)
        ts_button = QPushButton("Add Transition State")
        ts_button.setStyleSheet(stylesheets.QPushButton)

        top_layout.addWidget(intermediate_button)
        top_layout.addWidget(ts_button)

        table.setRowCount(len(csv_model["SMILES"]))
        table.setColumnCount(len(csv_model.keys()))
        table.setHorizontalHeaderLabels(csv_model.keys())
        table.verticalHeader().setDefaultSectionSize(120)
        table.horizontalHeader().setDefaultSectionSize(120)

        for row in range(len(csv_model["SMILES"])):
            for col, key in enumerate(csv_model.keys()):
                if key == "SMILES":
                    pixmap = smiles2pixmap(csv_model[key][row])
                    item = QTableWidgetItem()
                    item.setData(Qt.DecorationRole, pixmap)
                    table.setItem(row, col, item)
                else:
                    item = QTableWidgetItem(str(csv_model[key][row]))
                    table.setItem(row, col, item)
            
        for row in range(table.rowCount()):
            item = table.item(row, 0)
            if item:
                item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
        
        layout.addLayout(top_layout)
        layout.addWidget(table)
        group.setLayout(layout)
        main_layout.addWidget(group)
        self.setLayout(main_layout)
        self.table = table
        self.intermediate_button = intermediate_button
        self.ts_button = ts_button