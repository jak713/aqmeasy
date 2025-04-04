from PySide6.QtWidgets import QDialog, QTableWidget, QTableWidgetItem, QVBoxLayout
from PySide6.QtCore import Qt

class csv_table(QDialog):
    def __init__(self, csv_dictionary, smiles2pixmap, parent=None):
        super().__init__(parent)
        self.setWindowTitle("CSV Data")
        self.resize(1000, 500)
        self.setModal(False)
        
        layout = QVBoxLayout()
        table = QTableWidget()

        table.setRowCount(len(csv_dictionary["SMILES"]))
        table.setColumnCount(len(csv_dictionary.keys()))
        table.setHorizontalHeaderLabels(csv_dictionary.keys())
        table.verticalHeader().setDefaultSectionSize(120)
        table.horizontalHeader().setDefaultSectionSize(120)

        for row in range(len(self.csv_dictionary["SMILES"])):
            for col, key in enumerate(self.csv_dictionary.keys()):
                if key == "SMILES":
                    pixmap = smiles2pixmap(csv_dictionary[key][row])
                    item = QTableWidgetItem()
                    item.setData(Qt.DecorationRole, pixmap)
                    table.setItem(row, col, item)
                else:
                    item = QTableWidgetItem(str(self.csv_dictionary[key][row]))
                    table.setItem(row, col, item)
            
        for row in range(table.rowCount()):
            item = table.item(row, 0)
            if item:
                item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
        
        layout.addWidget(self.table_widget)
        self.table = table