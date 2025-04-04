from PySide6.QtWidgets import QDialog, QTableWidget, QTableWidgetItem, QVBoxLayout
from PySide6.QtGui import QPixmap
from PySide6.QtCore import Qt

import models.smiles2csv_model as model
import ui.dialogs.smiles2csv_table_dialog as csv_table
from utils import smiles2pixmap 

class csv_controller:
    """Controller for the csv model"""
    def __init__(self, model):
        self.model = model
        self.model.signals.updated.connect(self.update_model)

    def update_model(self):
        """Update the model when the signals are emitted"""
        # Add any additional logic needed when the model is updated
        pass

    def show_csv(self):
        csv = csv_table.csv_table(
            csv_dictionary=model.csv_dictionary,
            smiles2pixmap=smiles2pixmap,
            parent=self.smiles_to_csv)
        csv.table.itemChanged.connect(self.update_csv_dictionary)
        csv.show()
        
