from PySide6.QtWidgets import QDialog, QTableWidget, QTableWidgetItem, QVBoxLayout
from PySide6.QtGui import QPixmap
from PySide6.QtCore import Qt

from models.smiles2csv_model import csv_dictionary as csv_model
import ui.dialogs.smiles2csv_table_dialog as csv_table
from utils import smiles2pixmap, smiles2enumerate

class csv_controller:
    """Controller for the csv model"""
    def __init__(self, model):
        self.model = model
        self.model.signals.updated.connect(self.update_model)

    def update_model(self):
        """Update the model when the signals are emitted"""
    
        pass

    def show_csv(self):
        csv = csv_table.csv_table(
            csv_model=csv_model,
            smiles2pixmap=smiles2pixmap,
            parent=None
        )
        # csv.table.itemChanged.connect(self.update_model_from_table)
        # csv.intermediate_button.clicked.connect(self.add_intermediate)
        # csv.ts_button.clicked.connect(self.add_transition_state)
        csv.exec()
        
    def add_intermediate(self):
        """Add an intermediate to the csv table."""
        pass

    def add_transition_state(self):
        """Add a transition state to the csv table."""
        pass

    def update_model_from_table(self):
        """Update the model if changes within csv_table occur."""
        pass

    def update_smiles_csv(self, smiles, index):
        """Update the csv_dictionary with the current SMILES string.
        If constraints are present, update the csv_dictionary["SMILES"] with the enumerated SMILES string."""
        if not smiles:
            return
        enumerated_smiles = smiles2enumerate(smiles) 
        if csv_model["SMILES"][index - 1] != "" and csv_model["SMILES"][index- 1] == enumerated_smiles:
            return

        if csv_model["constraints_dist"][index - 1] != "":
            csv_model["SMILES"][index - 1] = enumerated_smiles

        elif csv_model["constraints_angle"][index - 1] != "":
            csv_model["SMILES"][index - 1] = enumerated_smiles

        elif csv_model["constraints_dihedral"][index - 1] != "":
            csv_model["SMILES"][index - 1] = enumerated_smiles

        else:
            csv_model["SMILES"][index - 1] = smiles
        
    def save_csv_file(self):
        """Save the csv_dictionary to a file."""
        pass