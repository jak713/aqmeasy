from models.smiles2csv_model import csv_dictionary as csv_model
import ui.dialogs.smiles2csv_table_dialog as csv_table
from utils import smiles2pixmap, smiles2enumerate, smiles2charge, smiles2multiplicity
from rdkit import Chem

class CsvController:
    """Controller for the CSV model"""
    def __init__(self, model):
        self.model = model
        
    def show_csv(self):
        csv = csv_table.csv_table(
            csv_model=self.model,
            smiles2pixmap=smiles2pixmap,
            parent=None
        )
        csv.table.itemChanged.connect(lambda item: self.update_model_table(item))
        # csv.intermediate_button.clicked.connect(
        #     lambda: self.add_intermediate(csv.table.currentItem())
        # )
        csv.exec()

    def update_model_table(self, item):
        """Update the model when the signals are emitted"""
        row = item.row()
        col = item.column()
        key = list(self.model.keys())[col]
        self.model[key][row] = item.text()
        self.model.signals.updated.emit()

    # def add_intermediate(self):
    #     """Add an intermediate to the CSV table by combining the SMILES strings of the selected items in the csv_table."""
    #     selected_items = self.table.selectedItems()
    #     if len(selected_items) < 2:
    #         return
    #     smiles = ""
    #     for item in selected_items:
    #         smiles += self.model["SMILES"][item.row()] + "."
    #     smiles = smiles[:-1]
    #     index = self.table.currentRow()
    #     self.model["SMILES"][index] = smiles
    #     self.model["code_name"][index] = "Intermediate"
    #     self.model["charge"][index] = smiles2charge(smiles)
    #     self.model["multiplicity"][index] = smiles2multiplicity(smiles)
    #     self.model["constraints_atoms"][index] = ""
    #     self.model["constraints_dist"][index] = ""
    #     self.model["constraints_angle"][index] = ""
    #     self.model["constraints_dihedral"][index] = ""

    def add_transition_state(self):
        """Add a transition state to the CSV table."""
        pass

    def update_smiles_model(self, smiles, index):
        """Update the model with the current SMILES string.
        If constraints are present, update the model["SMILES"] with the enumerated SMILES string.
        If the SMILES string changes, remove the constraints."""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return
        if not smiles:
            self.model["SMILES"][index - 1] = ""
            self.model["code_name"][index - 1] = ""
            self.model["charge"][index - 1] = ""
            self.model["multiplicity"][index - 1] = ""
            self.model["constraints_atoms"][index - 1] = ""
            self.model["constraints_dist"][index - 1] = ""
            self.model["constraints_angle"][index - 1] = ""
            self.model["constraints_dihedral"][index - 1] = ""
            self.model["complex_type"][index - 1] = ""
            self.model["geom"][index - 1] = ""
            self.model.signals.updated.emit()
            return
        enumerated_smiles = smiles2enumerate(smiles)
        if self.model["SMILES"][index - 1] != "" and self.model["SMILES"][index - 1] == enumerated_smiles:
            return
        if self.model["SMILES"][index - 1] != enumerated_smiles:
            self.model["constraints_atoms"][index - 1] = ""
            self.model["constraints_dist"][index - 1] = ""
            self.model["constraints_angle"][index - 1] = ""
            self.model["constraints_dihedral"][index - 1] = ""

        self.model["charge"][index - 1] = smiles2charge(smiles)
        self.model["multiplicity"][index - 1] = smiles2multiplicity(smiles)

        if self.model["constraints_dist"][index - 1] != "":
            self.model["SMILES"][index - 1] = enumerated_smiles

        elif self.model["constraints_angle"][index - 1] != "":
            self.model["SMILES"][index - 1] = enumerated_smiles

        elif self.model["constraints_dihedral"][index - 1] != "":
            self.model["SMILES"][index - 1] = enumerated_smiles

        else:
            self.model["SMILES"][index - 1] = smiles
        self.model.signals.updated.emit()
        
    def save_csv_file(self):
        """Save the csv_dictionary to a file."""
        pass


csv_controller = CsvController(csv_model)
