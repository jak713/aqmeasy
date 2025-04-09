import logging
import csv
from models.smiles2csv_model import csv_dictionary as csv_model
from models.smiles2csv_command import general_command_dictionary as gen_command
import ui.dialogs.smiles2csv_table_dialog as csv_table
from utils import smiles2pixmap, smiles2enumerate, smiles2charge, smiles2multiplicity
import rdkit.rdBase
from rdkit import Chem
from rdkit.Chem import Draw as rdMolDraw2D, rdDepictor
from PySide6.QtWidgets import QInputDialog, QMessageBox, QFileDialog
from PySide6.QtGui import QPixmap

class CsvController:
    """Controller for the CSV model"""
    def __init__(self, model):
        self.model = model
        self.current_index = 1
        self.total_index = self.get_total_index()

    def set_parent(self, parent):
        """Set the parent of the controller"""
        self.parent = parent

    def get_total_index(self):
        return len(self.model["SMILES"])
        
    def show_csv(self):
        csv = csv_table.csv_table(
            csv_model=self.model,
            smiles2pixmap=smiles2pixmap,
            parent=None 
        )
        csv.table.itemChanged.connect(lambda item: self.update_model_table(item))
        csv.intermediate_button.clicked.connect(lambda: csv.table.update if self.add_intermediate(csv.table.selectedItems()) else self.parent.failure("Please select two or more items to add an intermediate."))
        csv.ts_button.clicked.connect(lambda: csv.table.update if self.add_transition_state(csv.table.selectedItems()) else self.parent.failure("Please select one or more items to add a transition state."))
        csv.exec()

    def update_model_table(self, item):
        """Update the model when the signals are emitted"""
        row = item.row()
        col = item.column()
        key = list(self.model.keys())[col]
        self.model[key][row] = item.text()
        self.model.signals.updated.emit()

    def add_intermediate(self,items):
        """Add an intermediate to the CSV table by combining the SMILES strings of the selected items in the csv_table."""
        selected_items = [self.model["SMILES"][item.row()] for item in items]
        print(selected_items)
        if len(selected_items) < 2:
            return self.parent.failure("Please select two or more items to add an intermediate.")
        smiles = ".".join(selected_items)
        self.new_molecule()
        index = self.get_total_index() - 1
        self.model["SMILES"][index] = smiles
        self.model["code_name"][index] = "Intermediate"
        self.model["charge"][index] = smiles2charge(smiles)
        self.model["multiplicity"][index] = smiles2multiplicity(smiles)
        self.model["constraints_atoms"][index] = ""
        self.model["constraints_dist"][index] = ""
        self.model["constraints_angle"][index] = ""
        self.model["constraints_dihedral"][index] = ""
        self.model.signals.updated.emit()
        return True


    def add_transition_state(self, items):
        """Add a transition state to the CSV table."""
        selected_items = [self.model["SMILES"][item.row()] for item in items]
        print(selected_items)
        if len(selected_items) < 1:
            return self.parent.failure("Please select one or more items to add a transition state.")
        smiles = ".".join(selected_items)
        self.new_molecule()
        index = self.get_total_index() - 1
        self.model["SMILES"][index] = smiles2enumerate(smiles)
        self.model["code_name"][index] = "TS"
        self.model["charge"][index] = smiles2charge(smiles)
        self.model["multiplicity"][index] = smiles2multiplicity(smiles)
        self.model["constraints_atoms"][index] = ""
        self.model["constraints_dist"][index] = ""
        self.model["constraints_angle"][index] = ""
        self.model["constraints_dihedral"][index] = ""

        self.model.signals.updated.emit()
        return True

    def update_smiles_model(self, smiles):
        """Update the model with the current SMILES string.
        If constraints are present, update the model["SMILES"] with the enumerated SMILES string.
        If the SMILES string changes, remove the constraints."""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return
        if not smiles:
            self.model["SMILES"][self.current_index - 1] = ""
            self.model["code_name"][self.current_index - 1] = ""
            self.model["charge"][self.current_index - 1] = ""
            self.model["multiplicity"][self.current_index - 1] = ""
            self.model["constraints_atoms"][self.current_index - 1] = ""
            self.model["constraints_dist"][self.current_index - 1] = ""
            self.model["constraints_angle"][self.current_index - 1] = ""
            self.model["constraints_dihedral"][self.current_index - 1] = ""
            self.model["complex_type"][self.current_index - 1] = ""
            self.model["geom"][self.current_index - 1] = ""
            self.model.signals.updated.emit()
            if hasattr(self, 'each_dist') and self.current_index in self.each_dist:
                del self.each_dist[self.current_index]
            return
        enumerated_smiles = smiles2enumerate(smiles)
        if self.model["SMILES"][self.current_index - 1] != "" and self.model["SMILES"][self.current_index - 1] == enumerated_smiles:
            return
        if self.model["SMILES"][self.current_index - 1] != enumerated_smiles:
            self.model["constraints_atoms"][self.current_index - 1] = ""
            self.model["constraints_dist"][self.current_index - 1] = ""
            self.model["constraints_angle"][self.current_index - 1] = ""
            self.model["constraints_dihedral"][self.current_index - 1] = ""

        self.model["charge"][self.current_index - 1] = smiles2charge(smiles)
        self.model["multiplicity"][self.current_index - 1] = smiles2multiplicity(smiles)

        if self.model["constraints_dist"][self.current_index - 1] != "":
            self.model["SMILES"][self.current_index - 1] = enumerated_smiles

        elif self.model["constraints_angle"][self.current_index - 1] != "":
            self.model["SMILES"][self.current_index - 1] = enumerated_smiles

        elif self.model["constraints_dihedral"][self.current_index - 1] != "":
            self.model["SMILES"][self.current_index - 1] = enumerated_smiles

        else:
            self.model["SMILES"][self.current_index - 1] = smiles
        self.model.signals.updated.emit()
        
    def save_csv_file(self):
        """Save the csv_dictionary to a file."""
        try:
            for index, code_name in enumerate(csv_model["code_name"]):
                if code_name == "":
                    self.parent.failure("Please enter a code name for the molecule before saving.")
                    return False  # for the closing event

            file_name, _ = QFileDialog.getSaveFileName(self.parent, "Save CSV File", "", "CSV Files (*.csv)")
            if not file_name:
                self.parent.failure("Please enter a file name before saving.")
                return False  # for the closing event

            gen_command["input"] = file_name
            with open(gen_command["input"], 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(csv_model.keys())
                for i in range(len(csv_model["SMILES"])):
                    writer.writerow([csv_model[key][i] for key in csv_model.keys()])

            self.parent.success("CSV file saved successfully.")
            return True  # for the closing event

        except Exception as e:
            self.parent.failure(f"An error occurred while saving the file.")
            logging.error(f"Error saving CSV file: {str(e)}")

            return False  # for the closing event

    def display_molecule(self,checked = None):
        """Display the molecule in the molecule_label using rdkit.Chem.Draw module"""
        rdkit.rdBase.DisableLog('rdApp.*')
        rdDepictor.SetPreferCoordGen(True)

        smiles = csv_model["SMILES"][self.current_index - 1]
        if not smiles:
            self.parent.molecule_label.setText("No SMILES string provided.")
            self.atom_coords = None
            return

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid SMILES string.")
            
            mol = Chem.AddHs(mol)
            for i, atom in enumerate(mol.GetAtoms()):
                atom.SetAtomMapNum(i + 1)
            
            width = self.parent.molecule_label.width()
            height = self.parent.molecule_label.height()
            drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
            
            if checked:
                drawer.drawOptions().noAtomLabels = False
            else:
                drawer.drawOptions().noAtomLabels = True

            drawer.drawOptions().bondLineWidth = 1.5

            highlight_atoms = []
            highlight_colors = {}
            
            if hasattr(self, 'selected_atoms') and self.selected_atoms:
                highlight_atoms = [atom_idx - 1 for atom_idx in self.selected_atoms]
                highlight_colors = {idx: (0.5176, 0.6314, 0.9922) for idx in highlight_atoms}

            drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms, highlightAtomColors=highlight_colors)
            drawer.FinishDrawing()
            drawer.WriteDrawingText("/tmp/molecule.svg")
            
            self.atom_coords = [drawer.GetDrawCoords(i) for i in range(mol.GetNumAtoms())]

            pixmap = QPixmap("/tmp/molecule.svg")
            if self.parent.molecule_label is not None:
                self.parent.molecule_label.setPixmap(pixmap)
                
        except Exception as e:
            if self.parent.molecule_label is not None:
                self.parent.molecule_label.setText(f"Error displaying molecule: {str(e)}")

    def mousePressEvent(self, pos):
        """Handle mouse press events to select atoms and add constraints.
        The logic is to check if the mouse press event is within the molecule_label"""
        if self.parent.molecule_label and self.parent.molecule_label.geometry().contains(pos.toPoint()):
            x = pos.x() - self.parent.molecule_label.x()
            y = pos.y() - self.parent.molecule_label.y()
            selected_atom = self.get_atom_at_position(x, y)
            if selected_atom is not None:
                self.handle_atom_selection(selected_atom)
                self.display_molecule(self.parent.show_numbered_atoms_toggle.isChecked())  

    def get_atom_at_position(self, x, y):
        """Get the atom index at the given position by 
        checking the distance from the atom coordinates. 
        The atom coordinates are found using RDKit.
        The logic is to check if the distance between the mouse click
        and the atom coordinates is less than a threshold."""
        if not hasattr(self, 'atom_coords'):
            return None
        elif self.atom_coords is not None:
            for idx, coord in enumerate(self.atom_coords):
                if len(csv_model["SMILES"][self.current_index - 1]) <= 50: # small molecule = bigger click area
                    if (coord.x - x) ** 2 + (coord.y - y) ** 2 < 100: 
                        return idx + 1
                elif len(csv_model["SMILES"][self.current_index - 1]) > 50 : # big molecule = smaller click area BUT should probs keep playing around with these
                    if (coord.x - x) ** 2 + (coord.y - y) ** 2 < 40: 
                        return idx + 1
            return None

    def handle_atom_selection(self, atom_idx):
        """Handle the selection of an atom by adding constraints.
        If same atom pressed twice, deselect it.
        If two atoms are selected, add distance constraint.
        If three atoms are selected, add angle constraint.
        If four atoms are selected, add dihedral constraint."""
        if not hasattr(self, 'selected_atoms'):
            self.selected_atoms = []

        if atom_idx in self.selected_atoms:  # deselect if already selected
            self.selected_atoms.remove(atom_idx)
            self.display_molecule(self.parent.show_numbered_atoms_toggle.isChecked())  # update the highlighted atom
            return
        
        self.selected_atoms.append(atom_idx)
        self.display_molecule(self.parent.show_numbered_atoms_toggle.isChecked())  # update the highlighted atom
        atom1 = self.selected_atoms[0]
        mol = Chem.MolFromSmiles(csv_model["SMILES"][self.current_index - 1])
        mol = Chem.AddHs(mol) 
        atom1_type = mol.GetAtomWithIdx(atom1 - 1).GetSymbol()

        if len(self.selected_atoms) == 1:
            self.constraint_type, ok = QInputDialog.getItem(
                self.parent, 
                "Select Constraint Type", 
                f"Current selection: {atom1_type}{self.selected_atoms}.\nChoose constraint type:", 
                ["Distance", "Angle", "Dihedral"],
                0,
                False
            )
            if not ok:
                self.selected_atoms = []
                return
        
        if self.constraint_type == "Distance" and len(self.selected_atoms) == 2:
            self.add_distance_constraints()
        elif self.constraint_type == "Angle" and len(self.selected_atoms) == 3:
            self.add_angle_constraints()
        elif self.constraint_type == "Dihedral" and len(self.selected_atoms) == 4:
            self.add_dihedral_constraints()

    def add_distance_constraints(self):
        """Add distance constraints between two selected atoms."""
        if len(self.selected_atoms) == 2:
            atom1, atom2 = self.selected_atoms
            mol = Chem.MolFromSmiles(csv_model["SMILES"][self.current_index - 1])
            mol = Chem.AddHs(mol) 
            atom1_type = mol.GetAtomWithIdx(atom1 - 1).GetSymbol()
            atom2_type = mol.GetAtomWithIdx(atom2 - 1).GetSymbol()
            distance, ok = QInputDialog.getDouble(self.parent, "Distance Constraint", f"Enter distance between {atom1_type}:{atom1} and {atom2_type}:{atom2} (Å):")
            if ok:
                if distance <= 0:
                    QMessageBox.warning(self.parent, "Invalid Distance", "Distance must be greater than 0.", QMessageBox.StandardButton.Ok, QMessageBox.StandardButton.Ok)
                    return  self.add_distance_constraints()
                
                if not hasattr(self, 'each_dist'):
                    self.each_dist = {self.current_index: []}
                if self.current_index not in self.each_dist:
                    self.each_dist[self.current_index] = []
                self.each_dist[self.current_index].append([atom1, atom2, distance])
                csv_model["constraints_dist"][self.current_index - 1] = str(self.each_dist[self.current_index])
                self.parent.update_properties()
                self.selected_atoms = []
            else:
                self.selected_atoms = []
                return

    def add_angle_constraints(self):
        """Add angle constraints between three selected atoms."""
        if len(self.selected_atoms) == 3:
            atom1, atom2, atom3 = self.selected_atoms
            mol = Chem.MolFromSmiles(csv_model["SMILES"][self.current_index - 1])
            mol = Chem.AddHs(mol) 
            atom1_type = mol.GetAtomWithIdx(atom1 - 1).GetSymbol()
            atom2_type = mol.GetAtomWithIdx(atom2 - 1).GetSymbol()
            atom3_type = mol.GetAtomWithIdx(atom3 - 1).GetSymbol()
            angle, ok = QInputDialog.getDouble(self.parent, "Angle Constraint", f"Enter angle between {atom1_type}:{atom1}, {atom2_type}:{atom2}, and {atom3_type}:{atom3} (°):")
            if ok:
                if angle <= 0 or angle >= 360:
                    QMessageBox.warning(self.parent, "Invalid Angle", "Angle must be between 0 and 360 degrees.", QMessageBox.StandardButton.Ok, QMessageBox.StandardButton.Ok)
                    return self.addAngleConstraints()
                if not hasattr(self, 'each_angle'):
                    self.each_angle = {self.current_index: []}
                if self.current_index not in self.each_angle:
                    self.each_angle[self.current_index] = []
                self.each_angle[self.current_index].append([atom1, atom2, atom3, angle])
                csv_model["constraints_angle"][self.current_index - 1] = str(self.each_angle[self.current_index])
                self.parent.update_properties()
                self.selected_atoms = []
            else:
                self.selected_atoms = []
                return

    def add_dihedral_constraints(self):
        """Add dihedral constraints between four selected atoms."""
        if len(self.selected_atoms) == 4:
            atom1, atom2, atom3, atom4 = self.selected_atoms
            mol = Chem.MolFromSmiles(csv_model["SMILES"][self.current_index - 1])
            mol = Chem.AddHs(mol)
            atom1_type = mol.GetAtomWithIdx(atom1 - 1).GetSymbol()
            atom2_type = mol.GetAtomWithIdx(atom2 - 1).GetSymbol()
            atom3_type = mol.GetAtomWithIdx(atom3 - 1).GetSymbol()
            atom4_type = mol.GetAtomWithIdx(atom4 - 1).GetSymbol()
            dihedral, ok = QInputDialog.getDouble(self.parent, "Dihedral Constraint", f"Enter dihedral angle between {atom1_type}:{atom1}, {atom2_type}:{atom2}, {atom3_type}:{atom3}, and {atom4_type}:{atom4} (°):")
            if ok:
                if dihedral < -180 or dihedral > 180:
                    QMessageBox.warning(self.parent, "Invalid Dihedral", "Dihedral angle must be between -180 and 180 degrees.", QMessageBox.StandardButton.Ok, QMessageBox.StandardButton.Ok)
                    return self.addDihedralConstraints()
                if not hasattr(self, 'each_dihedral'):
                    self.each_dihedral = {self.current_index: []}
                if self.current_index not in self.each_dihedral:
                    self.each_dihedral[self.current_index] = []
                self.each_dihedral[self.current_index].append([atom1, atom2, atom3, atom4, dihedral])
                csv_model["constraints_dihedral"][self.current_index - 1] = str(self.each_dihedral[self.current_index])
                self.parent.update_properties()
                self.selected_atoms = []
            else:
                self.selected_atoms = []
                return
            
    def next_molecule(self):
        """Move to the next index in csv dictionary and update display."""
        if self.current_index == self.total_index == 1:
            return
        self.current_index = (self.current_index + 1) % (len(csv_model["SMILES"]) +1)
        if self.current_index == 0:
            self.current_index = 1
        self.parent.update_ui()

    def previous_molecule(self):
        """Move to the previous molecule and update display."""
        if self.current_index == self.total_index == 1:
            return
        self.current_index = (self.current_index - 1) % (len(csv_model["SMILES"])+1)
        if self.current_index == 0:
            self.current_index = len(csv_model["SMILES"])
        self.parent.update_ui()

    def new_molecule(self):
        """Create a new empty molecule entry in the csv dictionary and update the display."""
        for key in csv_model.keys():
            csv_model[key].append("")
        self.current_index = len(csv_model["SMILES"]) 
        self.total_index = self.get_total_index()
        csv_model.signals.updated.emit()
        self.parent.update_ui()

    def delete_molecule(self):
        """Delete the current molecule and all associated data (including constraints) from the csv dictionary and update the display."""
        if self.total_index == 1:
            for key, value in csv_model.items():
                value.pop(self.current_index - 1)
            self.new_molecule()
            return

        for key in csv_model.keys():
            csv_model[key].pop(self.current_index - 1)

        for constraint in ['each_dist', 'each_angle', 'each_dihedral']:
            if hasattr(self, constraint):
                constraint_dict = getattr(self, constraint)
                if self.current_index in constraint_dict:
                    del constraint_dict[self.current_index]

        if hasattr(self, 'selected_atoms'):
            self.selected_atoms = []

        self.total_index = self.get_total_index()
        self.current_index = (self.current_index - 1) % (self.total_index + 1)
        if self.current_index == 0:
            self.current_index = 1
        self.parent.update_ui()

csv_controller = CsvController(csv_model)