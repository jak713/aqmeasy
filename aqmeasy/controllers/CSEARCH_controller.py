import logging
import csv
import os, sys
from io import StringIO
import psutil
import time
import threading
import select
import subprocess
import tempfile
from aqmeasy.models.CSEARCH_model.CSEARCH_command import general_command_default, crest_command_model as crest_command
from aqmeasy.ui.CSEARCH_ui.CSEARCH_csvtable import csv_table
from aqmeasy.ui.CSEARCH_ui.CSEARCH_SmilesTutorialViewer import SmilesTutorialViewer
from aqmeasy.utils import smiles2enumerate, smiles2charge, smiles2multiplicity

from aqme.csearch import csearch

import rdkit.rdBase
from rdkit import Chem
from rdkit.Chem import Draw as rdMolDraw2D, rdDepictor
from PySide6.QtWidgets import QInputDialog, QMessageBox, QFileDialog
from PySide6.QtGui import QPixmap
from PySide6.QtCore import QObject, Signal, QThread, QTimer

ELECTRON_LABEL_X_OFFSET = 15
ELECTRON_LABEL_Y_OFFSET = 25
CLICK_RADIUS_SMALL_MOLECULE = 100
CLICK_RADIUS_LARGE_MOLECULE = 40
LARGE_MOLECULE_SMILES_LENGTH_THRESHOLD = 50

class CsvController(QObject):
    """Controller for the CSV model"""
    def __init__(self, model, gen_command_model):
        self.model = model
        self.gen_command_model = gen_command_model
        self.current_index = 1
        self.total_index = self.model.get_total_index()
        self.model.signals.updated.connect(self.update_table_from_model)

    def set_parent(self, parent) -> None:
        """Set the parent of the controller"""
        self.parent = parent

    def show_csv(self) -> None:
        self.csv = csv_table(
            csv_model=self.model,
        )
        self.csv.table.itemChanged.connect(lambda item: self.update_model_from_table(item))

        self.csv.intermediate_button.clicked.connect(lambda: self.csv.refresh_view() if self.add_intermediate(self.csv.table.selectedItems()) else ...)

        self.csv.ts_button.clicked.connect(lambda: self.csv.refresh_view() if self.add_transition_state(self.csv.table.selectedItems()) else self.parent.failure("Please select one or more items to add a transition state. Select multiple SMILES with Cmd/Ctrl+right click"))
        self.csv.show()

    def update_model_from_table(self, item) -> None:
        """Update the model when signal emitted from the csv table.
            Args: item (?)
            Returns None"""
        row = item.row()
        col = item.column()
        key = list(self.model.keys())[col]
        self.model.__getitem__(key)[row] = item.text()
        self.model.signals.updated.emit()

    def update_table_from_model(self) -> None:
        # Disconnect signals to prevent recursive updates
        if hasattr(self, 'csv'):
            self.csv.table.blockSignals(True)
            self.csv.refresh_view()
            self.csv.table.blockSignals(False)
        else:
            return

    def add_intermediate(self,items: list) -> bool:
        """Add an intermediate to the CSV table by combining the SMILES strings of the selected items in the csv_table.
        
            Args: items selected (list)
            Returns False if items selected are fewer than two or if changing the model fails, True otherwise"""
        selected_items = [self.model["SMILES"][item.row()] for item in items]
        selected_code_names = [self.model["code_name"][item.row()] for item in items]

        if len(selected_items) < 2:
            self.parent.failure("Please select two or more items to add an intermediate. Select multiple SMILES with Cmd/Ctrl+right click.")
            return False
        try:
            smiles = ".".join(selected_items)
            self.new_molecule()
            index = self.model.get_total_index() - 1
            self.model["SMILES"][index] = smiles
            self.model["code_name"][index] = "Intermediate_" + "+".join(selected_code_names)
            self.model["charge"][index] = smiles2charge(smiles)
            self.model["multiplicity"][index] = smiles2multiplicity(smiles)
            self.model["constraints_atoms"][index] = ""
            self.model["constraints_dist"][index] = ""
            self.model["constraints_angle"][index] = ""
            self.model["constraints_dihedral"][index] = ""
            self.model.signals.updated.emit()
            return True
        except Exception as e:
            logging.warning(f"Error encountered while adding intermediate: {e}")
            self.parent.failure("Something went wrong when trying to create an intermediate. Please try again.")
            return False

    def add_transition_state(self, items):
        """Add a transition state to the CSV table."""
        selected_items = [self.model["SMILES"][item.row()] for item in items]
        selected_code_names = [self.model["code_name"][item.row()] for item in items]
        if len(selected_items) < 1:
            return False
        
        smiles = ".".join(selected_items)
        self.new_molecule()
        index = self.model.get_total_index() - 1

        # to avoid duplicate names, need to check if the name already exists, if it does do magic trick of adding number at the end
        base_TS_name = "TS_" + "+".join(selected_code_names) 
        existing_names = self.model["code_name"]
        if base_TS_name in existing_names:
            count = sum(1 for name in existing_names if name.startswith(base_TS_name))
            unique_TS_name = f"{base_TS_name}_{count}"
        else:
            unique_TS_name = base_TS_name

        self.model["SMILES"][index] = smiles2enumerate(smiles)
        self.model["code_name"][index] = unique_TS_name
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
            self.parent.smiles_are_bad_bro(smiles)
            return
        if not smiles:
            self.model.set_item_at_index("SMILES", self.current_index - 1, "")
            self.model.set_item_at_index("code_name", self.current_index - 1, "")
            self.model.set_item_at_index("charge", self.current_index - 1, "")
            self.model.set_item_at_index("multiplicity", self.current_index - 1, "")
            self.model.set_item_at_index("constraints_atoms", self.current_index - 1, "")
            self.model.set_item_at_index("constraints_dist", self.current_index - 1, "")
            self.model.set_item_at_index("constraints_angle", self.current_index - 1, "")
            self.model.set_item_at_index("constraints_dihedral", self.current_index - 1, "")
            self.model.set_item_at_index("complex_type", self.current_index - 1, "")
            self.model.set_item_at_index("geom", self.current_index - 1, "")
            if self.model.get_total_index() == 1:
                self.gen_command_model["input"] = ""
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
    
        self.gen_command_model["input"] = None
        self.model.signals.updated.emit()
        
    def update_command(self, key: str, value) -> None:
        """Update the command model with the given key and value."""
        self.gen_command_model.__setitem__(key, value)

    def update_crest_command(self, key: str, value) -> None:
        """Update the crest command model with the given key and value."""
        crest_command.__setitem__(key, value)

    def import_file(self,file_name=None):
        """Import an SDF or ChemDraw file, extract SMILES, and display them. For CSV files, read the data and update the model."""

        if not file_name:
            file_name, _ = QFileDialog.getOpenFileName(self.parent, "Import File", "", "Accepted Files (*.cdxml *.sdf *.csv);;ChemDraw Files (*.cdx *.cdxml);;SDF Files (*.sdf);;CSV files (*.csv)")
        try:
            for key in self.model.keys():
                self.model[key].clear()
            smiles_list = []

            if file_name.endswith(".sdf"):
                mol_supplier = Chem.SDMolSupplier(file_name)
                for mol in mol_supplier:
                    if mol is not None:
                        smiles = Chem.MolToSmiles(mol)
                        smiles_list.append(smiles)

            elif file_name.endswith(".cdx"): 
                QMessageBox.warning(self.parent, "CDX File Import Warning","ChemDraw .cdx files are converted via Open Babel to SDF and then read. This may cause issues, and is generally not recommended. Consider using .cdxml files instead.")
                with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as temp_sdf:
                    temp_sdf_path = temp_sdf.name 

                try:
                    # Run Open Babel to convert CDX to SDF and save it to the temp file
                    subprocess.run(["obabel", file_name, "-O", temp_sdf_path], check=True)
                    # print(f"Conversion successful: {temp_sdf_path}")

                except subprocess.CalledProcessError as e:
                    print(f"Error during conversion: {e}")

                mol_supplier = Chem.SDMolSupplier(temp_sdf_path)
                for mol in mol_supplier:
                    if mol is not None:
                        smiles = Chem.MolToSmiles(mol)
                        smiles_list.append(smiles)

                if os.path.exists(temp_sdf_path):
                    os.remove(temp_sdf_path)

            elif file_name.endswith(".cdxml"):
                try:
                    mols = Chem.MolsFromCDXMLFile(file_name, sanitize=True, removeHs=True)
                    for mol in mols:
                        if mol is not None:
                            smiles = Chem.MolToSmiles(mol)
                            smiles_list.append(smiles)
                except Exception as e:
                    QMessageBox.critical(self, "CDXML Read Error", f"Failed to read {file_name}:\n{str(e)}")

            elif file_name.endswith(".csv"):
                with open(file_name, 'r') as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader: 
                        self.model.add_row(row)
                    
                self.total_index = self.model.get_total_index()
                if self.total_index > 0:
                    self.current_index = 1  # assuming indices start from 1
                    self.parent.update_properties()
                    self.parent.update_ui()
                
                self.gen_command_model.__setitem__("destination", os.path.dirname(file_name))
                return

            # If user cancels
            elif not file_name:
                self.model.add_row({}) # Make sure we are not left with zero rows
                return

            else:
                QMessageBox.warning(self.parent, "Error", "Unsupported file format. Please select a ChemDraw, SDF, or CSV file.")
                return
            if not smiles_list:
                QMessageBox.warning(self.parent, "Error", "No valid SMILES found in the file.")
                return

            for index, smiles in enumerate(smiles_list):
                if index < len(self.model["SMILES"]):
                    self.model["SMILES"][index] = smiles
                    self.model["charge"][index] = smiles2charge(smiles)
                    self.model["multiplicity"][index] = smiles2multiplicity(smiles)
                else:
                    for key in self.model.keys():
                        self.model[key].append("")
                    self.model["SMILES"][-1] = smiles
                    self.model["charge"][-1] = smiles2charge(smiles)
                    self.model["multiplicity"][-1] = smiles2multiplicity(smiles)
            self.total_index = self.model.get_total_index()

            for index in range(self.total_index):
                self.model["code_name"][index] = f"mol_{index + 1}"
                self.current_index = 1
                
            self.parent.update_ui()
            self.parent.update_properties()

        except ImportError as e:
            print(f"Error importing required module: {e}")
        except IOError as e:
            print(f"Error reading or writing file: {e}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

    def save_csv_file(self):
        """Save the csv_dictionary to a file."""
        try:
            for index, code_name in enumerate(self.model["code_name"]):
                if code_name == "":
                    self.parent.failure("Please enter a code name for the molecule before saving.")
                    return False  # for the closing event

            file_name, _ = QFileDialog.getSaveFileName(self.parent, "Save CSV File", "", "CSV Files (*.csv)")
            if not file_name:
                self.parent.failure("Please enter a file name before saving.")
                return False  # for the closing event

            self.gen_command_model.__setitem__("input", file_name)
            with open(self.gen_command_model.__getitem__("input"), 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(self.model.keys())
                for i in range(len(self.model["SMILES"])):
                    writer.writerow([self.model[key][i] for key in self.model.keys()])

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
        try:
            smiles = self.model["SMILES"][self.current_index - 1]
        except IndexError:
            smiles = ""
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
            
            width = self.parent.molecule_label.width() # this results in a bug when window is reopened atm
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
                highlight_colors = {idx: (0.565, 0.878, 0.937) for idx in highlight_atoms} #normalised RGB

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
            x = pos.x() - self.parent.molecule_label.x() + ELECTRON_LABEL_X_OFFSET
            y = pos.y() - self.parent.molecule_label.y() + ELECTRON_LABEL_Y_OFFSET
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
            print("Atom coordinates not available.")
            return None
        elif self.atom_coords is not None:
            for idx, coord in enumerate(self.atom_coords):
                if len(self.model["SMILES"][self.current_index - 1]) <= LARGE_MOLECULE_SMILES_LENGTH_THRESHOLD: 
                    if (coord.x - x) ** 2 + (coord.y - y) ** 2 < CLICK_RADIUS_SMALL_MOLECULE: 
                        return idx + 1
                elif len(self.model["SMILES"][self.current_index - 1]) > LARGE_MOLECULE_SMILES_LENGTH_THRESHOLD: 
                    if (coord.x - x) ** 2 + (coord.y - y) ** 2 < CLICK_RADIUS_LARGE_MOLECULE: 
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
        mol = Chem.MolFromSmiles(self.model["SMILES"][self.current_index - 1])
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
            mol = Chem.MolFromSmiles(self.model["SMILES"][self.current_index - 1])
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
                self.model["constraints_dist"][self.current_index - 1] = str(self.each_dist[self.current_index])
                self.parent.update_properties()
                self.selected_atoms = []
            else:
                self.selected_atoms = []
                return

    def add_angle_constraints(self):
        """Add angle constraints between three selected atoms."""
        if len(self.selected_atoms) == 3:
            atom1, atom2, atom3 = self.selected_atoms
            mol = Chem.MolFromSmiles(self.model["SMILES"][self.current_index - 1])
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
                self.model["constraints_angle"][self.current_index - 1] = str(self.each_angle[self.current_index])
                self.parent.update_properties()
                self.selected_atoms = []
            else:
                self.selected_atoms = []
                return

    def add_dihedral_constraints(self):
        """Add dihedral constraints between four selected atoms."""
        if len(self.selected_atoms) == 4:
            atom1, atom2, atom3, atom4 = self.selected_atoms
            mol = Chem.MolFromSmiles(self.model["SMILES"][self.current_index - 1])
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
                self.model["constraints_dihedral"][self.current_index - 1] = str(self.each_dihedral[self.current_index])
                self.parent.update_properties()
                self.selected_atoms = []
            else:
                self.selected_atoms = []
                return

    def next_molecule(self):
        """Move to the next index in csv dictionary and update display."""
        if self.current_index == self.total_index == 1:
            return
        self.current_index = (self.current_index + 1) % (len(self.model["SMILES"]) + 1)
        if self.current_index == 0:
            self.current_index = 1
        self.parent.update_ui()

    def previous_molecule(self):
        """Move to the previous molecule and update display."""
        if self.current_index == self.total_index == 1:
            return
        self.current_index = (self.current_index - 1) % (len(self.model["SMILES"])+1)
        if self.current_index == 0:
            self.current_index = len(self.model["SMILES"])
        self.parent.update_ui()

    def new_molecule(self):
        """Create a new empty molecule entry in the csv dictionary and update the display."""
        self.model.add_row({}) # add empty row
        self.current_index = self.model.get_total_index()
        self.total_index = self.model.get_total_index()
        self.parent.update_ui() 

    def delete_molecule(self) -> None:
        """Delete the current molecule and all associated data (including constraints) from the csv dictionary and update the display.
        
        Returns None"""
        if self.total_index == 1:
            self.model.delete_row(0)
            self.current_index = 1
            self.total_index = 1
            self.new_molecule()
            return

        self.model.delete_row(self.current_index - 1)

        for constraint in ['each_dist', 'each_angle', 'each_dihedral']:
            if hasattr(self, constraint):
                constraint_dict = getattr(self, constraint)
                if self.current_index in constraint_dict:
                    del constraint_dict[self.current_index]

        if hasattr(self, 'selected_atoms'):
            self.selected_atoms = []

        self.total_index = self.model.get_total_index()
        self.current_index = (self.current_index - 1) % (self.total_index + 1)
        if self.current_index == 0:
            self.current_index = 1
        self.model.signals.updated.emit()
        self.parent.update_ui()

    def open_smiles_tutorial(self):
        """Open the SMILES tutorial, uses QWebEngineView to open the URL."""
        self.smiles_tutorial_viewer = SmilesTutorialViewer()

    def generate_csearch_command(self):
        """Generate the CSEARCH command based on the current model settings and copy it to clipboard."""
        command_args = self.parent.parent.worker.collect_csearch_params()
        command_str = "python3 -m aqme --csearch " + " ".join(f"--{k} '{v}'" for k,v in command_args.items())
        return command_str

class StreamCapture:
    """Captures output and emits it via signal in real-time."""
    def __init__(self, signal):
        self.signal = signal
        self.buffer = StringIO()
    
    def write(self, text):
        if text and text.strip():  # Only emit non-empty lines
            self.signal.emit(text.rstrip())
        self.buffer.write(text)
    
    def flush(self):
        pass  # Required for file-like object
    
    def getvalue(self):
        return self.buffer.getvalue()
    

class QThreadLogHandler(logging.Handler):
    """Custom logging handler that emits logs via Qt signal."""
    def __init__(self, signal):
        super().__init__()
        self.signal = signal
    
    def emit(self, record):
        msg = self.format(record)
        self.signal.emit(msg)

class FDCapture(threading.Thread):
    """Capture output at file descriptor level."""
    def __init__(self, fd, signal, stop_event):
        super().__init__(daemon=True)
        self.fd = fd
        self.signal = signal
        self.stop_event = stop_event
        
    def run(self):
        """Read from file descriptor and emit via signal."""
        while not self.stop_event.is_set():
            try:
                # Use select with timeout to check if data is available
                if select.select([self.fd], [], [], 0.1)[0]:
                    data = os.read(self.fd, 4096)
                    if data:
                        text = data.decode('utf-8', errors='replace')
                        # Split by lines and emit each
                        for line in text.splitlines():
                            if line.strip():
                                self.signal.emit(line)
            except Exception as e:
                # File descriptor might be closed
                break

class CSEARCHThread(QThread):
    result = Signal(str)
    error = Signal(str)
    finished_signal = Signal(str)
    confirm = Signal(str, str)
    
    def __init__(self, parent, model):
        super().__init__()
        self.parent_widget = parent
        self.model = model
        self._stop_requested = False
        self._killer_active = False

    def request_stop(self):
        """Request the thread to stop and kill any running subprocesses."""
        self._stop_requested = True
        self._killer_active = True
        self.result.emit("Stopping CSEARCH and terminating subprocesses...")
        
        # Start aggressive subprocess killing in a loop
        QTimer.singleShot(0, self.aggressive_kill_loop)

    def aggressive_kill_loop(self):
        """Aggressively kill CSEARCH-related processes in the destination directory."""
        kill_duration = 5  # Kill for 5 seconds
        start_time = time.time()
        killed_count = 0
        
        # Get the destination directory to match against
        target_dir = os.path.abspath(self.model.__get__("destination"))
        
        while time.time() - start_time < kill_duration and self._killer_active:
            try:
                # Check all processes on the system
                for proc in psutil.process_iter(['pid', 'name', 'cmdline', 'cwd']):
                    try:
                        proc_name = proc.info['name'].lower()
                        
                        # Check if it's a process we care about
                        if any(name in proc_name for name in ['crest', 'xtb', 'obabel', 'rdkit']):
                            # Get the working directory of the process
                            try:
                                proc_cwd = proc.cwd()
                                
                                # Check if process is running in our target directory or subdirectory
                                if proc_cwd.startswith(target_dir):
                                    self.result.emit(f"Killing: {proc.info['name']} (PID: {proc.info['pid']}) in {proc_cwd}")
                                    psutil.Process(proc.info['pid']).kill()
                                    killed_count += 1
                            except (psutil.AccessDenied, psutil.NoSuchProcess):
                                # If we can't get cwd, try checking command line for directory hints
                                cmdline = proc.info.get('cmdline', [])
                                if cmdline and any(target_dir in str(arg) for arg in cmdline):
                                    self.result.emit(f"Killing (by cmdline): {proc.info['name']} (PID: {proc.info['pid']})")
                                    psutil.Process(proc.info['pid']).kill()
                                    killed_count += 1
                                    
                    except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                        pass
            
            except Exception as e:
                logging.error(f"Error in kill loop: {e}")
            
            time.sleep(0.2)  # Check every 200ms for new processes
        
        self.result.emit(f"Kill loop finished. Killed {killed_count} processes.")

    def run(self) -> None:
        """Runs AQME CSEARCH in the background."""
        
        if self.model["input"] is None or self.model["input"] == "":
            self.error.emit("Please save the CSV file before running AQME.")
            self.result.emit("CSEARCH aborted: No input CSV file specified.")
            return
        
        if self._stop_requested:
            self.result.emit("CSEARCH cancelled before starting.")
            return
        
        command_args = self.collect_csearch_params()
        self.result.emit(f"Starting CSEARCH with parameters: {command_args}")
        
        # Set up logging capture
        log_handler = QThreadLogHandler(self.result)
        log_handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(message)s')
        log_handler.setFormatter(formatter)
        
        root_logger = logging.getLogger()
        aqme_logger = logging.getLogger('aqme')
        
        root_logger.addHandler(log_handler)
        aqme_logger.addHandler(log_handler)
        
        # Redirect stdout/stderr at file descriptor level
        stdout_pipe_read, stdout_pipe_write = os.pipe()
        stderr_pipe_read, stderr_pipe_write = os.pipe()
        
        # Save original file descriptors
        original_stdout_fd = os.dup(sys.stdout.fileno())
        original_stderr_fd = os.dup(sys.stderr.fileno())
        
        # Create stop event for capture threads
        stop_event = threading.Event()
        
        # Start capture threads
        stdout_capture = FDCapture(stdout_pipe_read, self.result, stop_event)
        stderr_capture = FDCapture(stderr_pipe_read, self.result, stop_event)
        stdout_capture.start()
        stderr_capture.start()
        
        try:
            # Redirect file descriptors
            os.dup2(stdout_pipe_write, sys.stdout.fileno())
            os.dup2(stderr_pipe_write, sys.stderr.fileno())
            
            if self.model["program"] == "rdkit":
                self.result.emit("Running RDKit conformer search...")
                if not self._stop_requested:
                    csearch(**command_args)
                    
            elif self.model["program"] == "crest":
                self.result.emit("Running CREST conformer search (this might take a while)...")
                if not self._stop_requested:
                    csearch(**command_args, **crest_command)
            
            if self._stop_requested:
                self.result.emit("CSEARCH was stopped.")
                self._killer_active = False
                return
            
            self.result.emit("CSEARCH completed successfully!")
            self.finished_signal.emit("CSEARCH run finished.")
            self.confirm.emit(
                "Would you like to open QPREP with the generated SDF file(s)?",
                self.model["destination"]
            )
            
        except Exception as e:
            if self._stop_requested:
                self.result.emit("CSEARCH stopped by user.")
                self._killer_active = False
            else:
                logging.exception(f"Error encountered during CSEARCH: {e}")
                self.error.emit(f"Error during CSEARCH: {str(e)}")
        
        finally:
            # Restore original file descriptors
            os.dup2(original_stdout_fd, sys.stdout.fileno())
            os.dup2(original_stderr_fd, sys.stderr.fileno())
            
            # Close pipe file descriptors
            os.close(original_stdout_fd)
            os.close(original_stderr_fd)
            os.close(stdout_pipe_write)
            os.close(stderr_pipe_write)
            
            # Stop capture threads
            stop_event.set()
            stdout_capture.join(timeout=1)
            stderr_capture.join(timeout=1)
            
            # Close read pipes
            os.close(stdout_pipe_read)
            os.close(stderr_pipe_read)
            
            # Remove logging handlers
            root_logger.removeHandler(log_handler)
            aqme_logger.removeHandler(log_handler)

    def collect_csearch_params(self) -> dict:
        """Collects csearch parameters from the model."""
        csearch_params = {}

        if self.model["program"] == "":
            self.model["program"] = "rdkit"
        if self.model["destination"] is None or self.model["destination"] == "":
            if self.model["input"] != "" and self.model["input"] is not None:
                self.model["destination"] = self.model["input"].replace(".csv", "_aqme")

        try:
            params = self.model
            for key, value in params.items():
                if key in general_command_default:
                    if value != general_command_default[key]:
                        csearch_params[key] = value
        except Exception as e:
            logging.warning(f"Error at collect_csearch_params: {e}")

        return csearch_params