import csv 
import os
import subprocess
import tempfile

from PySide6.QtWidgets import QMainWindow, QLabel, QVBoxLayout, QWidget, QPushButton, QHBoxLayout, QTextBrowser,QLineEdit, QTextEdit, QCheckBox, QInputDialog, QMessageBox, QSizePolicy, QFileDialog, QDialog, QTableWidget, QTableWidgetItem, QHeaderView, QApplication
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QPixmap, QKeySequence, QShortcut, QMouseEvent

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D, rdDepictor

import PubChemPy as pcp
# from openbabel import pybel as pb

class smiles_to_csv(QWidget):
    def __init__(self):
        super().__init__()
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self.resize(900,800)
        self.setWindowTitle("smiles2csv")

        self.initialize_csv_dictionary()
        self.current_index = 1
        self.total_index = len(self.csv_dictionary["SMILES"]) 

        self.top_layout = QHBoxLayout()
        self.left_layout = QVBoxLayout()
        self.right_layout = QVBoxLayout()

        self.control1_layout = QHBoxLayout()
        self.control2_layout = QHBoxLayout()    
        self.control3_layout = QHBoxLayout()

        self.right_layout.addLayout(self.control1_layout)
        self.bottom_layout = QHBoxLayout()
        self.main_layout = QVBoxLayout()

        self.top_layout.addLayout(self.left_layout, 3)  
        self.top_layout.addLayout(self.right_layout, 1)  
        self.main_layout.addLayout(self.top_layout, 3)
        self.main_layout.addLayout(self.bottom_layout, 1)
        self.setLayout(self.main_layout)

        self.index_and_total_label = QLabel(f"{self.current_index}/{self.total_index}", self)
        self.control1_layout.addWidget(self.index_and_total_label)

        self.import_button = QPushButton("Import", self)
        self.import_button.clicked.connect(self.import_file)
        self.control1_layout.addWidget(self.import_button)

        self.new_molecule_button = QPushButton("New Molecule", self)
        self.new_molecule_button.clicked.connect(self.new_molecule)
        self.control1_layout.addWidget(self.new_molecule_button)

        self.show_all_button = QPushButton("Show All", self)
        self.show_all_button.clicked.connect(self.show_csv)
        self.control1_layout.addWidget(self.show_all_button)

        self.molecule_label = QLabel("Enter SMILES in the box on the right...", self)
        self.molecule_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.molecule_label.setStyleSheet("background-color: white; border: 1px solid black;")
        self.molecule_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)  
        self.molecule_label.setMinimumSize(250, 200)  
        self.left_layout.addWidget(self.molecule_label)

        self.atom_electron_label = QLabel(self.molecule_label)
        self.atom_electron_label.setFixedSize(150, 30)
        self.atom_electron_label.setStyleSheet("background-color: rgba(255, 255, 255, 0); border: none; font-size: 10px;")

        self.smiles_input = QTextEdit(self)
        self.smiles_input.setPlaceholderText("Enter SMILES here or search PubChem below...")
        self.smiles_input.setStyleSheet("padding: 5px;")
        self.smiles_input.setAutoFillBackground(True)
        self.smiles_input.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.right_layout.addWidget(self.smiles_input, 1)

        self.smiles_input.textChanged.connect(self.enumerate_smiles)
        self.smiles_input.textChanged.connect(self.display_molecule)
        self.smiles_input.textChanged.connect(lambda: self.update_smiles_csv())
        self.smiles_input.textChanged.connect(lambda: self.display_molecule(self.show_numbered_atoms_toggle.isChecked()))
        self.smiles_input.textChanged.connect(lambda: self.get_multiplicity())
        self.smiles_input.textChanged.connect(lambda: self.get_charge())
        self.smiles_input.textChanged.connect(lambda: self.atom_electron_label.setText(f" Electrons: {self.get_number_electrons()}\n Atoms: {self.get_number_atoms()}"))

        QShortcut(QKeySequence(Qt.Key_Escape), self.smiles_input, self.smiles_input.clearFocus)
        palette = self.smiles_input.palette()
        palette.setColor(self.smiles_input.backgroundRole(), self.palette().color(self.backgroundRole()))
        self.smiles_input.setPalette(palette)

        self.show_numbered_atoms_toggle = QCheckBox("Display atom labels", self)
        self.show_numbered_atoms_toggle.setChecked(False)
        self.show_numbered_atoms_toggle.stateChanged.connect(lambda: self.display_molecule(self.show_numbered_atoms_toggle.isChecked()))
        self.control2_layout.addWidget(self.show_numbered_atoms_toggle)

        self.search_pubchem_input = QLineEdit(self)
        self.search_pubchem_input.setPlaceholderText("Search PubChem...")
        self.search_pubchem_input.returnPressed.connect(self.find_smiles_from_PubChem)
        self.control2_layout.addWidget(self.search_pubchem_input)

        self.right_layout.addLayout(self.control2_layout)
        self.right_layout.addLayout(self.control3_layout)

        self.delete_button = QPushButton("Delete", self)
        self.delete_button.clicked.connect(self.delete_molecule)
        self.control3_layout.addWidget(self.delete_button)

        self.previous_button = QPushButton("Previous", self)
        self.previous_button.clicked.connect(self.previous_molecule)
        QShortcut(QKeySequence(Qt.Key_Left), self, self.previous_molecule)
        self.control3_layout.addWidget(self.previous_button)

        self.next_button = QPushButton("Next", self)
        self.next_button.clicked.connect(self.next_molecule)
        QShortcut(QKeySequence(Qt.Key_Right), self, self.next_molecule)
        self.control3_layout.addWidget(self.next_button)

        self.smiles_output = QTextBrowser(self)
        self.smiles_output.setStyleSheet("padding: 5px;")
        self.smiles_output.setAutoFillBackground(True)
        palette = self.smiles_output.palette()
        palette.setColor(self.smiles_output.backgroundRole(), self.palette().color(self.backgroundRole()))
        self.smiles_output.setPalette(palette)
        self.smiles_output.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.right_layout.addWidget(self.smiles_output, 1)

        self.properties_table = QTableWidget(self)
        self.properties_table.setRowCount(9)
        self.properties_table.setColumnCount(2)
        self.properties_table.horizontalHeader().setVisible(False)
        self.properties_table.verticalHeader().setVisible(False)
        self.properties_table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.properties_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.properties_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.properties_table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)

        # Apply styling to make the table visually appealing
        self.properties_table.setStyleSheet("""
            QTableWidget {
            background-color: #f9f9f9;
            border: 1px solid #dcdcdc;
            gridline-color: #dcdcdc;
            font-family: Arial, sans-serif;
            font-size: 12px;
            }
            QTableWidget::item {
            padding: 5px;
            }
            QTableWidget::item:selected {
            background-color: #cce7ff;
            color: #000;
            }
        """)
        self.right_layout.addWidget(self.properties_table, 2)
        self.update_properties()

        self.run_button = QPushButton("Run AQME", self)
        # self.run_button.clicked.connect(self.run_aqme)
        self.bottom_layout.addWidget(self.run_button)



# CSV RELATED FUNCTIONS

    def initialize_csv_dictionary(self):
        """Initialize the dictionary that will be used to create the CSV file, as per what aqme can read in. Key values are columns, values are lists where each index corresponds to a row (molecule)."""
        self.csv_dictionary = {
                "SMILES": [""],
                "code_name": [""],
                "charge": [""],
                "multiplicity": [""],
                "constraints_atoms": [""],
                "constraints_dist": [""],
                "constraints_angle": [""],
                "constraints_dihedral": [""],
                "complex_type": [""],
                "geom": [""]
        }

    def show_csv(self):
        """Display the csv_dictionary in table format in a new popup window. The user can edit the file, which automatically saves."""
        self.csv_table = QDialog(self)
        self.csv_table.setWindowTitle("CSV Data")
        self.csv_table.resize(800, 600)
        self.csv_table_layout = QVBoxLayout()
        self.csv_table.setLayout(self.csv_table_layout)

        self.table_widget = QTableWidget(self)
        self.table_widget.setRowCount(len(self.csv_dictionary["SMILES"]))
        self.table_widget.setColumnCount(len(self.csv_dictionary.keys()))
        self.table_widget.setHorizontalHeaderLabels(self.csv_dictionary.keys())

        for row in range(len(self.csv_dictionary["SMILES"])):
            for col, key in enumerate(self.csv_dictionary.keys()):
                item = QTableWidgetItem(str(self.csv_dictionary[key][row]))
                self.table_widget.setItem(row, col, item)

        self.table_widget.itemChanged.connect(self.update_csv_dictionary)

        self.csv_table_layout.addWidget(self.table_widget)
        self.csv_table.exec()

    def update_csv_dictionary(self, item):
        """Update the csv_dictionary when the user edits the table."""
        row = item.row()
        col = item.column()
        key = list(self.csv_dictionary.keys())[col]
        self.csv_dictionary[key][row] = item.text()

    def update_smiles_csv(self):
        """Update the csv_dictionary with the current enumerated SMILES string."""
        if self.smiles_input.toPlainText() != "":
            self.csv_dictionary["SMILES"][self.current_index - 1] = self.smiles_input.toPlainText() 
        else:
            return

    def save_csv_file(self):
        """Save the csv_dictionary to a file."""
        if not self.file_name:
            self.file_name, _ = QFileDialog.getSaveFileName(self, "Save CSV File", "", "CSV Files (*.csv)")
            if not self.file_name:
                return

        with open(self.file_name, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(self.csv_dictionary.keys())
            for i in range(len(self.csv_dictionary["SMILES"])):
                writer.writerow([self.csv_dictionary[key][i] for key in self.csv_dictionary.keys()])

# SMILES HANDILING FUNCTIONS
    def find_smiles_from_PubChem(self):
        """Opens a search bar and takes input to find canonical smiles
        from pubchem and pastes them into the input_smiles box."""
        if self.search_pubchem_input.text() != "":
            search_text = self.search_pubchem_input.text()
        else:
            search_text, ok = QInputDialog.getText(self, "PubChem Search", "Enter molecule name:")

        try:
            compound = pcp.get_compounds(search_text, 'name')[0]
            smiles = compound.isomeric_smiles
            current_text = self.smiles_input.toPlainText() + "." if self.smiles_input.toPlainText() != "" else ""
            new_text = current_text  + smiles
            self.smiles_input.setText(new_text)
            self.search_pubchem_input.clear()
        except Exception as e:
            print(f"Trouble finding SMILES from PubChem: {e}")
            QMessageBox.warning(self, "Error", "Could not find SMILES from PubChem.")

    def enumerate_smiles(self):
        """Enumerate all possible SMILES strings for a molecule."""
        if not self.smiles_input or self.smiles_input.toPlainText() == "":
            self.smiles_output.setText("Please enter SMILES in the box above.")
            return
            
        try:
            mol = Chem.MolFromSmiles(self.smiles_input.toPlainText())
            mol = Chem.AddHs(mol)
            for i, atom in enumerate(mol.GetAtoms()):
                atom.SetAtomMapNum(i + 1)
            self.smiles_output.setText(Chem.MolToSmiles(mol))
        except Exception as e:
            self.smiles_output.setText(f"Error numbering atoms: {str(e)}")

    def resizeEvent(self, event):
        """Handle window resize events to refresh the molecule display."""
        super().resizeEvent(event)
        self.display_molecule(self.show_numbered_atoms_toggle.isChecked())

    def display_molecule(self,checked = None):
        """Display the molecule in the molecule_label using rdkit.Chem.Draw module"""
        if not self.smiles_input or self.smiles_input.toPlainText() == "":
            self.molecule_label.setText("No molecule to display.")
            self.atom_coords = None
            return

        try:
            mol = Chem.MolFromSmiles(self.smiles_input.toPlainText())
            if mol is None:
                raise ValueError("Invalid SMILES string.")
            
            mol = Chem.AddHs(mol)
            for i, atom in enumerate(mol.GetAtoms()):
                atom.SetAtomMapNum(i + 1)
            
            width = self.molecule_label.width()
            height = self.molecule_label.height()
            drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
            
            if checked:
                drawer.drawOptions().noAtomLabels = False
            if not checked:
                drawer.drawOptions().noAtomLabels = True

            drawer.drawOptions().bondLineWidth = 1.5

            highlight_atoms = []
            highlight_colors = {}
            
            if hasattr(self, 'selected_atoms') and self.selected_atoms:
                highlight_atoms = [atom_idx - 1 for atom_idx in self.selected_atoms]
                highlight_colors = {idx: (1, 1, 0) for idx in highlight_atoms}

            drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms, highlightAtomColors=highlight_colors)
            drawer.FinishDrawing()
            drawer.WriteDrawingText("/tmp/molecule.svg")
            
            # Save atom coordinates
            self.atom_coords = [drawer.GetDrawCoords(i) for i in range(mol.GetNumAtoms())]

            pixmap = QPixmap("/tmp/molecule.svg")
            if self.molecule_label is not None:
                self.molecule_label.setPixmap(pixmap)
        except Exception as e:
            if self.molecule_label is not None:
                self.molecule_label.setText(f"Error displaying molecule: {str(e)}")

    def mousePressEvent(self, event: QMouseEvent):
        """Handle mouse press events to select atoms and add constraints.
        The logic is to check if the mouse press event is within the molecule_label"""
        if event.button() == Qt.MouseButton.LeftButton:
            pos = event.position()
            if self.molecule_label and self.molecule_label.geometry().contains(pos.toPoint()):
                x = pos.x() - self.molecule_label.x()
                y = pos.y() - self.molecule_label.y()
                selected_atom = self.get_atom_at_position(x, y)
                if selected_atom is not None:
                    self.handle_atom_selection(selected_atom)
                    self.display_molecule(self.show_numbered_atoms_toggle.isChecked())  

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
                if len(self.smiles_input.toPlainText()) < 30:
                    if (coord.x - x) ** 2 + (coord.y - y) ** 2 < 100: 
                        return idx + 1
                elif len(self.smiles_input.toPlainText()) > 30 :
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
            self.display_molecule(self.show_numbered_atoms_toggle.isChecked())  # update the highlighted atom
            return
        
        self.selected_atoms.append(atom_idx)
        self.display_molecule(self.show_numbered_atoms_toggle.isChecked())  # update the highlighted atom
        
        if len(self.selected_atoms) == 1:
            self.constraint_type, ok = QInputDialog.getItem(
                self, 
                "Select Constraint Type", 
                "Choose constraint type:", 
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
            mol = Chem.MolFromSmiles(self.smiles_input.toPlainText())
            mol = Chem.AddHs(mol) 
            atom1_type = mol.GetAtomWithIdx(atom1 - 1).GetSymbol()
            atom2_type = mol.GetAtomWithIdx(atom2 - 1).GetSymbol()
            distance, ok = QInputDialog.getDouble(self, "Distance Constraint", f"Enter distance between {atom1_type}:{atom1} and {atom2_type}:{atom2} (Å):")
            if ok:
                if distance <= 0:
                    QMessageBox.warning(self, "Invalid Distance", "Distance must be greater than 0.", QMessageBox.StandardButton.Ok, QMessageBox.StandardButton.Ok)
                    return self.add_distance_constraints()
                
            if not hasattr(self, 'each_dist'):
                self.each_dist = []
            self.each_dist.append([atom1, atom2, distance])
            while len(self.csv_dictionary["constraints_dist"]) <= self.current_index:
                self.csv_dictionary["constraints_dist"].append([])
            self.csv_dictionary["constraints_dist"][self.current_index - 1] = self.each_dist
            self.update_properties()
            self.selected_atoms = []

    def add_angle_constraints(self):
        """Add angle constraints between three selected atoms."""
        if len(self.selected_atoms) == 3:
            atom1, atom2, atom3 = self.selected_atoms
            mol = Chem.MolFromSmiles(self.smiles_input.toPlainText())
            mol = Chem.AddHs(mol) 
            atom1_type = mol.GetAtomWithIdx(atom1 - 1).GetSymbol()
            atom2_type = mol.GetAtomWithIdx(atom2 - 1).GetSymbol()
            atom3_type = mol.GetAtomWithIdx(atom3 - 1).GetSymbol()
            angle, ok = QInputDialog.getDouble(self, "Angle Constraint", f"Enter angle between {atom1_type}:{atom1}, {atom2_type}:{atom2}, and {atom3_type}:{atom3} (°):")
            if ok:
                if angle <= 0 or angle >= 360:
                    QMessageBox.warning(self, "Invalid Angle", "Angle must be between 0 and 360 degrees.", QMessageBox.StandardButton.Ok, QMessageBox.StandardButton.Ok)
                    return self.addAngleConstraints()

            if not hasattr(self, 'each_angle'):
                self.each_angle = []
            self.each_angle.append([atom1, atom2, atom3, angle])
            while len(self.csv_dictionary["constraints_angle"]) <= self.current_index:
                self.csv_dictionary["constraints_angle"].append([])
            self.csv_dictionary["constraints_angle"][self.current_index - 1] = self.each_angle
            self.update_properties()
            self.selected_atoms = []

    def add_dihedral_constraints(self):
        """Add dihedral constraints between four selected atoms."""
        if len(self.selected_atoms) == 4:
            atom1, atom2, atom3, atom4 = self.selected_atoms
            mol = Chem.MolFromSmiles(self.smiles_input.toPlainText())
            mol = Chem.AddHs(mol)
            atom1_type = mol.GetAtomWithIdx(atom1 - 1).GetSymbol()
            atom2_type = mol.GetAtomWithIdx(atom2 - 1).GetSymbol()
            atom3_type = mol.GetAtomWithIdx(atom3 - 1).GetSymbol()
            atom4_type = mol.GetAtomWithIdx(atom4 - 1).GetSymbol()
            dihedral, ok = QInputDialog.getDouble(self, "Dihedral Constraint", f"Enter dihedral angle between {atom1_type}:{atom1}, {atom2_type}:{atom2}, {atom3_type}:{atom3}, and {atom4_type}:{atom4} (°):")
            if ok:
                if dihedral < -180 or dihedral > 180:
                    QMessageBox.warning(self, "Invalid Dihedral", "Dihedral angle must be between -180 and 180 degrees.", QMessageBox.StandardButton.Ok, QMessageBox.StandardButton.Ok)
                    return self.addDihedralConstraints()

            if not hasattr(self, 'each_dihedral'):
                self.each_dihedral = []
            self.each_dihedral.append([atom1, atom2, atom3, atom4, dihedral])
            self.csv_dictionary["constraints_dihedral"][self.current_index - 1] = self.each_dihedral
            self.update_properties()
            self.selected_atoms = []

# PROPERTIES FUNCTIONS
    def update_properties(self):
        """Update the properties box with:
            "code_name",
            "charge",
            "multiplicity",
            "constraints_atoms",
            "constraints_dist",
            "constraints_angle",
            "constraints_dihedral",
            "complex_type",
            "geom" assigned to the current smiles index displayed. These values are stored in the csv_dictionary."""
        self.code_name = self.csv_dictionary["code_name"][self.current_index - 1]
        self.charge = self.csv_dictionary["charge"][self.current_index - 1]
        self.multiplicity = self.csv_dictionary["multiplicity"][self.current_index - 1]
        self.constraints_atoms = self.csv_dictionary["constraints_atoms"][self.current_index - 1]
        self.constraints_dist = self.csv_dictionary["constraints_dist"][self.current_index - 1]
        self.constraints_angle = self.csv_dictionary["constraints_angle"][self.current_index - 1]
        self.constraints_dihedral = self.csv_dictionary["constraints_dihedral"][self.current_index - 1]
        self.complex_type = self.csv_dictionary["complex_type"][self.current_index - 1]
        self.geom = self.csv_dictionary["geom"][self.current_index - 1]

        self.properties_table.setRowCount(9)  
        properties = [
            ("code_name", self.code_name),
            ("charge", self.charge),
            ("multiplicity", self.multiplicity),
            ("constraints_atoms", self.constraints_atoms),
            ("constraints_dist", self.constraints_dist),
            ("constraints_angle", self.constraints_angle),
            ("constraints_dihedral", self.constraints_dihedral),
            ("complex_type", self.complex_type),
            ("geom", self.geom),
        ]

        for row, (property_name, value) in enumerate(properties):
            self.properties_table.setItem(row, 0, QTableWidgetItem(property_name))
            self.properties_table.setItem(row, 1, QTableWidgetItem(str(value)))
        
    def get_multiplicity(self):
        """"""
        if self.csv_dictionary["multiplicity"][self.current_index - 1] != "":
            return
        elif self.csv_dictionary["SMILES"][self.current_index - 1] == "":
            return
        else:
            try:
                smiles = self.csv_dictionary["SMILES"][self.current_index - 1]
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    mult = Descriptors.NumRadicalElectrons(mol) + 1 
                    self.csv_dictionary["multiplicity"][self.current_index - 1] = mult
            except Exception:
                return 

    def get_charge(self):
        """"""
        if self.csv_dictionary["charge"][self.current_index - 1] != "":
            return
        elif self.csv_dictionary["SMILES"][self.current_index - 1] == "":
            return
        else:
            try:
                smiles = self.csv_dictionary["SMILES"][self.current_index - 1]
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    mol = Chem.AddHs(mol)
                    charge = Chem.GetFormalCharge(mol)
                    self.csv_dictionary["charge"][self.current_index - 1] = charge
            except Exception:
                return 

    def get_number_atoms(self):
        """Get number of atoms in the molecule (highest enumerated number with hydrogens)"""
        if self.smiles_input.toPlainText() == "":
                return 0
        try:
            mol = Chem.MolFromSmiles(self.smiles_input.toPlainText())
            mol = Chem.AddHs(mol)
            return mol.GetNumAtoms()
        except Exception:
            return 
        
    def get_number_electrons(self):
        """Electrons counted using RDKit Chem module:
        electrons = sum of atomic numbers - formal charge"""
        if self.smiles_input.toPlainText() == "":
                return 0
        try:
            mol = Chem.MolFromSmiles(self.smiles_input.toPlainText())
            mol = Chem.AddHs(mol)
            charge = Chem.GetFormalCharge(mol)
            electrons = sum([atom.GetAtomicNum() for atom in mol.GetAtoms()]) - charge
            return electrons
        except Exception:
            return 

# NAVIGATION FUNCTIONS
    def next_molecule(self):
        """Move to the next index in csv dictionary and update display."""
        if self.current_index == self.total_index == 1:
            return
        self.current_index = (self.current_index + 1) % (len(self.csv_dictionary["SMILES"]) +1)
        if self.current_index == 0:
            self.current_index = 1
        self.smiles_input.clear()
        self.update_display()

    def previous_molecule(self):
        """Move to the previous molecule and update display."""
        if self.current_index == self.total_index == 1:
            return
        self.current_index = (self.current_index - 1) % (len(self.csv_dictionary["SMILES"])+1)
        if self.current_index == 0:
            self.current_index = len(self.csv_dictionary["SMILES"])
        self.smiles_input.clear()
        self.update_display()

    def new_molecule(self):
        """Create a new empty molecule entry in the csv dictionary and update the display."""
        for key in self.csv_dictionary.keys():
            self.csv_dictionary[key].append("")
        self.current_index = len(self.csv_dictionary["SMILES"])  
        self.total_index = len(self.csv_dictionary["SMILES"])
        self.update_display()

    def delete_molecule(self):
        """Delete the current molecule from the csv dictionary and update the display."""
        if self.total_index == 1:
            return
        for key in self.csv_dictionary.keys():
            self.csv_dictionary[key].pop(self.current_index - 1)
        self.total_index = len(self.csv_dictionary["SMILES"])
        self.current_index = (self.current_index - 1) % (len(self.csv_dictionary["SMILES"]) + 1)
        if self.current_index == 0:
            self.current_index = 1
        self.update_display()

# UI ELEMENTS UPDATE FUNCTIONS
    def index_and_total_label_update(self):
        self.index_and_total_label.setText(f"{self.current_index}/{self.total_index}")

    def update_display(self):
        """Update all UI elements to reflect the current molecule's data."""
        self.smiles_input.clear()
        self.smiles_input.setText(self.csv_dictionary["SMILES"][self.current_index - 1])
        self.index_and_total_label_update()
        self.get_multiplicity()
        self.get_charge()
        self.update_properties()


    # def closeEvent(self, event):
    #     """Handle the close event to save the csv_dictionary to a file."""
    #     self.save_csv_file()
    #     event.accept()
    

# CHEMDRAW FUNCTIONS
    def import_file(self):
        """Import an SDF or ChemDraw file, extract SMILES, and display them."""
        self.initialize_csv_dictionary()
        file_name, _ = QFileDialog.getOpenFileName(self, "Import File", "", "ChemDraw Files (*.cdx);;SDF Files (*.sdf)")
        if not file_name:
            return

        try:
            smiles_list = []

            if file_name.endswith(".sdf"):
                mol_supplier = Chem.SDMolSupplier(file_name)
                for mol in mol_supplier:
                    if mol is not None:
                        smiles = Chem.MolToSmiles(mol)
                        smiles_list.append(smiles)

            elif file_name.endswith(".cdx"): # there is little to no support for cdxml files in obabel
                
                # This does not seem to work and is a known (?) issue with pybel as far as my research goes.
                # mol = next(pb.readfile("cdx", file_name))
                # with tempfile.NamedTemporaryFile(suffix='.sdf', delete=False) as temp_sdf:
                #     sdf_path = temp_sdf.name
                # mol.write("sdf", sdf_path, overwrite=True)

                # if not os.path.exists(sdf_path):
                #     print(f"Failed to write SDF file to {sdf_path}")
                #     return


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

            for index, smiles in enumerate(smiles_list):
                if index < len(self.csv_dictionary["SMILES"]):
                    self.csv_dictionary["SMILES"][index] = smiles
                else:
                    for key in self.csv_dictionary.keys():
                        self.csv_dictionary[key].append("")
                    self.csv_dictionary["SMILES"][-1] = smiles
            self.total_index = len(self.csv_dictionary["SMILES"])
            self.update_display()

        except ImportError as e:
            print(f"Error importing required module: {e}")
        except IOError as e:
            print(f"Error reading or writing file: {e}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

# AQME RUN FUNCTIONS


# MISC FOR LATER
    def find_metal_atom(self):
            """Find metal atoms in the SMILES string, display a statement in the log_console_box and 
            enable complex_type option if any are found."""
            smiles = self.smiles_input.toPlainText()
            mol = Chem.MolFromSmiles(smiles)
            metal_atoms = [] # for batch jobs such as CSV inputs with many SMILES
            transition_metals = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo',
                                'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
                                'Hg', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
            for atom in mol.GetAtoms():
                if atom.GetSymbol() in transition_metals:
                    metal_atoms.append(atom.GetSymbol())
            if len(metal_atoms) > 0:
                self.log_console_box.setText(f"Transition metal atoms detected in SMILES: {metal_atoms}. Select complex_type.")
                self.complex_type_combobox.setEnabled(True)
            return 

    def update_command(self):
        """Update the command based on smiles, constraints and chosen options.
        
        Requires expanding to include all options."""
        smiles = self.smiles_output.toPlainText()

        name = self.code_name_label.text().split(":")[-1].strip()
        nprocs = 12
        program = "crest" if self.program_crest_cb.isChecked() else "rdkit"
        constraints = []
        if hasattr(self, 'constraints_distance'):
            constraints_distance = self.constraints_distance.get(self.current_smiles_index, [])
            constraints.append(f'--constraints_dist "{constraints_distance}"')
        else:
            pass
        if hasattr(self, 'constraints_angle'):
            constraints_angle = self.constraints_angle.get(self.current_smiles_index, [])
            constraints.append(f'--constraints_angle "{constraints_angle}"')
        else:
            pass
        if hasattr(self, 'constraints_dihedral'):
            constraints_dihedral = self.constraints_dihedral.get(self.current_smiles_index, [])
            constraints.append(f'--constraints_dihedral "{constraints_dihedral}"')
        else:
            pass
        if self.complex_type_combobox != "":
            complex_type = self.complex_type_combobox.currentText()

        
        constraints_str = " ".join(constraints)
        if "[" in smiles:
            command = f'python3 -m aqme --csearch --name "{name}" --nprocs {nprocs} --program {program} {constraints_str} --smi "{smiles}"'
            if complex_type:
                command = f'python3 -m aqme --csearch --name "{name}" --nprocs {nprocs} --program {program} {constraints_str} --smi "{smiles}" --complex_type {complex_type}'
        else:
            command = "No valid SMILES to generate command."

        if self.csv_toggle.isChecked():
            file_name = self.file_name if self.file_name else "No csv detected"
            command_csv = f'python3 -m aqme --csearch --program {program} --nprocs {nprocs} --input {file_name} '
            self.command_line.setText(command_csv)
        else:
            self.command_line.setText(command)

    def copy_command_to_clipboard(self):
        """Copy the generated AQME command to the clipboard."""
        clipboard = QApplication.clipboard()
        clipboard.setText(self.command_line.text())
        QMessageBox.information(self, "Command Copied", "Command copied to clipboard.")

    def redo_coordinates(self):
        """unsure fully how this works in rdEditor but i liked the 
        function so i am trying to replicate it """
        rdDepictor.SetPreferCoordGen(True)
        self.display_molecule(self.show_numbered_atoms_toggle.isChecked())
