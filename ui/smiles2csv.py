import logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
import os
import csv
import subprocess
import tempfile

from PySide6.QtWidgets import  QLabel, QVBoxLayout, QWidget, QPushButton, QHBoxLayout, QTextBrowser,QLineEdit, QTextEdit, QCheckBox, QInputDialog, QMessageBox, QSizePolicy, QFileDialog, QTableWidget, QTableWidgetItem, QHeaderView, QApplication, QComboBox, QSpinBox, QStyle, QTableWidgetItem, QFrame, QGridLayout, QDoubleSpinBox
from PySide6.QtCore import Qt, QProcess
from PySide6.QtGui import QPixmap, QKeySequence, QShortcut, QMouseEvent, QIcon, QDoubleValidator, QTextCursor

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D, rdDepictor
import rdkit.rdBase

from utils import  pubchem2smiles, smiles2enumerate, smiles2numatoms, smiles2numelectrons, smiles2findmetal, smiles2charge, smiles2multiplicity
import ui.resources.icons as icons

from models.smiles2csv_model import csv_dictionary as csv_model
from models.smiles2csv_command import (
    general_command_dictionary as gen_command,
    summ_command_dictionary as summ_command,
    fullmonte_command_dictionary as fullmonte_command,
    crest_command_dictionary as crest_command
    )
from controllers.smiles2csv_controller import csv_controller as control

class smiles_to_csv(QWidget):
    def __init__(self):
        super().__init__()
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self.resize(900,800)
        self.setWindowTitle("smiles2csv")
        QShortcut(QKeySequence(Qt.Key_Escape), self, lambda: self.clear_focus_on_inputs())
        self.current_index = 1
        self.file_name = None
        self.total_index = len(csv_model["SMILES"]) 
        self.smiles_w_metal = [] 

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
        self.import_button.setFixedWidth(65)
        self.import_button.clicked.connect(lambda: (logging.debug("at import_button >>> self.import_file()"), self.import_file()))
        self.control1_layout.addWidget(self.import_button)

        self.new_molecule_button = QPushButton("New Molecule", self)
        self.new_molecule_button.clicked.connect(lambda: (logging.debug("at new_molecule_button >>> self.new_molecule()"), self.new_molecule()))
        self.control1_layout.addWidget(self.new_molecule_button)

        self.show_all_button = QPushButton("Show All", self)
        self.show_all_button.setFixedWidth(65)
        self.show_all_button.clicked.connect(control.show_csv)
        self.control1_layout.addWidget(self.show_all_button)

        self.save_csv_button = QPushButton(self)
        self.save_csv_button.setToolTip("Save CSV file")
        self.save_csv_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_DialogSaveButton))
        self.save_csv_button.setStyleSheet("font-size: 12px; color: black;")
        self.save_csv_button.clicked.connect(lambda: (logging.debug("at save_csv_button >>> save_csv"), self.save_csv_file()))
        self.control1_layout.addWidget(self.save_csv_button)

        self.how_to_label = QLabel("Click to select atoms and add constraints, click again to deselect.", self)
        self.how_to_label.setStyleSheet("background-color: rgba(255, 255, 255, 0); border: none; font-size: 10px; color: black;")
        self.how_to_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.left_layout.addWidget(self.how_to_label)

    # THIS IS WHERE THE MOLECULE IS DISPLAYED
        self.molecule_label = QLabel("Enter SMILES in the box on the right...", self)
        self.molecule_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.molecule_label.setStyleSheet("background-color: #f8f8f8; border: 1px solid black; color: black; font-size: 12px;")
        self.molecule_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)  
        self.molecule_label.setMinimumSize(250, 200)  
        self.left_layout.addWidget(self.molecule_label)

        self.atom_electron_label = QLabel(self.molecule_label)
        self.atom_electron_label.setFixedSize(150, 30)
        self.atom_electron_label.setStyleSheet("background-color: rgba(255, 255, 255, 0); border: none; font-size: 10px; color: black;")

        self.log_box_label = QLabel(self.molecule_label)
        self.log_box_label.setFixedSize(550, 20)
        self.log_box_label.setStyleSheet("background-color: rgba(255, 255, 255, 0); border: none; font-size: 10px; color: black;")
        self.molecule_label.resizeEvent = lambda event: self.log_box_label.move(5, self.molecule_label.height() - self.log_box_label.height())

    # VARIOUS INPUTS
        self.smiles_input = QTextEdit(self)
        self.smiles_input.setPlaceholderText("Enter SMILES here or search PubChem below...")
        self.smiles_input.setStyleSheet("border: 1px solid #dcdcdc;")
        self.smiles_input.setAutoFillBackground(True)
        self.smiles_input.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.right_layout.addWidget(self.smiles_input, 1)

        self.smiles_input.textChanged.connect(lambda: control.update_smiles_model(self.smiles_input.toPlainText(), self.current_index))

        self.show_numbered_atoms_toggle = QCheckBox("Show atom labels", self)
        self.show_numbered_atoms_toggle.setChecked(False)
        self.show_numbered_atoms_toggle.stateChanged.connect(lambda: self.display_molecule(self.show_numbered_atoms_toggle.isChecked()))
        self.control2_layout.addWidget(self.show_numbered_atoms_toggle)

        self.search_pubchem_input = QLineEdit(self)
        self.search_pubchem_input.setPlaceholderText("Search PubChem...")
        self.search_pubchem_input.returnPressed.connect(self.smiles_from_pubchem)
        self.control2_layout.addWidget(self.search_pubchem_input)

        self.search_pubchem_advanced_button = QPushButton(self)
        self.search_pubchem_advanced_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_FileDialogDetailedView))
        # this is for Getting a full results list for common compound names (https://pubchempy.readthedocs.io/en/latest/guide/searching.html)
        # probs wont incorporate advanced search types but idk yet
        # might also add searching for SIDs along with name and CID IDK yet need to ask Juanvi
        self.search_pubchem_advanced_button.setToolTip("Advanced PubChem Search")
        self.control2_layout.addWidget(self.search_pubchem_advanced_button)

        self.right_layout.addLayout(self.control2_layout)
        self.right_layout.addLayout(self.control3_layout)

    # CONTROL BUTTONS
        self.delete_button = QPushButton("Delete", self)
        self.delete_button.clicked.connect(self.delete_molecule)
        self.control3_layout.addWidget(self.delete_button)

        self.previous_button = QPushButton("Previous", self)
        self.previous_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_ArrowLeft))
        self.previous_button.clicked.connect(lambda: (logging.debug("at previous_button >>> self.previous_molecule"), self.previous_molecule()))
        QShortcut(QKeySequence(Qt.Key_Left), self, self.previous_molecule)
        self.control3_layout.addWidget(self.previous_button)

        self.next_button = QPushButton("Next ", self)
        self.next_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_ArrowRight))
        self.next_button.setLayoutDirection(Qt.LayoutDirection.RightToLeft)
        self.next_button.clicked.connect(lambda: (logging.debug("at next_button >>> self.next_molecule"), self.next_molecule()))
        QShortcut(QKeySequence(Qt.Key_Right), self, self.next_molecule)
        self.control3_layout.addWidget(self.next_button)

    # OUTPUTS
        self.smiles_output = QTextBrowser(self)
        self.smiles_output.setStyleSheet("border: 1px solid #dcdcdc;")
        self.smiles_output.setAutoFillBackground(True)
        self.smiles_output.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.right_layout.addWidget(self.smiles_output, 1)

        self.properties_table = QTableWidget(self)
        self.properties_table.setRowCount(9)
        self.properties_table.setColumnCount(2)
        self.properties_table.horizontalHeader().setVisible(False)
        self.properties_table.verticalHeader().setVisible(False)
        self.properties_table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        self.is_programmatic_update = False
        self.update_properties()
        logging.debug("at __init__ >>> self.update_properties()")
        
        self.properties_table.itemChanged.connect(lambda item: (
            logging.debug("at properties_table >>> handle_property_change"),
            self.handle_property_change(item) if not self.is_programmatic_update 
            and item.flags() & Qt.ItemFlag.ItemIsEditable  
            and item.column() == 1 else None
        ))

        self.properties_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.properties_table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.properties_table.setStyleSheet("""
            QTableWidget {
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

        self.shell_output = QTextBrowser(self)
        self.shell_output.setStyleSheet("background-color: black; color: white; padding: 10px; border: 1px solid #444;")
        self.shell_output.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.shell_output.setReadOnly(True)
        self.bottom_layout.addWidget(self.shell_output, 2)
    
        # AQME SETUP
        self.aqme_setup_grid = QGridLayout()
        self.bottom_layout.addLayout(self.aqme_setup_grid, 1)

        # Row 0 - Program selection
        self.program_label = QLabel("Select CSEARCH program:", self)
        self.program_label.setStyleSheet("font-size: 12px; color: black;")
        self.aqme_setup_grid.addWidget(self.program_label, 0, 0)
        self.program_combo = QComboBox(self)
        self.program_combo.addItems(["RDKit", "CREST", "GOAT*"])
        self.program_combo.setStyleSheet("font-size: 12px; color: black;")
        self.aqme_setup_grid.addWidget(self.program_combo, 0, 1)

        # Row 1 - Number of processors
        self.nprocs_label = QLabel("Number of processors:", self)
        self.nprocs_label.setStyleSheet("font-size: 12px; color: black;")
        self.aqme_setup_grid.addWidget(self.nprocs_label, 1, 0)
        self.nprocs_input = QSpinBox(self)
        self.nprocs_input.setRange(1, 40)
        self.nprocs_input.setValue(8)
        self.nprocs_input.setStyleSheet("font-size: 12px; color: black;")
        self.nprocs_input.setToolTip("Number of processors to use for the calculation.\nOnly relevant for CREST and GOAT*.")
        self.nprocs_input.setEnabled(self.program_combo.currentText() == "CREST")
        self.aqme_setup_grid.addWidget(self.nprocs_input, 1, 1)
        self.program_combo.currentTextChanged.connect(lambda text: self.nprocs_input.setEnabled(text == "CREST"))

        # Row 2 - Stack size
        self.stacksize_label = QLabel("Stack size (GB):", self)
        self.stacksize_label.setStyleSheet("font-size: 12px; color: black;")
        self.aqme_setup_grid.addWidget(self.stacksize_label, 2, 0)
        self.stacksize_input = QSpinBox(self)
        self.stacksize_input.setRange(1, 8)
        self.stacksize_input.setValue(1)
        self.stacksize_input.setStyleSheet("font-size: 12px; color: black;")
        self.aqme_setup_grid.addWidget(self.stacksize_input, 2, 1)

        # Row 3 - Output directory
        self.output_dir_label = QLabel("Output directory:", self)
        self.output_dir_label.setStyleSheet("font-size: 12px; color: black;")
        self.aqme_setup_grid.addWidget(self.output_dir_label, 3, 0)
        self.output_dir_input = QLineEdit(self)
        self.output_dir_input.setPlaceholderText("Select output directory...")
        self.output_dir_button = QPushButton(self)
        dir_icon = self.style().standardIcon(QStyle.StandardPixmap.SP_DirIcon)
        self.output_dir_button.setIcon(dir_icon)
        self.output_dir_button.clicked.connect(self.select_output_directory)
        output_dir_layout = QHBoxLayout()
        output_dir_layout.addWidget(self.output_dir_input)
        output_dir_layout.addWidget(self.output_dir_button)
        self.aqme_setup_grid.addLayout(output_dir_layout, 3, 1)

        # Row 4 - copy command/save csv
        self.copy_command_button = QPushButton("Copy Command", self)
        self.copy_command_button.setStyleSheet("font-size: 12px; color: black;")
        self.copy_command_button.clicked.connect(lambda: (logging.debug("at copy_command_button >>> copy_command"), self.copy_command_to_clipboard()))
        self.aqme_setup_grid.addWidget(self.copy_command_button, 4,0)

        # Row 5 - Run button
        self.run_button = QPushButton("Run AQME", self)
        self.run_button.setStyleSheet("font-size: 12px; color: black; background-color: lightblue;")
        self.run_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay))
        self.run_button.setFixedHeight(55)
        self.run_button.clicked.connect(lambda: (logging.debug("at run_button >>> self.run_aqme()"), self.run_aqme()))
        self.aqme_setup_grid.addWidget(self.run_button, 4, 1, 2, 1)

        for widget in [self.smiles_input, self.smiles_output, self.properties_table,  self.shell_output, self.molecule_label, self.atom_electron_label]:
            palette = widget.palette()
            palette.setColor(widget.backgroundRole(), self.palette().color(self.backgroundRole()))
            widget.setPalette(palette)

        # ADVANCED SETTINGS 

        self.advanced_settings_button = QPushButton("Show Advanced Settings", self)
        self.advanced_settings_button.setCheckable(True)
        self.advanced_settings_button.clicked.connect(lambda: (logging.debug("at advanced_settings_button >>> toggle_panel"), self.toggle_panel()))
        self.advanced_settings_button.setStyleSheet("font-size: 12px; color: black;")
        self.aqme_setup_grid.addWidget(self.advanced_settings_button, 5, 0)

        self.advanced_panel = QFrame()
        self.advanced_panel.setFixedHeight(0)
        self.main_layout.addWidget(self.advanced_panel)

        advanced_layout = QGridLayout()
        self.advanced_panel.setLayout(advanced_layout)

        self.sample_size_label = QLabel("Sample size:", self)
        self.sample_size_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.sample_size_label, 0, 0)
        self.sample_size_input = QSpinBox(self)
        self.sample_size_input.setRange(1, 500)
        self.sample_size_input.setValue(25)
        self.sample_size_input.setStyleSheet("font-size: 12px; color: black;")
        self.sample_size_input.setToolTip("Number of conformers to keep after the initial RDKit sampling.\nThey are selected using a combination of RDKit energies and Butina clustering.")
        advanced_layout.addWidget(self.sample_size_input, 0, 1)

        self.auto_sample_label = QLabel("Auto sample level:", self)
        self.auto_sample_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.auto_sample_label, 1, 0)
        self.auto_sample_combo = QComboBox(self)
        self.auto_sample_combo.addItems(["low", "mid", "high", "false"])
        self.auto_sample_combo.setCurrentText("mid")
        self.auto_sample_combo.setStyleSheet("font-size: 12px; color: black;")
        self.auto_sample_combo.setToolTip("Apply automatic calculation of the number of conformers generated initially with RDKit. \nThis number of conformers is initially generated and then reduced to the number specified in --sample with different filters. \nOptions:\n• Low: Base multiplier = 5, max confs = 100\n• Mid: Base multiplier = 10, max confs = 250\n• High: Base multiplier = 20, max confs = 500\n• False: Use the number specified in --sample")
        advanced_layout.addWidget(self.auto_sample_combo, 1, 1)

        self.energy_window_label = QLabel("Energy window (kcal/mol):", self)
        self.energy_window_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.energy_window_label, 2, 0)
        self.energy_window_input = QDoubleSpinBox(self)
        self.energy_window_input.setRange(1.0, 50.0)
        self.energy_window_input.setValue(5.0)
        self.energy_window_input.setStyleSheet("font-size: 12px; color: black;")
        self.energy_window_input.setToolTip("Energy window in kcal/mol to discard conformers\n(i.e. if a conformer is more than the E window compared to the most stable conformer).")
        advanced_layout.addWidget(self.energy_window_input, 2, 1)

        self.initial_energy_threshold_label = QLabel("Initial E threshold (kcal/mol):", self)
        self.initial_energy_threshold_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.initial_energy_threshold_label, 3, 0)
        self.initial_energy_threshold_input = QLineEdit(self)
        self.initial_energy_threshold_input.setText("0.0001")
        self.initial_energy_threshold_input.setValidator(QDoubleValidator())
        self.initial_energy_threshold_input.setStyleSheet("font-size: 12px; color: black;")
        self.initial_energy_threshold_input.setToolTip("Energy difference in kcal/mol between unique conformers for the first filter of only E.")
        advanced_layout.addWidget(self.initial_energy_threshold_input, 3, 1)

        self.energy_threshold_label = QLabel("Energy threshold (kcal/mol):", self)
        self.energy_threshold_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.energy_threshold_label, 4, 0)

        self.energy_threshold_input = QLineEdit(self)
        self.energy_threshold_input.setText("0.25")
        self.energy_threshold_input.setValidator(QDoubleValidator())
        self.energy_threshold_input.setStyleSheet("font-size: 12px; color: black;")
        self.energy_threshold_input.setToolTip("Energy difference in kcal/mol between unique conformers for the second filter of E + RMS")
        advanced_layout.addWidget(self.energy_threshold_input, 4, 1)

        self.rms_threshold_label = QLabel("RMS threshold:", self) #(kcal/mol)?
        self.rms_threshold_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.rms_threshold_label, 0, 3)

        self.rms_threshold_input = QLineEdit(self)
        self.rms_threshold_input.setText("0.25")
        self.rms_threshold_input.setValidator(QDoubleValidator())
        self.rms_threshold_input.setStyleSheet("font-size: 12px; color: black;")
        self.rms_threshold_input.setToolTip("RMS difference between unique conformers for the second filter of E + RMS")
        advanced_layout.addWidget(self.rms_threshold_input, 0, 4)

        self.opt_steps_rdkit_label = QLabel("RDKit opt steps:", self)
        self.opt_steps_rdkit_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.opt_steps_rdkit_label, 1, 3)

        self.opt_steps_rdkit_input = QSpinBox(self)
        self.opt_steps_rdkit_input.setRange(1, 10000)
        self.opt_steps_rdkit_input.setValue(1000)
        self.opt_steps_rdkit_input.setStyleSheet("font-size: 12px; color: black;")
        self.opt_steps_rdkit_input.setToolTip("Max cycles used in RDKit optimizations. ")
        advanced_layout.addWidget(self.opt_steps_rdkit_input, 1, 4)


        csv_model.signals.updated.connect(self.update_ui)
        self.update_ui()

    def toggle_panel(self):
        expanded_height = 130
        if self.advanced_settings_button.isChecked():
            self.advanced_panel.setFixedHeight(expanded_height)
            self.resize(self.width(), self.height() + expanded_height)
            self.advanced_settings_button.setText("Hide Advanced Settings")
        else:
            self.advanced_panel.setFixedHeight(0)
            self.resize(self.width(), self.height() - expanded_height)
            self.advanced_settings_button.setText("Show Advanced Settings")


    # def handle_smiles_change(self):
    #     # logging.debug("at smiles_input >>> handling SMILES text change")
    #     smiles = self.smiles_input.toPlainText()
    #     enumerated_smiles = smiles2enumerate(smiles)
    #     if enumerated_smiles is not None:
    #         self.smiles_output.setText(enumerated_smiles)
    #     self.display_molecule(self.show_numbered_atoms_toggle.isChecked())
    #     self.find_metal_atom(smiles)


    def handle_property_change(self, item):
        """Handle changes to editable properties and update the csv_dictionary."""
        logging.debug("at handle_property_change >>> handling property change")
        if item.row() == 0:  
            csv_model["code_name"][self.current_index - 1] = item.text()
        elif item.row() == 1:  
            csv_model["charge"][self.current_index - 1] = item.text()
            # Track user-defined charge
            if not hasattr(self, 'user_defined_charge'):
                self.user_defined_charge = {}
            self.user_defined_charge[self.current_index] = item.text()
        elif item.row() == 2:  
            csv_model["multiplicity"][self.current_index - 1] = item.text()
            # Track user-defined multiplicity
            if not hasattr(self, 'user_defined_multiplicity'):
                self.user_defined_multiplicity = {}
            self.user_defined_multiplicity[self.current_index] = item.text()









# SMILES HANDILING FUNCTIONS (this to certain extent is also UI handling)
    def smiles_from_pubchem(self):
        """Searches PubChem for a compound using its CID or name and inserts the canonical SMILES into the input box."""
        code_name = self.search_pubchem_input.text()
        if not code_name:
            QMessageBox.warning(self, "Error", "Please enter a valid CID, CAS or compound name.")
        try:
            smiles = pubchem2smiles(code_name)
        except IndexError:
            pixmap = QPixmap(icons.red_path)
            icon = QIcon(pixmap)
            msgBox = QMessageBox(self)
            msgBox.setWindowTitle("Error")
            msgBox.setText("No compound found for the given input. Please check the CID or name.")
            msgBox.setWindowIcon(icon)
            msgBox.setIconPixmap(pixmap)
            msgBox.setStandardButtons(QMessageBox.Ok)
            msgBox.exec()
        except Exception as e:
            QMessageBox.warning(self, "Error", f"An error occurred: {str(e)}")

        current_text = csv_model["SMILES"][self.current_index - 1]
        new_text = f"{current_text}.{smiles}" if current_text else smiles
        self.smiles_input.setText(new_text)

        if not csv_model["code_name"][self.current_index - 1]:
            csv_model["code_name"][self.current_index - 1] = code_name
        self.update_properties()
        self.update_ui()
        self.search_pubchem_input.clear()

    def resizeEvent(self, event):
        """Handle window resize events to refresh the molecule display."""
        super().resizeEvent(event)
        self.display_molecule(self.show_numbered_atoms_toggle.isChecked())

    def display_molecule(self,checked = None):
        """Display the molecule in the molecule_label using rdkit.Chem.Draw module"""
        logging.debug("at display_molecule >>> displaying molecule")
        rdkit.rdBase.DisableLog('rdApp.*')
        rdDepictor.SetPreferCoordGen(True)

        smiles = csv_model["SMILES"][self.current_index - 1]
        if not smiles:
            self.molecule_label.setText("No molecule to display.")
            self.atom_coords = None
            return

        try:
            mol = Chem.MolFromSmiles(smiles)
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
                highlight_colors = {idx: (0.5176, 0.6314, 0.9922) for idx in highlight_atoms}

            drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms, highlightAtomColors=highlight_colors)
            drawer.FinishDrawing()
            drawer.WriteDrawingText("/tmp/molecule.svg")
            
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
        logging.debug("at mousePressEvent >>> handling mouse press event")
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
        logging.debug("at get_atom_at_position >>> getting atom at position")
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
        logging.debug("at handle_atom_selection >>> handling atom selection")
        if not hasattr(self, 'selected_atoms'):
            self.selected_atoms = []

        if atom_idx in self.selected_atoms:  # deselect if already selected
            self.selected_atoms.remove(atom_idx)
            self.display_molecule(self.show_numbered_atoms_toggle.isChecked())  # update the highlighted atom
            return
        
        self.selected_atoms.append(atom_idx)
        self.display_molecule(self.show_numbered_atoms_toggle.isChecked())  # update the highlighted atom
        atom1 = self.selected_atoms[0]
        mol = Chem.MolFromSmiles(csv_model["SMILES"][self.current_index - 1])
        mol = Chem.AddHs(mol) 
        atom1_type = mol.GetAtomWithIdx(atom1 - 1).GetSymbol()

        if len(self.selected_atoms) == 1:
            self.constraint_type, ok = QInputDialog.getItem(
                self, 
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
            distance, ok = QInputDialog.getDouble(self, "Distance Constraint", f"Enter distance between {atom1_type}:{atom1} and {atom2_type}:{atom2} (Å):")
            if ok:
                if distance <= 0:
                    QMessageBox.warning(self, "Invalid Distance", "Distance must be greater than 0.", QMessageBox.StandardButton.Ok, QMessageBox.StandardButton.Ok)
                    return  self.add_distance_constraints()
                
                if not hasattr(self, 'each_dist'):
                    self.each_dist = {self.current_index: []}
                if self.current_index not in self.each_dist:
                    self.each_dist[self.current_index] = []
                self.each_dist[self.current_index].append([atom1, atom2, distance])
                csv_model["constraints_dist"][self.current_index - 1] = str(self.each_dist[self.current_index])
                self.update_properties()
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
            angle, ok = QInputDialog.getDouble(self, "Angle Constraint", f"Enter angle between {atom1_type}:{atom1}, {atom2_type}:{atom2}, and {atom3_type}:{atom3} (°):")
            if ok:
                if angle <= 0 or angle >= 360:
                    QMessageBox.warning(self, "Invalid Angle", "Angle must be between 0 and 360 degrees.", QMessageBox.StandardButton.Ok, QMessageBox.StandardButton.Ok)
                    return self.addAngleConstraints()
                if not hasattr(self, 'each_angle'):
                    self.each_angle = {self.current_index: []}
                if self.current_index not in self.each_angle:
                    self.each_angle[self.current_index] = []
                self.each_angle[self.current_index].append([atom1, atom2, atom3, angle])
                csv_model["constraints_angle"][self.current_index - 1] = str(self.each_angle[self.current_index])
                self.update_properties()
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
            dihedral, ok = QInputDialog.getDouble(self, "Dihedral Constraint", f"Enter dihedral angle between {atom1_type}:{atom1}, {atom2_type}:{atom2}, {atom3_type}:{atom3}, and {atom4_type}:{atom4} (°):")
            if ok:
                if dihedral < -180 or dihedral > 180:
                    QMessageBox.warning(self, "Invalid Dihedral", "Dihedral angle must be between -180 and 180 degrees.", QMessageBox.StandardButton.Ok, QMessageBox.StandardButton.Ok)
                    return self.addDihedralConstraints()
                if not hasattr(self, 'each_dihedral'):
                    self.each_dihedral = {self.current_index: []}
                if self.current_index not in self.each_dihedral:
                    self.each_dihedral[self.current_index] = []
                self.each_dihedral[self.current_index].append([atom1, atom2, atom3, atom4, dihedral])
                csv_model["constraints_dihedral"][self.current_index - 1] = str(self.each_dihedral[self.current_index])
                self.update_properties()
                self.selected_atoms = []
            else:
                self.selected_atoms = []
                return




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
        old_value = self.is_programmatic_update
        self.is_programmatic_update = True

        try:
            self.properties_table.setRowCount(9)  
            properties = [
            ("code_name", csv_model["code_name"][self.current_index - 1]),
            ("charge", csv_model["charge"][self.current_index - 1]),
            ("multiplicity", csv_model["multiplicity"][self.current_index - 1]),
            ("constraints_atoms", csv_model["constraints_atoms"][self.current_index - 1]),
            ("constraints_dist", csv_model["constraints_dist"][self.current_index - 1]),
            ("constraints_angle", csv_model["constraints_angle"][self.current_index - 1]),
            ("constraints_dihedral", csv_model["constraints_dihedral"][self.current_index - 1]),
            ("complex_type", csv_model["complex_type"][self.current_index - 1]),
            ("geom", csv_model["geom"][self.current_index - 1]),
            ]
            # Temporarily disconnect the signal to avoid any callbacks
            self.properties_table.blockSignals(True)
            
            for row, (property_name, value) in enumerate(properties):
                self.properties_table.setItem(row, 0, QTableWidgetItem(property_name))
                self.properties_table.setItem(row, 1, QTableWidgetItem(str(value)))
            
            for row in range(self.properties_table.rowCount()):
                item = self.properties_table.item(row, 0)
                if item:
                    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)

            for row in range(self.properties_table.rowCount()):
                item = self.properties_table.item(row, 1)
                if item:
                    if row in [0, 1, 2]:
                        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsEditable)
                    else:  
                        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                        
            self.properties_table.blockSignals(False)
            self.file_name = None  # Reset file name to indicate unsaved changes idk if this works as intended tho 
        finally:
            self.is_programmatic_update = old_value



# NAVIGATION FUNCTIONS 
    def next_molecule(self):
        """Move to the next index in csv dictionary and update display."""
        if self.current_index == self.total_index == 1:
            return
        self.current_index = (self.current_index + 1) % (len(csv_model["SMILES"]) +1)
        if self.current_index == 0:
            self.current_index = 1
        self.update_ui()

    def previous_molecule(self):
        """Move to the previous molecule and update display."""
        if self.current_index == self.total_index == 1:
            return
        self.current_index = (self.current_index - 1) % (len(csv_model["SMILES"])+1)
        if self.current_index == 0:
            self.current_index = len(csv_model["SMILES"])
        self.update_ui()

    def new_molecule(self):
        """Create a new empty molecule entry in the csv dictionary and update the display."""
        for key in csv_model.keys():
            csv_model[key].append("")
        self.current_index = len(csv_model["SMILES"]) 
        self.total_index = len(csv_model["SMILES"])
        self.update_ui()

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

        self.total_index = len(csv_model["SMILES"])
        self.current_index = (self.current_index - 1) % (self.total_index + 1)
        if self.current_index == 0:
            self.current_index = 1
        self.update_ui()


# UI ELEMENTS UPDATE FUNCTIONS (This will definitely stay but might need to rewrite quite heavily)

    def update_ui(self):
        """Update all UI elements to reflect the current molecule's data."""
        smiles = csv_model["SMILES"][self.current_index - 1]
        self.smiles_input.blockSignals(True)
        self.smiles_input.setText(smiles)
        self.smiles_output.setText(smiles2enumerate(smiles))

        cursor = self.smiles_input.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        self.smiles_input.setTextCursor(cursor)

        self.smiles_input.blockSignals(False)

        self.display_molecule(self.show_numbered_atoms_toggle.isChecked())

        self.index_and_total_label_update()

        try:
            num_atoms = smiles2numatoms(smiles)
            num_electrons = smiles2numelectrons(smiles)
            self.atom_electron_label.setText(f" Electrons: {num_electrons}\n Atoms: {num_atoms}")
        except Exception as e:
            logging.error(f"Error calculating atom/electron counts: {e}")
            self.atom_electron_label.setText(" Electrons: 0\n Atoms: 0")

        self.update_properties()
        




    def index_and_total_label_update(self):
        self.index_and_total_label.setText(f"{self.current_index}/{self.total_index}")

    def clear_focus_on_inputs(self):
        self.smiles_input.clearFocus()
        self.search_pubchem_input.clearFocus()

    def closeEvent(self, event):
        """Handle the close event to prompt saving the CSV file, with the added icon."""
        pixmap = QPixmap(icons.red_path)
        icon = QIcon(pixmap)
        
        if self.file_name is None: 
            msgBox = QMessageBox(self)
            msgBox.setWindowTitle("Save CSV")
            msgBox.setText("Would you like to save the CSV file before exiting?")
            msgBox.setWindowIcon(icon)
            msgBox.setIconPixmap(pixmap)
            msgBox.setStandardButtons(QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel)
            reply = msgBox.exec()
            
            if reply == QMessageBox.Save:
                if not self.save_csv_file():
                    event.ignore()  # Prevent closing if saving fails
                else:
                    event.accept()
            elif reply == QMessageBox.Discard:
                event.accept()
            else:
                event.ignore()
        else:
            event.accept()







# need to move this one i fear
    def save_csv_file(self):
        """Save the csv_dictionary to a file."""
        error_pixmap = QPixmap(icons.red_path)
        error_icon = QIcon(error_pixmap)
        success_pixmap = QPixmap(icons.green_path)
        success_icon = QIcon(success_pixmap)
        for index in range(len(csv_model["code_name"])):
            if csv_model["code_name"][index] == "":
                msgBox = QMessageBox(self)
                msgBox.setWindowTitle("Error")
                msgBox.setText("Please enter a code name for all molecules before saving.")
                msgBox.setWindowIcon(error_icon)
                msgBox.setIconPixmap(error_pixmap)
                msgBox.setStandardButtons(QMessageBox.Ok)
                msgBox.exec()
                return False  # for the closing event
        file_name, _ = QFileDialog.getSaveFileName(self, "Save CSV File", "", "CSV Files (*.csv)")
        self.file_name = file_name 
        if not self.file_name:
            self.file_name = None
            return False  # for the closing event
        with open(self.file_name, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(csv_model.keys())
            for i in range(len(csv_model["SMILES"])):
                writer.writerow([csv_model[key][i] for key in csv_model.keys()])   
        msgBox = QMessageBox(self)
        msgBox.setWindowTitle("File Saved")
        msgBox.setText("CSV file saved successfully.")
        msgBox.setWindowIcon(success_icon)
        msgBox.setIconPixmap(success_pixmap)
        msgBox.setStandardButtons(QMessageBox.Ok)
        msgBox.exec()
        return True  # again for the closing event








# CHEMDRAW FUNCTIONS
    def import_file(self):
        """Import an SDF or ChemDraw file, extract SMILES, and display them."""
        file_name, _ = QFileDialog.getOpenFileName(self, "Import File", "", "ChemDraw Files (*.cdx *.cdxml);;SDF Files (*.sdf);;CSV files (*.csv)")
        if not file_name:
            return
        try:
            for key in csv_model.keys():
                csv_model[key].clear()
            smiles_list = []

            if file_name.endswith(".sdf"):
                mol_supplier = Chem.SDMolSupplier(file_name)
                for mol in mol_supplier:
                    if mol is not None:
                        smiles = Chem.MolToSmiles(mol)
                        smiles_list.append(smiles)

            elif file_name.endswith(".cdx"): 
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
                # fuckass solution but it works for now i guess (TODO: make it better) also idk how to handle rows i did not anticipate as of rn
                for key in csv_model.keys():
                    csv_model[key].pop(0) 
                with open(file_name, 'r') as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader: 
                        csv_model["SMILES"].append(row["SMILES"])
                        csv_model["code_name"].append(row["code_name"])
                        csv_model["charge"].append(row["charge"])
                        csv_model["multiplicity"].append(row["multiplicity"])
                        csv_model["constraints_atoms"].append(row["constraints_atoms"])
                        csv_model["constraints_dist"].append(row["constraints_dist"])
                        csv_model["constraints_angle"].append(row["constraints_angle"])
                        csv_model["constraints_dihedral"].append(row["constraints_dihedral"])
                        csv_model["complex_type"].append(row["complex_type"])
                        csv_model["geom"].append(row["geom"])
                        csv_model.signals.updated.emit()
                    
                self.total_index = len(csv_model["SMILES"])
                if self.total_index > 0:
                    self.current_index = 1  # assuming indices start from 1
                    self.update_properties()
                    self.update_ui()
                return

            else:
                QMessageBox.warning(self, "Error", "Unsupported file format. Please select a ChemDraw, SDF, or CSV file.")
                return
            if not smiles_list:
                QMessageBox.warning(self, "Error", "No valid SMILES found in the file.")
                return

            for index, smiles in enumerate(smiles_list):
                if index < len(csv_model["SMILES"]):
                    csv_model["SMILES"][index] = smiles
                    csv_model["charge"][index] = smiles2charge(smiles)
                    csv_model["multiplicity"][index] = smiles2multiplicity(smiles)
                else:
                    for key in csv_model.keys():
                        csv_model[key].append("")
                    csv_model["SMILES"][-1] = smiles
                    csv_model["charge"][-1] = smiles2charge(smiles)
                    csv_model["multiplicity"][-1] = smiles2multiplicity(smiles)
            self.total_index = len(csv_model["SMILES"])

            for index in range(self.total_index):
                csv_model["code_name"][index] = f"mol_{index + 1}"
                self.current_index = 1
                
            csv_model.signals.updated.emit()
            self.update_properties()

        except ImportError as e:
            print(f"Error importing required module: {e}")
        except IOError as e:
            print(f"Error reading or writing file: {e}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

# AQME RUN SETUP FUNCTIONS
    def copy_command_to_clipboard(self):
        """Copy the generated AQME command to the clipboard."""
        clipboard = QApplication.clipboard()
        command = self.aqme_rungen()
        clipboard.setText(command)
        pixmap = QPixmap(icons.green_path)
        icon = QIcon(pixmap)
        msg = QMessageBox(self)
        msg.setWindowTitle("Command Copied")
        msg.setText("Command copied to clipboard.")
        msg.setWindowIcon(icon)
        msg.setIconPixmap(pixmap)
        msg.exec()

    def select_output_directory(self):
        """Select the output directory for AQME results."""
        self.output_directory = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if self.output_directory:
            self.output_dir_input.setText(f"{self.output_directory}")
        else:
            return

# AQME RUN FUNCTIONS (still under construction, based almnost solely on the pyside6 book example)
    def run_aqme(self):
        """Run AQME with the generated command using QProcess in the file_name directory."""
        for value in csv_model["SMILES"]:
            if value == "":
                msgBox = QMessageBox(self)
                msgBox.setIconPixmap(QPixmap(icons.red_path))
                msgBox.setWindowTitle("Warning")
                msgBox.setText("No SMILES found. Please enter a SMILES string.")
                msgBox.setStandardButtons(QMessageBox.StandardButton.Ok)
                msgBox.exec()
                return
            else:
                pass

        changed_charge_multiplicity = []
        if hasattr(self, 'user_defined_multiplicity'):
            for key in self.user_defined_charge.keys():
                changed_charge_multiplicity.append(key)
        if hasattr(self, 'user_defined_charge'):
            for key in self.user_defined_multiplicity.keys():
                changed_charge_multiplicity.append(key)
        if not self.smiles_w_metal:
            pass
        else:
            for index in self.smiles_w_metal:
                if index not in changed_charge_multiplicity:
                    continue
            msgBox = QMessageBox(self)
            msgBox.setIconPixmap(QPixmap(icons.red_path))
            msgBox.setWindowTitle("Warning")
            msgBox.setText(f"Please check the charge and multiplicity for the transition metal complex(es).")
            msgBox.setStandardButtons(QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            msgBox.setButtonText(QMessageBox.StandardButton.Yes, "Proceed Anyway")
            msgBox.setButtonText(QMessageBox.StandardButton.No, "Go Back")
            result = msgBox.exec()
            if result == QMessageBox.StandardButton.No:
                return
            else:
                pass

        if self.file_name is None:
            self.save_csv_file() 
            if not self.file_name:
                return
        command = self.aqme_rungen()
        file_directory = os.path.dirname(self.file_name)

        # self.shell_output.clear()
        self.process = QProcess(self)
        self.process.setWorkingDirectory(file_directory)

        self.process.readyReadStandardOutput.connect(self.handle_stdout)
        self.process.readyReadStandardError.connect(self.handle_stderr)
        self.process.stateChanged.connect(self.handle_state)
        self.process.finished.connect(self.process_finished)

        if os.name == 'posix':
            self.process.start("bash", ["-c", f"{command}"])
        else:
            self.process.start(command)

    def handle_stdout(self):
        """Append standard output from the process to the shell_output."""
        output = bytes(self.process.readAllStandardOutput()).decode("utf8")
        if output:
            self.shell_output.append(output.strip())

    def handle_stderr(self):
        """Append standard error from the process to the shell_output."""
        error = bytes(self.process.readAllStandardError()).decode("utf8")
        if error:
            self.shell_output.append(f"Error: {error.strip()}")

    def handle_state(self, state):
        """Handle the state change of the process."""
        if state == QProcess.ProcessState.NotRunning:
            self.shell_output.append("AQME process finished.")
            self.run_button.setText("Run AQME")
            self.run_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay))
            self.run_button.setStyleSheet("font-size: 12px; color: black; background-color: lightblue;")
            self.run_button.clicked.disconnect()
            self.run_button.clicked.connect(self.run_aqme)
        elif state == QProcess.ProcessState.Running:
            self.shell_output.append("AQME process running...")
            self.run_button.setText("Stop AQME")
            self.run_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_MediaStop))
            self.run_button.setStyleSheet("font-size: 12px; color: black; background-color: #E06666;")
            self.run_button.clicked.disconnect()
            self.run_button.clicked.connect(self.stop_aqme)
        elif state == QProcess.ProcessState.Starting:
            self.shell_output.append("AQME process starting...")

    def stop_aqme(self):
        """Stop the AQME process if it is running."""
        if self.process and self.process.state() == QProcess.ProcessState.Running:
            self.process.kill()
            self.shell_output.append("AQME process stopped.")
            self.run_button.setText("Run AQME")
            self.run_button.setStyleSheet("font-size: 12px; color: black; background-color: lightblue;")
            self.run_button.clicked.disconnect()
            self.run_button.clicked.connect(self.run_aqme)
        else:
            self.shell_output.append("AQME process is not running.")

    def process_finished(self, exitCode):
        """Handle the process finish event and display a message box."""
        if exitCode == 0:
            pixmap = QPixmap(icons.green_path)
            icon = QIcon(pixmap)
            msgBox = QMessageBox(self)
            msgBox.setWindowTitle("AQME Run Completed")
            msgBox.setText("AQME run completed successfully.")
            msgBox.setWindowIcon(icon)
            msgBox.setIconPixmap(pixmap)
            msgBox.setStandardButtons(QMessageBox.Ok)
            msgBox.exec()
        else:
            pixmap = QPixmap(icons.red_path)
            icon = QIcon(pixmap)
            msgBox = QMessageBox(self)
            msgBox.setWindowTitle("AQME Run Failed")
            msgBox.setText("AQME run failed.")
            msgBox.setWindowIcon(icon)
            msgBox.setIconPixmap(pixmap)
            msgBox.setStandardButtons(QMessageBox.Ok)
            msgBox.exec()

    def aqme_rungen(self):
        """Update command w/ csv file."""
        program = self.program_combo.currentText()
        aqme_rungen = f'python -u -m aqme --csearch --program {program} --input {self.file_name} '
        if self.output_dir_input == "":
            aqme_rungen += f'--destination {self.output_dir_input.text()} '
        if self.nprocs_input:
            aqme_rungen += f'--nprocs {int(self.nprocs_input.text())} '
        if self.stacksize_input:
            aqme_rungen += f'--stacksize {int(self.stacksize_input.text())}GB '
        if self.sample_size_input:
            aqme_rungen += f'--sample {int(self.sample_size_input.text())} '
        if self.auto_sample_combo:
            aqme_rungen += f'--auto_sample {self.auto_sample_combo.currentText()} '
        if self.energy_window_input:
            aqme_rungen += f'--ewin_csearch {float(self.energy_window_input.text())} '
        if self.initial_energy_threshold_input:
            aqme_rungen += f'--initial_energy_threshold {float(self.initial_energy_threshold_input.text())} '
        if self.energy_threshold_input:
            aqme_rungen += f'--energy_threshold {float(self.energy_threshold_input.text())} '
        if self.rms_threshold_input:
            aqme_rungen += f'--rms_threshold {float(self.rms_threshold_input.text())} '
        if self.opt_steps_rdkit_input:
            aqme_rungen += f'--opt_steps_rdkit {int(self.opt_steps_rdkit_input.text())} '
        print(aqme_rungen)
        return aqme_rungen

# FOR LATER

# To do:
# - Add all the parameters to the actual command ... !!!!!
# - Add the complex_type option when TM is found !
# - Fix run aqme command  !!!!
# - expand the pubchem search  ... !

# add intermedaite plus transition state option in the show all 

# fucking dark mode man !!!!
