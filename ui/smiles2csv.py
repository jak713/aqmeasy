import logging
import csv
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
import os
import subprocess
import tempfile

from PySide6.QtWidgets import  QLabel, QVBoxLayout, QWidget, QPushButton, QHBoxLayout, QTextBrowser,QLineEdit, QTextEdit, QCheckBox, QMessageBox, QSizePolicy, QFileDialog, QTableWidget, QTableWidgetItem, QHeaderView, QApplication, QComboBox, QSpinBox, QStyle, QTableWidgetItem, QFrame, QGridLayout, QDoubleSpinBox
from PySide6.QtCore import Qt, QProcess, QEvent
from PySide6.QtGui import QPixmap, QKeySequence, QShortcut, QMouseEvent, QIcon, QDoubleValidator, QTextCursor, QIntValidator

from rdkit import Chem

from utils import  pubchem2smiles, smiles2enumerate, smiles2numatoms, smiles2numelectrons, smiles2charge, smiles2multiplicity, command2clipboard
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
        control.set_parent(self)

        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self.resize(900,800)
        self.setWindowTitle("smiles2csv")
        QShortcut(QKeySequence(Qt.Key_Escape), self, lambda: self.clear_focus_on_inputs())

        self.smiles_w_metal = [] # will take this out later

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

        self.index_and_total_label = QLabel(f"{control.current_index}/{control.total_index}", self)
        self.control1_layout.addWidget(self.index_and_total_label)

        self.import_button = QPushButton("Import", self)
        self.import_button.setFixedWidth(65)
        self.import_button.clicked.connect(lambda: (logging.debug("at import_button >>> self.import_file()"), self.import_file()))
        self.control1_layout.addWidget(self.import_button)

        self.new_molecule_button = QPushButton("New Molecule", self)
        self.new_molecule_button.clicked.connect(lambda: (logging.debug("at new_molecule_button >>> self.new_molecule()"), control.new_molecule()))
        self.control1_layout.addWidget(self.new_molecule_button)

        self.show_all_button = QPushButton("Show All", self)
        self.show_all_button.setFixedWidth(65)
        self.show_all_button.clicked.connect(control.show_csv)
        self.control1_layout.addWidget(self.show_all_button)

        self.save_csv_button = QPushButton(self)
        self.save_csv_button.setToolTip("Save CSV file")
        self.save_csv_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_DialogSaveButton))
        self.save_csv_button.setStyleSheet("font-size: 12px; color: black;")
        self.save_csv_button.clicked.connect(lambda: (logging.debug("at save_csv_button >>> save_csv"), control.save_csv_file()))
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

        # Pass the molecule_label to the controller
        self.left_layout.addWidget(self.molecule_label)

        self.atom_electron_label = QLabel(self.molecule_label)
        self.atom_electron_label.setFixedSize(150, 30)
        self.atom_electron_label.setStyleSheet("background-color: rgba(255, 255, 255, 0); border: none; font-size: 10px; color: black;")

        self.log_box_label = QLabel(self.molecule_label)
        self.log_box_label.setFixedSize(550, 20)
        self.log_box_label.setStyleSheet("background-color: rgba(255, 255, 255, 0); border: none; font-size: 10px; color: black;")
        self.molecule_label.resizeEvent = lambda event: self.log_box_label.move(5, self.molecule_label.height() - self.log_box_label.height() - 2)

    # VARIOUS INPUTS
        self.smiles_input = QTextEdit(self)
        self.smiles_input.setPlaceholderText("Enter SMILES here or search PubChem below...")
        self.smiles_input.setStyleSheet("border: 1px solid #dcdcdc;")
        self.smiles_input.setAutoFillBackground(True)
        self.smiles_input.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.smiles_input.installEventFilter(self)
        self.right_layout.addWidget(self.smiles_input, 1)

        self.smiles_input.textChanged.connect(lambda: control.update_smiles_model(self.smiles_input.toPlainText()))

        self.show_numbered_atoms_toggle = QCheckBox("Show atom labels", self)
        self.show_numbered_atoms_toggle.setChecked(False)
        self.show_numbered_atoms_toggle.stateChanged.connect(lambda: control.display_molecule(self.show_numbered_atoms_toggle.isChecked()))
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
        self.delete_button.clicked.connect(control.delete_molecule)
        self.control3_layout.addWidget(self.delete_button)

        self.previous_button = QPushButton("Previous", self)
        self.previous_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_ArrowLeft))
        self.previous_button.clicked.connect(lambda: (logging.debug("at previous_button >>> self.previous_molecule"), control.previous_molecule()))
        QShortcut(QKeySequence(Qt.Key_Left), self, control.previous_molecule)
        self.control3_layout.addWidget(self.previous_button)

        self.next_button = QPushButton("Next ", self)
        self.next_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_ArrowRight))
        self.next_button.setLayoutDirection(Qt.LayoutDirection.RightToLeft)
        self.next_button.clicked.connect(lambda: (logging.debug("at next_button >>> self.next_molecule"), control.next_molecule()))
        QShortcut(QKeySequence(Qt.Key_Right), self, control.next_molecule)
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
        self.program_combo.model().item(2).setEnabled(False)
        self.program_combo.currentTextChanged.connect(lambda text: gen_command.update({"program": text}))
        self.program_combo.setStyleSheet("font-size: 12px; color: black;")
        self.aqme_setup_grid.addWidget(self.program_combo, 0, 1)

        # Row 1 - Number of processors
        self.nprocs_label = QLabel("Number of processors:", self)
        self.nprocs_label.setStyleSheet("font-size: 12px; color: black;")
        self.aqme_setup_grid.addWidget(self.nprocs_label, 1, 0)
        
        self.nprocs_input = QSpinBox(self)
        self.nprocs_input.setRange(1, 40)
        self.nprocs_input.setValue(8)
        self.nprocs_input.valueChanged.connect(lambda: gen_command.update({"nprocs": self.nprocs_input.value()}))
        self.nprocs_input.setStyleSheet("font-size: 12px; color: black;")
        self.nprocs_input.setToolTip("Number of processors to use for the calculation.\nOnly relevant for CREST and GOAT*.")
        self.nprocs_input.setEnabled(self.program_combo.currentText() == "CREST")
        self.aqme_setup_grid.addWidget(self.nprocs_input, 1, 1)
        self.program_combo.currentTextChanged.connect(lambda text: self.nprocs_input.setEnabled(text == "CREST"))

        # Row 2 - Stack size
        self.stacksize_label = QLabel("Stack size:", self)
        self.stacksize_label.setStyleSheet("font-size: 12px; color: black;")
        self.aqme_setup_grid.addWidget(self.stacksize_label, 2, 0)
        self.stacksize_input = QSpinBox(self)
        self.stacksize_input.setRange(1, 8)
        self.stacksize_input.setValue(1)
        self.stacksize_input.valueChanged.connect(lambda: gen_command.update({"stacksize": f"{self.stacksize_input.value()}"}))
        self.stacksize_input.setSuffix(" GB")
        self.stacksize_input.setStyleSheet("font-size: 12px; color: black;")
        self.aqme_setup_grid.addWidget(self.stacksize_input, 2, 1)

        # Row 3 - Output directory
        self.output_dir_label = QLabel("Output directory:", self)
        self.output_dir_label.setStyleSheet("font-size: 12px; color: black;")
        self.aqme_setup_grid.addWidget(self.output_dir_label, 3, 0)
        self.output_dir_input = QLineEdit(self)
        self.output_dir_input.setPlaceholderText("Select output directory...")
        self.output_dir_input.setStyleSheet("font-size: 12px; color: black;")
        self.output_dir_input.textChanged.connect(lambda: gen_command.update({"destination": self.output_dir_input.text()}))

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
        self.copy_command_button.clicked.connect(lambda: self.success("Command copied to clipboard.") if command2clipboard(self.aqme_rungen()) == True else self.failure("Failed to copy command to clipboard."))
        self.aqme_setup_grid.addWidget(self.copy_command_button, 4,0)

        # Row 5 - Run button
        self.run_button = QPushButton("Run AQME", self)
        self.run_button.setStyleSheet("font-size: 12px; color: black; background-color: lightblue;")
        self.run_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay))
        self.run_button.setFixedHeight(55)
        self.run_button.clicked.connect(control.run_aqme)
        self.aqme_setup_grid.addWidget(self.run_button, 4, 1, 2, 1)

        for widget in [self.smiles_input, self.smiles_output, self.properties_table,  self.shell_output, self.molecule_label, self.atom_electron_label]:
            palette = widget.palette()
            palette.setColor(widget.backgroundRole(), self.palette().color(self.backgroundRole()))
            widget.setPalette(palette)

        # ADVANCED SETTINGS 

        # Column 1 & 2
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
        self.sample_size_input.valueChanged.connect(lambda: gen_command.update({"sample": self.sample_size_input.value()}))
        self.sample_size_input.setStyleSheet("font-size: 12px; color: black;")
        self.sample_size_input.setToolTip("Number of conformers to keep after the initial RDKit sampling.\nThey are selected using a combination of RDKit energies and Butina clustering.")
        advanced_layout.addWidget(self.sample_size_input, 0, 1)

        self.auto_sample_label = QLabel("Auto sample level:", self)
        self.auto_sample_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.auto_sample_label, 1, 0)
        self.auto_sample_combo = QComboBox(self)
        self.auto_sample_combo.addItems(["low", "mid", "high", "false"])
        self.auto_sample_combo.setCurrentText("mid")
        self.auto_sample_combo.currentTextChanged.connect(lambda text: gen_command.update({"auto_sample": text}))
        self.auto_sample_combo.setStyleSheet("font-size: 12px; color: black;")
        self.auto_sample_combo.setToolTip("Apply automatic calculation of the number of conformers generated initially with RDKit. \nThis number of conformers is initially generated and then reduced to the number specified in --sample with different filters. \nOptions:\n• Low: Base multiplier = 5, max confs = 100\n• Mid: Base multiplier = 10, max confs = 250\n• High: Base multiplier = 20, max confs = 500\n• False: Use the number specified in --sample")
        advanced_layout.addWidget(self.auto_sample_combo, 1, 1)

        self.energy_window_label = QLabel("Energy window:", self)
        self.energy_window_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.energy_window_label, 2, 0)
        self.energy_window_input = QDoubleSpinBox(self)
        self.energy_window_input.setRange(1.0, 50.0)
        self.energy_window_input.setValue(5.0)
        self.energy_window_input.valueChanged.connect(lambda: gen_command.update({"ewin_csearch": self.energy_window_input.value()}))
        self.energy_window_input.setStyleSheet("font-size: 12px; color: black;")
        self.energy_window_input.setToolTip("Energy window in kcal/mol to discard conformers\n(i.e. if a conformer is more than the E window compared to the most stable conformer).")
        advanced_layout.addWidget(self.energy_window_input, 2, 1)

        self.initial_energy_threshold_label = QLabel("Initial E threshold:", self)
        self.initial_energy_threshold_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.initial_energy_threshold_label, 3, 0)
        self.initial_energy_threshold_input = QLineEdit(self)
        self.initial_energy_threshold_input.setText("0.0001")
        self.initial_energy_threshold_input.setValidator(QDoubleValidator())
        self.initial_energy_threshold_input.textChanged.connect(lambda: gen_command.update({"initial_energy_threshold": self.initial_energy_threshold_input.text()}))
        self.initial_energy_threshold_input.setStyleSheet("font-size: 12px; color: black;")
        self.initial_energy_threshold_input.setToolTip("Energy difference in kcal/mol between unique conformers for the first filter of only E.")
        advanced_layout.addWidget(self.initial_energy_threshold_input, 3, 1)

        self.energy_threshold_label = QLabel("Energy threshold:", self)
        self.energy_threshold_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.energy_threshold_label, 4, 0)

        self.energy_threshold_input = QLineEdit(self)
        self.energy_threshold_input.setText("0.25")
        self.energy_threshold_input.setValidator(QDoubleValidator())
        self.energy_threshold_input.textChanged.connect(lambda: gen_command.update({"energy_threshold": self.energy_threshold_input.text()}))
        self.energy_threshold_input.setStyleSheet("font-size: 12px; color: black;")
        self.energy_threshold_input.setToolTip("Energy difference in kcal/mol between unique conformers for the second filter of E + RMS")
        advanced_layout.addWidget(self.energy_threshold_input, 4, 1)

        self.rms_threshold_label = QLabel("RMS threshold:", self) #(kcal/mol)?
        self.rms_threshold_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.rms_threshold_label, 0, 3)

        self.rms_threshold_input = QLineEdit(self)
        self.rms_threshold_input.setText("0.25")
        self.rms_threshold_input.setValidator(QDoubleValidator())
        self.rms_threshold_input.textChanged.connect(lambda: gen_command.update({"rms_threshold": self.rms_threshold_input.text()}))
        self.rms_threshold_input.setStyleSheet("font-size: 12px; color: black;")
        self.rms_threshold_input.setToolTip("RMS difference between unique conformers for the second filter of E + RMS")
        advanced_layout.addWidget(self.rms_threshold_input, 0, 4)

        self.opt_steps_rdkit_label = QLabel("RDKit opt steps:", self)
        self.opt_steps_rdkit_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.opt_steps_rdkit_label, 1, 3)

        self.opt_steps_rdkit_input = QSpinBox(self)
        self.opt_steps_rdkit_input.setRange(1, 10000)
        self.opt_steps_rdkit_input.setValue(1000)
        self.opt_steps_rdkit_input.valueChanged.connect(lambda: gen_command.update({"opt_steps_rdkit": self.opt_steps_rdkit_input.value()}))
        self.opt_steps_rdkit_input.setStyleSheet("font-size: 12px; color: black;")
        self.opt_steps_rdkit_input.setToolTip("Max cycles used in RDKit optimizations. ")
        advanced_layout.addWidget(self.opt_steps_rdkit_input, 1, 4)

        # Column 3 & 4
        self.max_matches_rmsd_label = QLabel("Max matches RMSD:", self)
        self.max_matches_rmsd_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.max_matches_rmsd_label, 2, 3)

        self.max_matches_rmsd_input = QSpinBox(self)
        self.max_matches_rmsd_input.setRange(1, 10000)
        self.max_matches_rmsd_input.setValue(1000)
        self.max_matches_rmsd_input.valueChanged.connect(lambda: gen_command.update({"max_matches_rmsd": self.max_matches_rmsd_input.value()}))
        self.max_matches_rmsd_input.setStyleSheet("font-size: 12px; color: black;")
        self.max_matches_rmsd_input.setToolTip("Max matches during RMS calculations for filtering \n(maxMatches option in the Chem.rdMolAlign.GetBestRMS() RDKit function)")
        advanced_layout.addWidget(self.max_matches_rmsd_input, 2, 4)

        self.max_mol_wt_label = QLabel("Max mol weight:", self)
        self.max_mol_wt_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.max_mol_wt_label, 3, 3)

        self.max_mol_wt_input = QSpinBox(self)
        self.max_mol_wt_input.setRange(0, 10000)
        self.max_mol_wt_input.setValue(0)
        self.max_mol_wt_input.valueChanged.connect(lambda: gen_command.update({"max_mol_wt": self.max_mol_wt_input.value()}))
        self.max_mol_wt_input.setSuffix(" g/mol")
        self.max_mol_wt_input.setStyleSheet("font-size: 12px; color: black;")
        self.max_mol_wt_input.setToolTip("Discard systems with molecular weights higher than this parameter (in g/mol). \nIf 0 is set, this filter is off.")
        advanced_layout.addWidget(self.max_mol_wt_input, 3, 4)

        self.max_torsions_label = QLabel("Max torsions:", self)
        self.max_torsions_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.max_torsions_label, 4, 3)

        self.max_torsions_input = QSpinBox(self)
        self.max_torsions_input.setRange(0, 1000)
        self.max_torsions_input.setValue(0)
        self.max_torsions_input.valueChanged.connect(lambda: gen_command.update({"max_torsions": self.max_torsions_input.value()}))
        self.max_torsions_input.setStyleSheet("font-size: 12px; color: black;")
        self.max_torsions_input.setToolTip("Discard systems with more than this many torsions (relevant to avoid molecules with many rotatable bonds). \nIf 0 is set, this filter is off.")
        advanced_layout.addWidget(self.max_torsions_input, 4, 4)

        # Column  5 & 6
        self.heavyonly_checkbox = QCheckBox("Heavy Only", self)
        self.heavyonly_checkbox.setChecked(True)
        self.heavyonly_checkbox.stateChanged.connect(lambda: gen_command.update({"heavyonly": self.heavyonly_checkbox.isChecked()}))
        self.heavyonly_checkbox.setStyleSheet("font-size: 12px; color: black;")
        self.heavyonly_checkbox.setToolTip("Only consider heavy atoms during RMS calculations for filtering \n(in the Chem.rdMolAlign.GetBestRMS() RDKit function)")
        advanced_layout.addWidget(self.heavyonly_checkbox, 0, 6)

        self.auto_metal_atoms_checkbox = QCheckBox("Auto Metal Atoms", self)
        self.auto_metal_atoms_checkbox.setChecked(True)
        self.auto_metal_atoms_checkbox.stateChanged.connect(lambda: gen_command.update({"auto_metal_atoms": self.auto_metal_atoms_checkbox.isChecked()}))
        self.auto_metal_atoms_checkbox.setStyleSheet("font-size: 12px; color: black;")
        self.auto_metal_atoms_checkbox.setToolTip("Automatically detect metal atoms for the RDKit conformer generation. \nCharge and mult should be specified as well since the automatic charge and mult detection might not be precise.")
        advanced_layout.addWidget(self.auto_metal_atoms_checkbox, 0, 5)

        self.seed_label = QLabel("Seed:", self)
        self.seed_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.seed_label, 1, 5)

        self.seed_input = QLineEdit(self)
        self.seed_input.setValidator(QIntValidator())  
        self.seed_input.setText("62609")
        self.seed_input.textChanged.connect(lambda: gen_command.update({"seed": self.seed_input.text()}))
        self.seed_input.setStyleSheet("font-size: 12px; color: black;")
        self.seed_input.setToolTip("Random seed used during RDKit embedding \n(in the Chem.rdDistGeom.EmbedMultipleConfs() RDKit function)")
        advanced_layout.addWidget(self.seed_input, 1, 6)

        self.bond_thres_label = QLabel("Bond threshold:", self)
        self.bond_thres_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.bond_thres_label, 2, 5)

        self.bond_thres_input = QDoubleSpinBox(self)
        self.bond_thres_input.setRange(0.0, 10.0)
        self.bond_thres_input.setValue(0.2)
        self.bond_thres_input.valueChanged.connect(lambda: gen_command.update({"bond_thres": self.bond_thres_input.value()}))
        self.bond_thres_input.setStyleSheet("font-size: 12px; color: black;")
        self.bond_thres_input.setToolTip("Threshold used to discard bonds in the geom option (+-0.2 A)")
        advanced_layout.addWidget(self.bond_thres_input, 2, 6)

        self.angle_thres_label = QLabel("Angle threshold:", self)
        self.angle_thres_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.angle_thres_label, 3, 5)

        self.angle_thres_input = QDoubleSpinBox(self)
        self.angle_thres_input.setRange(0.0, 360.0)
        self.angle_thres_input.setValue(30.0)
        self.angle_thres_input.valueChanged.connect(lambda: gen_command.update({"angle_thres": self.angle_thres_input.value()}))
        self.angle_thres_input.setStyleSheet("font-size: 12px; color: black;")
        self.angle_thres_input.setToolTip("Threshold used to discard angles in the geom option (+-30 degrees)")
        advanced_layout.addWidget(self.angle_thres_input, 3, 6)

        self.dihedral_thres_label = QLabel("Dihedral threshold:", self)
        self.dihedral_thres_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.dihedral_thres_label, 4, 5)

        self.dihedral_thres_input = QDoubleSpinBox(self)
        self.dihedral_thres_input.setRange(0.0, 360.0)
        self.dihedral_thres_input.setValue(30.0) 
        self.dihedral_thres_input.valueChanged.connect(lambda: gen_command.update({"dihedral_thres": self.dihedral_thres_input.value()}))
        self.dihedral_thres_input.setStyleSheet("font-size: 12px; color: black;")
        self.dihedral_thres_input.setToolTip("Threshold used to discard dihedrals in the geom option (+-30 degrees)")
        advanced_layout.addWidget(self.dihedral_thres_input, 4, 6)

        self.crest_force_label = QLabel("CREST force:", self)
        self.crest_force_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.crest_force_label, 0, 7)

        self.crest_force_input = QDoubleSpinBox(self)
        self.crest_force_input.setRange(0.0, 10.0)
        self.crest_force_input.setValue(0.5)
        self.crest_force_input.setSingleStep(0.1)
        self.crest_force_input.valueChanged.connect(lambda: crest_command.update({"crest_force": self.crest_force_input.value()}))
        self.crest_force_input.setStyleSheet("font-size: 12px; color: black;")
        self.crest_force_input.setToolTip("CREST ONLY: Force constant for constraints in the .xcontrol.sample file for CREST jobs.")
        advanced_layout.addWidget(self.crest_force_input, 0, 8)

        self.crest_keywords_input = QLineEdit(self)
        self.crest_keywords_input.setPlaceholderText("CREST keywords...")
        self.crest_keywords_input.textChanged.connect(lambda: crest_command.update({"crest_keywords": self.crest_keywords_input.text()}))
        self.crest_keywords_input.setStyleSheet("font-size: 12px; color: black;")
        self.crest_keywords_input.setToolTip("CREST ONLY: KDefine additional keywords to use in CREST that are not included in --chrg, --uhf, -T and -cinp. For example: '--alpb ch2cl2 --nci --cbonds 0.5'.")
        advanced_layout.addWidget(self.crest_keywords_input, 1, 7)

        self.xtb_keywords_input = QLineEdit(self)
        self.xtb_keywords_input.setPlaceholderText("xTB keywords...")
        self.xtb_keywords_input.textChanged.connect(lambda: crest_command.update({"xtb_keywords": self.xtb_keywords_input.text()}))
        self.xtb_keywords_input.setStyleSheet("font-size: 12px; color: black;")
        self.xtb_keywords_input.setToolTip("CREST ONLY: Define additional keywords to use in the xTB pre-optimization that are not included in -c, --uhf, -P and --input. For example: '--alpb ch2cl2 --gfn 1'.")
        advanced_layout.addWidget(self.xtb_keywords_input, 1, 8)

        self.cregen_checkbox = QCheckBox("CREGEN", self)
        self.cregen_checkbox.setChecked(True)
        self.cregen_checkbox.stateChanged.connect(lambda: crest_command.update({"cregen": self.cregen_checkbox.isChecked()}))
        self.cregen_checkbox.setStyleSheet("font-size: 12px; color: black;")
        self.cregen_checkbox.setToolTip("If True, perform a CREGEN analysis after CREST.")
        advanced_layout.addWidget(self.cregen_checkbox, 2, 7)

        self.cregen_keywords_input = QLineEdit(self)
        self.cregen_keywords_input.setPlaceholderText("CREGEN keywords...")
        self.cregen_keywords_input.textChanged.connect(lambda: crest_command.update({"cregen_keywords": self.cregen_keywords_input.text()}))
        self.cregen_keywords_input.setStyleSheet("font-size: 12px; color: black;")
        self.cregen_keywords_input.setToolTip("Additional keywords for CREGEN (i.e. cregen_keywords='--ethr 0.02').")
        advanced_layout.addWidget(self.cregen_keywords_input, 2, 8)

        self.crest_runs_label = QLabel("CREST runs:", self)
        self.crest_runs_label.setStyleSheet("font-size: 12px; color: black;")
        advanced_layout.addWidget(self.crest_runs_label, 3, 7)

        self.crest_runs_input = QSpinBox(self)
        self.crest_runs_input.setRange(1, 100)
        self.crest_runs_input.setValue(1)
        self.crest_runs_input.setSingleStep(1)
        self.crest_runs_input.valueChanged.connect(lambda: crest_command.update({"crest_runs": self.crest_runs_input.value()}))
        self.crest_runs_input.setStyleSheet("font-size: 12px; color: black;")
        self.crest_runs_input.setToolTip("Specify as number of runs if multiple starting points from RDKit starting points is required.")
        advanced_layout.addWidget(self.crest_runs_input, 3, 8)


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


    def handle_property_change(self, item):
        """Handle changes to editable properties and update the csv_dictionary."""
        logging.debug("at handle_property_change >>> handling property change")
        if item.row() == 0:  
            csv_model["code_name"][control.current_index - 1] = item.text()
        elif item.row() == 1:  
            csv_model["charge"][control.current_index - 1] = item.text()
            # Track user-defined charge
            if not hasattr(self, 'user_defined_charge'):
                self.user_defined_charge = {}
            self.user_defined_charge[control.current_index] = item.text()
        elif item.row() == 2:  
            csv_model["multiplicity"][control.current_index - 1] = item.text()
            # Track user-defined multiplicity
            if not hasattr(self, 'user_defined_multiplicity'):
                self.user_defined_multiplicity = {}
            self.user_defined_multiplicity[control.current_index] = item.text()


# SMILES HANDILING FUNCTIONS ??
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

        current_text = csv_model["SMILES"][control.current_index - 1]
        new_text = f"{current_text}.{smiles}" if current_text else smiles
        self.smiles_input.setText(new_text)

        if not csv_model["code_name"][control.current_index - 1]:
            csv_model["code_name"][control.current_index - 1] = code_name
        self.update_properties()
        self.update_ui()
        self.search_pubchem_input.clear()

    def metal_atom_detected(self, metals):
        """Detects metal atoms in the current molecule and updates the UI accordingly."""
        self.log_box_label.setText(f"Metal atoms detected: {', '.join(metals)}. Check the charge and multiplicity.")
        self.log_box_label.setStyleSheet("background-color: rgba(255, 255, 255, 0); border: none; font-size: 10px; color: black;")
        

    def smiles_are_bad_bro(self, smiles):
        """Got the signal the smiles are bad so gotta let the people kno ig"""
        self.log_box_label.setText(f"Bad SMILES: {smiles}")
        self.log_box_label.setStyleSheet("color: red; font-size: 8px; border: none; background-color: white;")
        self.log_box_label.setWordWrap(True)


    def resizeEvent(self, event):
        """Handle window resize events to refresh the molecule display."""
        super().resizeEvent(event)
        control.display_molecule(self.show_numbered_atoms_toggle.isChecked())

    def mousePressEvent(self, event: QMouseEvent):
        if event.button() == Qt.MouseButton.LeftButton:
            pos = event.position()
            control.mousePressEvent(pos)

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
            ("code_name", csv_model["code_name"][control.current_index - 1]),
            ("charge", csv_model["charge"][control.current_index - 1]),
            ("multiplicity", csv_model["multiplicity"][control.current_index - 1]),
            ("constraints_atoms", csv_model["constraints_atoms"][control.current_index - 1]),
            ("constraints_dist", csv_model["constraints_dist"][control.current_index - 1]),
            ("constraints_angle", csv_model["constraints_angle"][control.current_index - 1]),
            ("constraints_dihedral", csv_model["constraints_dihedral"][control.current_index - 1]),
            ("complex_type", csv_model["complex_type"][control.current_index - 1]),
            ("geom", csv_model["geom"][control.current_index - 1]),
            ]
            
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
        finally:
            self.is_programmatic_update = old_value





# UI ELEMENTS UPDATE FUNCTIONS 
    def update_ui(self):
        """Update all UI elements to reflect the current molecule's data."""
        smiles = csv_model["SMILES"][control.current_index - 1]
        self.smiles_input.blockSignals(True)
        self.smiles_input.setText(smiles)
        self.smiles_output.setText(smiles2enumerate(smiles))

        cursor = self.smiles_input.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        self.smiles_input.setTextCursor(cursor)

        self.smiles_input.blockSignals(False)

        control.display_molecule(self.show_numbered_atoms_toggle.isChecked())

        self.index_and_total_label_update()
        self.log_box_label.clear()

        try:
            num_atoms = smiles2numatoms(smiles)
            num_electrons = smiles2numelectrons(smiles)
            self.atom_electron_label.setText(f" Electrons: {num_electrons}\n Atoms: {num_atoms}")
        except Exception as e:
            logging.error(f"Error calculating atom/electron counts: {e}")
            self.atom_electron_label.setText(" Electrons: 0\n Atoms: 0")

        self.update_properties()
        
    def index_and_total_label_update(self):
        control.total_index = control.get_total_index()    
        self.index_and_total_label.setText(f"{control.current_index}/{control.total_index}")

    def clear_focus_on_inputs(self):
        self.smiles_input.clearFocus()
        self.search_pubchem_input.clearFocus()

    def closeEvent(self, event):
        """Handle the close event to prompt saving the CSV file, with the added icon."""
        pixmap = QPixmap(icons.red_path)
        icon = QIcon(pixmap)
        
        if gen_command["input"] is None: 
            msgBox = QMessageBox(self)
            msgBox.setWindowTitle("Save CSV")
            msgBox.setText("Would you like to save the CSV file before exiting?")
            msgBox.setWindowIcon(icon)
            msgBox.setIconPixmap(pixmap)
            msgBox.setStandardButtons(QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel)
            reply = msgBox.exec()
            
            if reply == QMessageBox.Save:
                if not control.save_csv_file():
                    event.ignore()  # Prevent closing if saving fails
                else:
                    event.accept()
            elif reply == QMessageBox.Discard:
                event.accept()
            else:
                event.ignore()
        else:
            event.accept()

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
                # for key in csv_model.keys():
                    # csv_model[key].pop(0) 
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
                    
                control.total_index = control.get_total_index()
                if control.total_index > 0:
                    control.current_index = 1  # assuming indices start from 1
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
            control.total_index = control.get_total_index()

            for index in range(control.total_index):
                csv_model["code_name"][index] = f"mol_{index + 1}"
                control.current_index = 1
            gen_command["input"] = None
            csv_model.signals.updated.emit()
            self.update_properties()

        except ImportError as e:
            print(f"Error importing required module: {e}")
        except IOError as e:
            print(f"Error reading or writing file: {e}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

# AQME RUN SETUP FUNCTIONS
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

        if gen_command["input"] is None or not gen_command["input"]:
            msgBox = QMessageBox(self)
            msgBox.setIconPixmap(QPixmap(icons.red_path))
            msgBox.setWindowTitle("Warning")
            msgBox.setText("Please save the input file before running AQME.")
            msgBox.setStandardButtons(QMessageBox.StandardButton.Ok)
            msgBox.exec()
            control.save_csv_file() 
            return
        command = self.aqme_rungen()
        file_directory = os.path.dirname(gen_command["input"])

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
        input_file, program, destination, stacksize, sample, auto_sample, ewin_csearch, \
        initial_energy_threshold, energy_threshold, rms_threshold, \
        opt_steps_rdkit, heavyonly, max_matches_rmsd, max_mol_wt, \
        max_torsions, seed, geom, bond_thres, angle_thres, dihedral_thres, auto_metal_atoms = (
            gen_command[key] for key in [
            "input","program", "destination", "stacksize", "sample", "auto_sample",
            "ewin_csearch", "initial_energy_threshold", "energy_threshold",
            "rms_threshold", "opt_steps_rdkit", "heavyonly",
            "max_matches_rmsd", "max_mol_wt", "max_torsions", "seed", "geom",
            "bond_thres", "angle_thres", "dihedral_thres", "auto_metal_atoms"
            ]
        )
        if not input_file:
            self.failure("Please save the file before running AQME.")
            return 
        if not os.path.exists(input_file):
            self.failure("Input file does not exist.")
            return 
        if not destination:
            destination = os.path.dirname(input_file)
        if not os.path.exists(destination):
            os.makedirs(destination)
        if not os.path.exists(destination):
            self.failure("Output directory does not exist.")
            return 
        
        aqme_rungen = f'python -u -m aqme --csearch --program {program} --input {input_file}  --destination {destination} --stacksize {stacksize} --sample {sample} --auto_sample {auto_sample} --ewin_csearch {ewin_csearch} --initial_energy_threshold {initial_energy_threshold} --energy_threshold {energy_threshold} --rms_threshold {rms_threshold} --opt_steps_rdkit {opt_steps_rdkit} --heavyonly {heavyonly} --max_matches_rmsd {max_matches_rmsd} --max_mol_wt {max_mol_wt} --max_torsions {max_torsions} --seed {seed} --geom "{geom}" --bond_thres {bond_thres} --angle_thres {angle_thres} --dihedral_thres {dihedral_thres} --auto_metal_atoms {auto_metal_atoms}'

        if program == "CREST":
            nprocs, crest_force, crest_keywords, cregen, cregen_keywords, xtb_keywords, crest_runs = (
                crest_command[key] for key in [
                "nprocs", "crest_force", "crest_keywords", "cregen", "cregen_keywords",
                "xtb_keywords", "crest_runs"
                ]
            )
            aqme_rungen += f' --nprocs {nprocs} --crest_force {crest_force}--cregen {cregen} --crest_runs {crest_runs}'

            if crest_keywords != None:
                aqme_rungen += f' --crest_keywords "{crest_keywords}"'
            if cregen_keywords != None:
                aqme_rungen += f' --cregen_keywords "{cregen_keywords}"'
            if xtb_keywords != None:
                aqme_rungen += f' --xtb_keywords "{xtb_keywords}"'
        return aqme_rungen

# random 
    def success(self, message):
        pixmap = QPixmap(icons.green_path)
        icon = QIcon(pixmap)
        msg = QMessageBox(self)
        msg.setWindowTitle("Success")
        msg.setText(message)
        msg.setWindowIcon(icon)
        msg.setIconPixmap(pixmap)
        msg.exec()

    def failure(self, message):
        pixmap = QPixmap(icons.red_path)
        icon = QIcon(pixmap)
        msg = QMessageBox(self)
        msg.setWindowTitle("Failure")
        msg.setText(message)
        msg.setWindowIcon(icon)
        msg.setIconPixmap(pixmap)
        msg.exec()

    def eventFilter(self, obj, event):
        # Check if the event is for the QTextEdit and is a key press
        if obj == self.smiles_input and event.type() == QEvent.KeyPress:
            # Ignore Space and Enter keys
            if event.key() in (Qt.Key_Space, Qt.Key_Return, Qt.Key_Enter):
                return True  # Block the event
        return super().eventFilter(obj, event) 

# TODO:
# - Add the complex_type option when TM is found !
# - Fix run aqme command  !!!!
# - expand the pubchem search  ... !
# ---- fix the warning where the smiles entered are not valid (rn nothing really happens but the change is only regirestered when the smiles are valid (or thats how it should be)), similarly changing smiles removes the view of the constraints but doesnt remove them completely at the index level, which I think should be done?
# add intermedaite plus transition state option in the show all 

# fucking dark mode man !!!!
