import logging

from PySide6.QtWidgets import  QLabel, QVBoxLayout, QWidget, QPushButton, QHBoxLayout, QTextBrowser,QLineEdit, QPlainTextEdit, QCheckBox, QMessageBox, QSizePolicy, QFileDialog, QTableWidget, QTableWidgetItem, QHeaderView, QApplication, QComboBox, QSpinBox, QStyle, QTableWidgetItem, QFrame, QGridLayout, QDoubleSpinBox, QGroupBox, QProgressBar
from PySide6.QtCore import Qt, Slot
from PySide6.QtGui import QPixmap, QKeySequence, QShortcut, QMouseEvent, QIcon, QDoubleValidator, QTextCursor, QIntValidator    

from aqmeasy.utils import  pubchem2smiles, smiles2enumerate, smiles2numatoms, smiles2numelectrons, smiles2findmetal, smiles2ismetalcomplex,command2clipboard, load_svg_as_pixmap 

from aqmeasy.controllers.CSEARCH_controller import CsvController
from aqmeasy.ui.stylesheets import stylesheets
from aqmeasy.ui.icons import Icons


class CSEARCHWidget(QWidget):
    def __init__(self, parent, model, general_command_model):
        super().__init__()
        self.parent = parent
        self.csv_model = model
        self.control = CsvController(self.csv_model, general_command_model)
        self.control.set_parent(self)

        self.csv_model.signals.updated.connect(self.update_ui)
        self.metal_detected_in_smiles = False

        self.setAcceptDrops(True)
        self.dragEnterEvent = self._drag_enter
        self.dragLeaveEvent = self._drag_leave 
        self.dropEvent = self._drop_event 

        self.file_name = None # this actually not needed but i need to fix the logic when nothing has been changed in the UI and we want to close the app (i.e. no need to save)
        
        QShortcut(QKeySequence(Qt.Key.Key_Escape), self, lambda: self.clear_focus_on_inputs())

        self.top_layout = QHBoxLayout()
        self.left_layout = QVBoxLayout()
        self.right_layout = QVBoxLayout()

        self.control1_layout = QHBoxLayout()
        self.control2_layout = QHBoxLayout()    
        self.control3_layout = QHBoxLayout()

        self.middle_layout = QHBoxLayout()

        self.right_layout.addLayout(self.control1_layout)
        self.bottom_layout = QHBoxLayout()
        self.main_layout = QVBoxLayout()

        self.top_layout.addLayout(self.left_layout, 3)  
        self.top_layout.addLayout(self.right_layout, 1)  
        self.main_layout.addLayout(self.top_layout, 3)
        self.main_layout.addLayout(self.bottom_layout, 1)
        self.main_layout.addLayout(self.middle_layout)
        self.setLayout(self.main_layout)

        self.index_and_total_label = QLabel(f"{self.control.current_index}/{self.control.total_index}")
        self.index_and_total_label.setFixedWidth(60)
        self.control1_layout.addWidget(self.index_and_total_label)

        self.import_button = QPushButton("Import", self)
        self.import_button.setIcon(QIcon(load_svg_as_pixmap(Icons.add_file)))
        self.import_button.clicked.connect(lambda: (logging.debug("at import_button >>> self.control.import_file()"), self.control.import_file()))
        self.control1_layout.addWidget(self.import_button)

        self.show_all_button = QPushButton("Show All", self)
        self.show_all_button.setIcon(QIcon(load_svg_as_pixmap(Icons.external_link)))
        self.show_all_button.clicked.connect(self.control.show_csv)
        self.control1_layout.addWidget(self.show_all_button)

        self.save_csv_button = QPushButton(self)
        self.save_csv_button.setIcon(QIcon(load_svg_as_pixmap(Icons.save)))
        self.save_csv_button.setToolTip("Save CSV file")
        self.save_csv_button.clicked.connect(lambda: (logging.debug("at save_csv_button >>> save_csv"), self.control.save_csv_file()))
        self.control1_layout.addWidget(self.save_csv_button)

    # THIS IS WHERE THE MOLECULE IS DISPLAYED
        molecule_group = QGroupBox("Click to select atoms and add constraints, click again to deselect.")

        self.molecule_label = QLabel()
        self.molecule_label.setStyleSheet(stylesheets.MoleculeLabel)
        self.molecule_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.molecule_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)  
        self.molecule_label.setMinimumSize(250, 200)  
        
        molecule_layout = QVBoxLayout(molecule_group)
        molecule_layout.addWidget(self.molecule_label)
        self.molecule_label.setCursor(Qt.CursorShape.PointingHandCursor)
        self.molecule_label.mousePressEvent = self._molecule_label_mouse_press
        self.left_layout.addWidget(molecule_group)

        self.atom_electron_label = QLabel(self.molecule_label)
        self.atom_electron_label.setStyleSheet("background-color: rgba(255, 255, 255, 0); border: none; font-size: 10px; color: black;")
        self.atom_electron_label.setFixedSize(150, 30)

        self.log_box_label = QLabel(self.molecule_label)
        self.log_box_label.setStyleSheet("background-color: rgba(255, 255, 255, 0); border: none; font-size: 10px; color: black;")        
        self.log_box_label.setFixedSize(550, 20)
        self.molecule_label.resizeEvent = lambda event: self.log_box_label.move(5, self.molecule_label.height() - self.log_box_label.height())

    # VARIOUS INPUTS
        self.smiles_input = QPlainTextEdit(self)
        self.smiles_input.setPlaceholderText("Enter SMILES, search PubChem below, drop in a ChemDraw/CSV/SDF file or click on the icon below for a SMILES tutorial...")
        self.smiles_input.setAutoFillBackground(True)
        self.smiles_input.setAcceptDrops(False)
        self.right_layout.addWidget(self.smiles_input, 1)

        self.smiles_input.textChanged.connect(lambda: self.control.update_smiles_model(self.smiles_input.toPlainText()))

        self.show_numbered_atoms_toggle = QCheckBox("Show atom labels", self)
        self.show_numbered_atoms_toggle.setChecked(False)
        self.show_numbered_atoms_toggle.stateChanged.connect(lambda: self.control.display_molecule(self.show_numbered_atoms_toggle.isChecked()))
        self.control2_layout.addWidget(self.show_numbered_atoms_toggle)

        self.search_pubchem_input = QLineEdit(self)
        self.search_pubchem_input.setPlaceholderText("Search PubChem...")
        self.search_pubchem_input.returnPressed.connect(self.smiles_from_pubchem)
        self.control2_layout.addWidget(self.search_pubchem_input)

        self.smiles_tutorial_button = QPushButton(self)
        self.smiles_tutorial_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_FileDialogDetailedView))
        self.smiles_tutorial_button.setToolTip("SMILES Tutorial")
        self.smiles_tutorial_button.clicked.connect(self.control.open_smiles_tutorial)
        self.control2_layout.addWidget(self.smiles_tutorial_button)

        self.right_layout.addLayout(self.control2_layout)
        self.right_layout.addLayout(self.control3_layout)

    # CONTROL BUTTONS
        self.delete_button = QPushButton()
        self.delete_button.setIcon(QIcon(load_svg_as_pixmap(Icons.trash)))
        self.delete_button.clicked.connect(self.control.delete_molecule)
        self.control3_layout.addWidget(self.delete_button)

        self.new_molecule_button = QPushButton()
        self.new_molecule_button.setIcon(QIcon(load_svg_as_pixmap(Icons.plus)))
        self.new_molecule_button.clicked.connect(lambda: (logging.debug("at new_molecule_button >>> self.new_molecule()"), self.control.new_molecule()))
        self.control3_layout.addWidget(self.new_molecule_button)

        self.previous_button = QPushButton()
        self.previous_button.setIcon(QIcon(load_svg_as_pixmap(Icons.chevron_double_left)))
        self.previous_button.clicked.connect(lambda: (logging.debug("at previous_button >>> self.previous_molecule"), self.control.previous_molecule()))
        QShortcut(QKeySequence(Qt.Key.Key_Left), self, self.control.previous_molecule)
        self.control3_layout.addWidget(self.previous_button)

        self.next_button = QPushButton()
        self.next_button.setIcon(QIcon(load_svg_as_pixmap(Icons.chevron_double_right)))
        self.next_button.setLayoutDirection(Qt.LayoutDirection.RightToLeft)
        self.next_button.clicked.connect(lambda: (logging.debug("at next_button >>> self.next_molecule"), self.control.next_molecule()))
        QShortcut(QKeySequence(Qt.Key.Key_Right), self, self.control.next_molecule)
        self.control3_layout.addWidget(self.next_button)

    # OUTPUTS
        self.smiles_output = QTextBrowser(self)
        self.smiles_output.setStyleSheet(stylesheets.QTextBrowser)
        self.smiles_output.setAutoFillBackground(True)
        self.smiles_output.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.right_layout.addWidget(self.smiles_output, 1)

        self.properties_table = QTableWidget(self)
        self.properties_table.setRowCount(9)
        self.properties_table.setColumnCount(2)
        self.properties_table.horizontalHeader().setVisible(False)
        self.properties_table.verticalHeader().setVisible(False)
        self.properties_table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        complex_type_combo = QComboBox()
        complex_type_combo.addItems(["", "squareplanar", "squarepyramidal", "linear", "trigonalplanar"])
        complex_type_combo.currentTextChanged.connect(lambda text: self.handle_combobox_change(8, text))
        self.properties_table.setCellWidget(8,1,complex_type_combo)
        self.properties_table.cellWidget(8,1).setEnabled(False)


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

        self.right_layout.addWidget(self.properties_table, 2)

        # Progress bar in middle layout
        self.progress_bar = QProgressBar(self)
        self.progress_bar.setFixedHeight(10)
        self.middle_layout.addWidget(self.progress_bar)


        shell_group = QGroupBox("Shell Output")

        self.shell_output = QTextBrowser(self)
        self.shell_output.setStyleSheet(stylesheets.ShellOutput)
        self.shell_output.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.shell_output.setReadOnly(True)
        shell_layout = QVBoxLayout(shell_group)
        shell_layout.addWidget(self.shell_output)

        self.bottom_layout.addWidget(shell_group, 2)

        # AQME SETUP
        self.aqme_setup_grid = QGridLayout()
        self.bottom_layout.addLayout(self.aqme_setup_grid, 1)

        # Row 0 - Program selection
        self.program_label = QLabel("Select CSEARCH program:", self)
        self.program_label.setStyleSheet(stylesheets.QLabel)
        self.aqme_setup_grid.addWidget(self.program_label, 0, 0)

        self.program_combo = QComboBox(self)
        self.program_combo.addItems(["RDKit", "CREST", "GOAT*"])
        self.program_combo.model().item(2).setEnabled(False)  # type: ignore
        self.control.update_command("program", "rdkit") # set default program
        self.program_combo.currentTextChanged.connect(lambda text: self.control.update_command("program", text.lower()))
        self.aqme_setup_grid.addWidget(self.program_combo, 0, 1)

        # Row 1 - Number of processors
        self.nprocs_label = QLabel("Number of processors:", self)
        self.aqme_setup_grid.addWidget(self.nprocs_label, 1, 0)
        
        self.nprocs_input = QSpinBox(self)
        self.nprocs_input.setRange(1, 40)
        self.nprocs_input.setValue(8)
        self.nprocs_input.valueChanged.connect(lambda: self.control.update_command("nprocs", self.nprocs_input.value()))
        self.nprocs_input.setToolTip("Number of processors to use for the calculation.\nOnly relevant for CREST and GOAT*.")
        self.nprocs_input.setEnabled(self.program_combo.currentText() == "CREST")
        self.aqme_setup_grid.addWidget(self.nprocs_input, 1, 1)
        self.program_combo.currentTextChanged.connect(lambda text: self.nprocs_input.setEnabled(text == "CREST"))

        # Row 2 - Stack size
        self.stacksize_label = QLabel("Stack size (GB):", self)
        self.aqme_setup_grid.addWidget(self.stacksize_label, 2, 0)
        self.stacksize_input = QSpinBox(self)
        self.stacksize_input.setRange(1, 8)
        self.stacksize_input.setValue(1)
        self.stacksize_input.valueChanged.connect(lambda: self.control.update_command("stacksize", f"{self.stacksize_input.value()}GB"))
        self.aqme_setup_grid.addWidget(self.stacksize_input, 2, 1)

        # Row 3 - Output directory
        self.output_dir_label = QLabel("Output directory:", self)
        self.aqme_setup_grid.addWidget(self.output_dir_label, 3, 0)
        self.output_dir_input = QLineEdit(self)
        self.output_dir_input.setPlaceholderText("Select output directory...")
        self.output_dir_input.textChanged.connect(lambda: self.control.update_command("destination", self.output_dir_input.text()))

        self.output_dir_button = QPushButton(self)

        self.output_dir_button.setIcon(QIcon(load_svg_as_pixmap(Icons.folder_open)))
        self.output_dir_button.clicked.connect(self.select_output_directory)
        output_dir_layout = QHBoxLayout()
        output_dir_layout.addWidget(self.output_dir_input)
        output_dir_layout.addWidget(self.output_dir_button)
        self.aqme_setup_grid.addLayout(output_dir_layout, 3, 1)

        # Row 4 - copy command/save csv
        self.copy_command_button = QPushButton("Copy Command", self)
        self.copy_command_button.setIcon(QIcon(load_svg_as_pixmap(Icons.command_line)))
        self.copy_command_button.clicked.connect(lambda: self.success("Command copied to clipboard.") if command2clipboard(self.control.generate_csearch_command()) == True else self.failure("Failed to copy command to clipboard."))
        self.aqme_setup_grid.addWidget(self.copy_command_button, 4,0)

        # Row 5 - Run button
        self.run_button = QPushButton("Run CSEARCH", self)
        self.run_button.setStyleSheet(stylesheets.RunButton)
        self.run_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay))

        self.run_button.clicked.connect(self.toggle_csearch)

        self.parent.worker.error.connect(self.failure)
        self.parent.worker.result.connect(self.shell_output.append)
        self.parent.worker.finished_signal.connect(self.success)
        self.parent.worker.confirm.connect(lambda msg: self.handle_qprep_confirm(msg, self.parent.worker.model["destination"]))
        self.parent.worker.finished.connect(self.on_thread_complete)

        self.aqme_setup_grid.addWidget(self.run_button, 4, 1, 2, 1)

        for widget in [self.smiles_input, self.smiles_output, self.properties_table,  self.shell_output, self.molecule_label, self.atom_electron_label]:
            palette = widget.palette()
            palette.setColor(widget.backgroundRole(), self.palette().color(self.backgroundRole()))
            widget.setPalette(palette)

        # ADVANCED SETTINGS

        # Column 1 & 2
        self.advanced_settings_button = QPushButton("Advanced Settings", self)
        self.advanced_settings_button.setIcon(QIcon(Icons.eye))
        self.advanced_settings_button.setCheckable(True)
        self.advanced_settings_button.clicked.connect(lambda: self.toggle_panel(self.height(), self.width()))
        self.aqme_setup_grid.addWidget(self.advanced_settings_button, 5, 0)

        self.advanced_panel = QFrame()
        self.advanced_panel.setFixedHeight(0)
        self.main_layout.addWidget(self.advanced_panel)

        advanced_panel_layout = QVBoxLayout()
        self.advanced_panel.setLayout(advanced_panel_layout)
        advanced_settings_group = QGroupBox("Advanced Settings")
        advanced_panel_layout.addWidget(advanced_settings_group)
        advanced_layout = QGridLayout()
        advanced_settings_group.setLayout(advanced_layout)

        self.sample_size_label = QLabel("Sample size:", self)
        advanced_layout.addWidget(self.sample_size_label, 0, 0)

        self.sample_size_input = QSpinBox(self)
        self.sample_size_input.setRange(1, 500)
        self.sample_size_input.setValue(25)
        self.sample_size_input.valueChanged.connect(lambda: self.control.update_command("sample", self.sample_size_input.value()))
        self.sample_size_input.setToolTip("Number of conformers to keep after the initial RDKit sampling.\nThey are selected using a combination of RDKit energies and Butina clustering.")
        advanced_layout.addWidget(self.sample_size_input, 0, 1)

        self.auto_sample_label = QLabel("Auto sample level:", self)
        advanced_layout.addWidget(self.auto_sample_label, 1, 0)
        self.auto_sample_combo = QComboBox(self)
        self.auto_sample_combo.addItems(["low", "mid", "high", "false"])
        self.auto_sample_combo.setCurrentText("mid")
        self.auto_sample_combo.currentTextChanged.connect(lambda text: self.control.update_command("auto_sample", text))
        self.auto_sample_combo.setToolTip("Apply automatic calculation of the number of conformers generated initially with RDKit. \nThis number of conformers is initially generated and then reduced to the number specified in --sample with different filters. \nOptions:\n• Low: Base multiplier = 5, max confs = 100\n• Mid: Base multiplier = 10, max confs = 250\n• High: Base multiplier = 20, max confs = 500\n• False: Use the number specified in --sample")
        advanced_layout.addWidget(self.auto_sample_combo, 1, 1)

        self.energy_window_label = QLabel("Energy window:", self)
        advanced_layout.addWidget(self.energy_window_label, 2, 0)
        self.energy_window_input = QDoubleSpinBox(self)
        self.energy_window_input.setRange(1.0, 50.0)
        self.energy_window_input.setValue(5.0)
        self.energy_window_input.valueChanged.connect(lambda: self.control.update_command("ewin_csearch", self.energy_window_input.value()))
        self.energy_window_input.setToolTip("Energy window in kcal/mol to discard conformers\n(i.e. if a conformer is more than the E window compared to the most stable conformer).")
        advanced_layout.addWidget(self.energy_window_input, 2, 1)

        self.initial_energy_threshold_label = QLabel("Initial E threshold:", self)
        advanced_layout.addWidget(self.initial_energy_threshold_label, 3, 0)
        self.initial_energy_threshold_input = QLineEdit(self)
        self.initial_energy_threshold_input.setText("0.0001")
        self.initial_energy_threshold_input.setValidator(QDoubleValidator())
        self.initial_energy_threshold_input.textChanged.connect(lambda: self.control.update_command("initial_energy_threshold", self.initial_energy_threshold_input.text()))
        self.initial_energy_threshold_input.setToolTip("Energy difference in kcal/mol between unique conformers for the first filter of only E.")
        advanced_layout.addWidget(self.initial_energy_threshold_input, 3, 1)

        self.energy_threshold_label = QLabel("Energy threshold:", self)
        advanced_layout.addWidget(self.energy_threshold_label, 4, 0)

        self.energy_threshold_input = QLineEdit(self)
        self.energy_threshold_input.setText("0.25")
        self.energy_threshold_input.setValidator(QDoubleValidator())
        self.energy_threshold_input.textChanged.connect(lambda: self.control.update_command("energy_threshold", self.energy_threshold_input.text()))
        self.energy_threshold_input.setToolTip("Energy difference in kcal/mol between unique conformers for the second filter of E + RMS")
        advanced_layout.addWidget(self.energy_threshold_input, 4, 1)

        self.rms_threshold_label = QLabel("RMS threshold:", self) #(kcal/mol)?
        advanced_layout.addWidget(self.rms_threshold_label, 0, 3)

        self.rms_threshold_input = QLineEdit(self)
        self.rms_threshold_input.setText("0.25")
        self.rms_threshold_input.setValidator(QDoubleValidator())
        self.rms_threshold_input.textChanged.connect(lambda: self.control.update_command("rms_threshold", self.rms_threshold_input.text()))
        self.rms_threshold_input.setToolTip("RMS difference between unique conformers for the second filter of E + RMS")
        advanced_layout.addWidget(self.rms_threshold_input, 0, 4)

        self.opt_steps_rdkit_label = QLabel("RDKit opt steps:", self)
        advanced_layout.addWidget(self.opt_steps_rdkit_label, 1, 3)

        self.opt_steps_rdkit_input = QSpinBox(self)
        self.opt_steps_rdkit_input.setRange(1, 10000)
        self.opt_steps_rdkit_input.setValue(1000)
        self.opt_steps_rdkit_input.valueChanged.connect(lambda: self.control.update_command("opt_steps_rdkit", self.opt_steps_rdkit_input.value()))
        self.opt_steps_rdkit_input.setToolTip("Max cycles used in RDKit optimizations. ")
        advanced_layout.addWidget(self.opt_steps_rdkit_input, 1, 4)

        # Column 3 & 4
        self.max_matches_rmsd_label = QLabel("Max matches RMSD:", self)
        advanced_layout.addWidget(self.max_matches_rmsd_label, 2, 3)

        self.max_matches_rmsd_input = QSpinBox(self)
        self.max_matches_rmsd_input.setRange(1, 10000)
        self.max_matches_rmsd_input.setValue(1000)
        self.max_matches_rmsd_input.valueChanged.connect(lambda: self.control.update_command("max_matches_rmsd", self.max_matches_rmsd_input.value()))
        self.max_matches_rmsd_input.setToolTip("Max matches during RMS calculations for filtering \n(maxMatches option in the Chem.rdMolAlign.GetBestRMS() RDKit function)")
        advanced_layout.addWidget(self.max_matches_rmsd_input, 2, 4)

        self.max_mol_wt_label = QLabel("Max mol weight:", self)
        advanced_layout.addWidget(self.max_mol_wt_label, 3, 3)

        self.max_mol_wt_input = QSpinBox(self)
        self.max_mol_wt_input.setRange(0, 10000)
        self.max_mol_wt_input.setValue(0)
        self.max_mol_wt_input.valueChanged.connect(lambda: self.control.update_command("max_mol_wt", self.max_mol_wt_input.value()))
        self.max_mol_wt_input.setSuffix(" g/mol")
        self.max_mol_wt_input.setToolTip("Discard systems with molecular weights higher than this parameter (in g/mol). \nIf 0 is set, this filter is off.")
        advanced_layout.addWidget(self.max_mol_wt_input, 3, 4)

        self.max_torsions_label = QLabel("Max torsions:", self)
        advanced_layout.addWidget(self.max_torsions_label, 4, 3)

        self.max_torsions_input = QSpinBox(self)
        self.max_torsions_input.setRange(0, 1000)
        self.max_torsions_input.setValue(0)
        self.max_torsions_input.valueChanged.connect(lambda: self.control.update_command("max_torsions", self.max_torsions_input.value()))
        self.max_torsions_input.setToolTip("Discard systems with more than this many torsions (relevant to avoid molecules with many rotatable bonds). \nIf 0 is set, this filter is off.")
        advanced_layout.addWidget(self.max_torsions_input, 4, 4)

        # Column  5 & 6
        self.heavyonly_checkbox = QCheckBox("Heavy Only", self)
        self.heavyonly_checkbox.setChecked(True)
        self.heavyonly_checkbox.stateChanged.connect(lambda: self.control.update_command("heavyonly", self.heavyonly_checkbox.isChecked()))
        self.heavyonly_checkbox.setToolTip("Only consider heavy atoms during RMS calculations for filtering \n(in the Chem.rdMolAlign.GetBestRMS() RDKit function)")
        advanced_layout.addWidget(self.heavyonly_checkbox, 0, 6)

        self.auto_metal_atoms_checkbox = QCheckBox("Auto Metal Atoms", self)
        self.auto_metal_atoms_checkbox.setChecked(True)
        self.auto_metal_atoms_checkbox.stateChanged.connect(lambda: self.control.update_command("auto_metal_atoms", self.auto_metal_atoms_checkbox.isChecked()))
        self.auto_metal_atoms_checkbox.setToolTip("Automatically detect metal atoms for the RDKit conformer generation. \nCharge and mult should be specified as well since the automatic charge and mult detection might not be precise.")
        advanced_layout.addWidget(self.auto_metal_atoms_checkbox, 0, 5)

        self.seed_label = QLabel("Seed:", self)
        advanced_layout.addWidget(self.seed_label, 1, 5)

        self.seed_input = QLineEdit(self)
        self.seed_input.setValidator(QIntValidator())  
        self.seed_input.setText("62609")
        self.seed_input.textChanged.connect(lambda: self.control.update_command("seed", self.seed_input.text()))
        self.seed_input.setToolTip("Random seed used during RDKit embedding \n(in the Chem.rdDistGeom.EmbedMultipleConfs() RDKit function)")
        advanced_layout.addWidget(self.seed_input, 1, 6)

        self.bond_thres_label = QLabel("Bond threshold:", self)
        advanced_layout.addWidget(self.bond_thres_label, 2, 5)

        self.bond_thres_input = QDoubleSpinBox(self)
        self.bond_thres_input.setRange(0.0, 10.0)
        self.bond_thres_input.setValue(0.2)
        self.bond_thres_input.valueChanged.connect(lambda: self.control.update_command("bond_thres", self.bond_thres_input.value()))
        self.bond_thres_input.setToolTip("Threshold used to discard bonds in the geom option (+-0.2 A)")
        advanced_layout.addWidget(self.bond_thres_input, 2, 6)

        self.angle_thres_label = QLabel("Angle threshold:", self)
        advanced_layout.addWidget(self.angle_thres_label, 3, 5)

        self.angle_thres_input = QDoubleSpinBox(self)
        self.angle_thres_input.setRange(0.0, 360.0)
        self.angle_thres_input.setValue(30.0)
        self.angle_thres_input.valueChanged.connect(lambda: self.control.update_command("angle_thres", self.angle_thres_input.value()))
        self.angle_thres_input.setToolTip("Threshold used to discard angles in the geom option (+-30 degrees)")
        advanced_layout.addWidget(self.angle_thres_input, 3, 6)

        self.dihedral_thres_label = QLabel("Dihedral threshold:", self)
        advanced_layout.addWidget(self.dihedral_thres_label, 4, 5)

        self.dihedral_thres_input = QDoubleSpinBox(self)
        self.dihedral_thres_input.setRange(0.0, 360.0)
        self.dihedral_thres_input.setValue(30.0) 
        self.dihedral_thres_input.valueChanged.connect(lambda: self.control.update_command("dihedral_thres", self.dihedral_thres_input.value()))
        self.dihedral_thres_input.setToolTip("Threshold used to discard dihedrals in the geom option (+-30 degrees)")
        advanced_layout.addWidget(self.dihedral_thres_input, 4, 6)

        self.crest_force_label = QLabel("CREST force:", self)
        advanced_layout.addWidget(self.crest_force_label, 0, 7)

        self.crest_force_input = QDoubleSpinBox(self)

        self.crest_force_input.setRange(0.0, 10.0)
        self.crest_force_input.setValue(0.5)
        self.crest_force_input.setSingleStep(0.1)
        self.crest_force_input.valueChanged.connect(lambda: self.control.update_crest_command("crest_force", self.crest_force_input.value()))
        self.crest_force_input.setToolTip("CREST ONLY: Force constant for constraints in the .xcontrol.sample file for CREST jobs.")
        advanced_layout.addWidget(self.crest_force_input, 0, 8)

        self.crest_keywords_input = QLineEdit(self)
        self.crest_keywords_input.setPlaceholderText("CREST keywords...")
        self.crest_keywords_input.textChanged.connect(lambda: self.control.update_crest_command("crest_keywords", self.crest_keywords_input.text()))
        self.crest_keywords_input.setToolTip("CREST ONLY: Define additional keywords to use in CREST that are not included in --chrg, --uhf, -T and -cinp. For example: '--alpb ch2cl2 --nci --cbonds 0.5'.")
        advanced_layout.addWidget(self.crest_keywords_input, 1, 7)

        self.xtb_keywords_input = QLineEdit(self)
        self.xtb_keywords_input.setPlaceholderText("xTB keywords...")
        self.xtb_keywords_input.textChanged.connect(lambda: self.control.update_crest_command("xtb_keywords", self.xtb_keywords_input.text()))
        self.xtb_keywords_input.setToolTip("CREST ONLY: Define additional keywords to use in the xTB pre-optimization that are not included in -c, --uhf, -P and --input. For example: '--alpb ch2cl2 --gfn 1'.")
        advanced_layout.addWidget(self.xtb_keywords_input, 1, 8)

        self.cregen_checkbox = QCheckBox("CREGEN", self)
        self.cregen_checkbox.setChecked(True)
        self.cregen_checkbox.stateChanged.connect(lambda: self.control.update_crest_command("cregen", self.cregen_checkbox.isChecked()))
        self.cregen_checkbox.setToolTip("If True, perform a CREGEN analysis after CREST.")
        advanced_layout.addWidget(self.cregen_checkbox, 2, 7)

        self.cregen_keywords_input = QLineEdit(self)
        self.cregen_keywords_input.setPlaceholderText("CREGEN keywords...")
        self.cregen_keywords_input.textChanged.connect(lambda: self.control.update_crest_command("cregen_keywords", self.cregen_keywords_input.text()))
        self.cregen_keywords_input.setToolTip("Additional keywords for CREGEN (i.e. cregen_keywords='--ethr 0.02').")
        advanced_layout.addWidget(self.cregen_keywords_input, 2, 8)

        self.crest_runs_label = QLabel("CREST runs:", self)
        advanced_layout.addWidget(self.crest_runs_label, 3, 7)

        self.crest_runs_input = QSpinBox(self)
        self.crest_runs_input.setRange(1, 100)
        self.crest_runs_input.setValue(1)
        self.crest_runs_input.setSingleStep(1)
        self.crest_runs_input.valueChanged.connect(lambda: self.control.update_crest_command("crest_runs", self.crest_runs_input.value()))
        self.crest_runs_input.setToolTip("Specify as number of runs if multiple starting points from RDKit starting points is required.")
        advanced_layout.addWidget(self.crest_runs_input, 3, 8)


        self.csv_model.signals.updated.connect(self.update_ui)
        self.update_ui()

    def _molecule_label_mouse_press(self,event):
        if event.button() == Qt.MouseButton.LeftButton:
            pos = event.position()
            self.control.mousePressEvent(pos)    
        QLabel.mousePressEvent(self.molecule_label, event)

    def toggle_panel(self, height, width):
        expanded_height = 200
        if self.advanced_settings_button.isChecked():
            self.advanced_panel.setFixedHeight(expanded_height)
            self.resize(width, height + expanded_height)
            self.advanced_settings_button.setIcon(QIcon(Icons.eye_crossed))
        else:
            self.advanced_panel.setFixedHeight(0)
            self.resize(width, height- expanded_height)
            self.advanced_settings_button.setIcon(QIcon(Icons.eye))

    def handle_property_change(self, item):
        """Handle changes to editable properties and update the csv_dictionary."""
        logging.debug("at handle_property_change >>> handling property change")
        if item.row() == 1:  
            self.csv_model["code_name"][self.control.current_index - 1] = item.text()
        elif item.row() == 2:  
            self.csv_model["charge"][self.control.current_index - 1] = item.text()
            # Track user-defined charge
            if not hasattr(self, 'user_defined_charge'):
                self.user_defined_charge = {}
            self.user_defined_charge[self.control.current_index] = item.text()
        elif item.row() == 3:  
            self.csv_model["multiplicity"][self.control.current_index - 1] = item.text()
            # Track user-defined multiplicity
            if not hasattr(self, 'user_defined_multiplicity'):
                self.user_defined_multiplicity = {}
            self.user_defined_multiplicity[self.control.current_index] = item.text()

    def smiles_from_pubchem(self):
        """
        Fetches the canonical SMILES string for a compound from PubChem using its CID, CAS, or name, and updates the input field and model accordingly.
        This function is triggered when the user presses Enter in the search field.
        It performs the following steps:
        - Retrieves the user input (CID, CAS, or compound name) from the search field.
        - Attempts to fetch the SMILES string from PubChem
        - Updates the SMILES field in the model and input box, appending if necessary.
        - Sets the compound name in the model if not already present.
        """
        code_name = self.search_pubchem_input.text()
        if not code_name:
            self.failure("Please enter a CID, CAS, or compound name to search PubChem.")
            return
        try:
            smiles = pubchem2smiles(code_name)
            current_text = self.csv_model.__getitem__("SMILES")[self.control.current_index - 1]
            new_text = f"{current_text}.{smiles}" if current_text else smiles # Allow for multople molecules 
            if isinstance(new_text, str): # To keep pylance happy
                self.smiles_input.setPlainText(new_text)

        except IndexError:
            self.failure("No compound found for the given input. Please check the CID or name.")
        except Exception as e:
            self.failure(f"An error occurred for {code_name}: {str(e)}")
            
        if not self.csv_model.__getitem__("code_name")[self.control.current_index - 1]:
            self.csv_model.set_item_at_index("code_name", self.control.current_index - 1, code_name.replace(" ", "_"))

        self.update_properties()
        self.search_pubchem_input.clear()

    def resizeEvent(self, event):
        """Handle window resize events to refresh the molecule display."""
        super().resizeEvent(event)
        self.control.display_molecule(self.show_numbered_atoms_toggle.isChecked())

    def mousePressEvent(self, event: QMouseEvent):
        if event.button() == Qt.MouseButton.LeftButton:
            pos = event.position()
            self.control.mousePressEvent(pos)

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
            current_index = self.control.current_index - 1 # Adjust for 0-based index
            try:
                properties = self.csv_model.get_row_as_list_of_tuples(current_index)
            except IndexError:
                logging.debug("IndexError in update_properties: current_index out of range.")
                properties = []
            
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
                    if row in  [1, 2, 3]:
                        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsEditable)
                    else:  
                        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)

            combo = self.properties_table.cellWidget(8, 1)
            if combo:  # Check if the widget exists
                if self.metal_detected_in_smiles:
                    combo.setEnabled(True)
                    complex_type = self.csv_model.__getitem__("complex_type")[self.control.current_index - 1]
                    index = combo.findText(complex_type)
                    if index != -1:
                        combo.setCurrentIndex(index)
                else:
                    combo.setEnabled(False)

            self.properties_table.blockSignals(False)
        finally:
            self.is_programmatic_update = old_value
    
    def handle_combobox_change(self, row, text):
        """Handle changes in the complex_type combobox and update the model."""
        if row == 8:
            self.csv_model.set_item_at_index("complex_type", self.control.current_index - 1, text)

# UI ELEMENTS UPDATE FUNCTIONS (This will definitely stay but might need to rewrite quite heavily)

    @Slot()
    def update_ui(self):
        """Update all UI elements to reflect the current molecule's data."""
        try:
            smiles = self.csv_model.__getitem__("SMILES")[self.control.current_index - 1]
        except IndexError:
            smiles = "" 
        self.smiles_input.blockSignals(True)
        self.smiles_input.setPlainText(smiles)

        try:
            self.smiles_output.setText(smiles2enumerate(smiles))
        except ValueError:
            self.smiles_output.setText("")

        cursor = self.smiles_input.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        self.smiles_input.setTextCursor(cursor)
        self.smiles_input.blockSignals(False)

        self.control.display_molecule(self.show_numbered_atoms_toggle.isChecked())

        self.index_and_total_label_update()

        try:
            num_atoms = smiles2numatoms(smiles)
            num_electrons = smiles2numelectrons(smiles)
            self.atom_electron_label.setText(f" Electrons: {num_electrons}\n Atoms: {num_atoms}")
        except Exception as e:
            print(f"Error calculating atom/electron counts: {e}")
            self.atom_electron_label.setText(" Electrons: 0\n Atoms: 0")

        self.metal_detected_in_smiles = smiles2ismetalcomplex(smiles)
        if self.metal_detected_in_smiles:
            self.log_box_label.setText(f"Metal atom(s) detected in SMILES: {smiles2findmetal(smiles)} Please specify the complex type.")
        else:
            self.log_box_label.setText("")
        self.update_properties()
        
    def smiles_are_bad_bro(self, smiles: str):
        self.log_box_label.setText(f"Invalid SMILES string: {smiles}. ")

    def index_and_total_label_update(self):
        self.control.total_index = self.csv_model.get_total_index() 
        self.index_and_total_label.setText(f"{self.control.current_index}/{self.control.total_index}")

    def clear_focus_on_inputs(self):
        self.smiles_input.clearFocus()
        self.search_pubchem_input.clearFocus()

    def closeEvent(self, event):
        """Handle the close event to stop CSEARCH if running and prompt saving the CSV file."""
        pixmap = QPixmap(Icons.blue)
        icon = QIcon(pixmap)
        
        # First, check if CSEARCH is running
        if hasattr(self, 'worker') and self.parent.worker.isRunning():
            msgBox = QMessageBox(self)
            msgBox.setWindowTitle("CSEARCH Running")
            msgBox.setText("CSEARCH is still running. Stop it and continue closing?")
            msgBox.setWindowIcon(icon)
            msgBox.setIconPixmap(icon.pixmap(64, 64))
            msgBox.setStandardButtons(QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            reply = msgBox.exec()
            
            if reply == QMessageBox.StandardButton.Yes:
                # Stop the worker
                self.parent.worker.request_stop()
                self.parent.worker.wait(5000)  # Wait up to 5 seconds
                if self.parent.worker.isRunning():
                    self.parent.worker.terminate()
                    self.parent.worker.wait()
                self.shell_output.append("CSEARCH stopped due to application closure.")
            else:
                event.ignore()
                return
        
        # Then, check if CSV needs saving
        if self.file_name is None: 
            msgBox = QMessageBox(self)
            msgBox.setWindowTitle("Save CSV")
            msgBox.setText("Would you like to save the CSV file before exiting?")
            msgBox.setWindowIcon(icon)
            msgBox.setIconPixmap(icon.pixmap(64, 64))
            msgBox.setStandardButtons(QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel)
            reply = msgBox.exec()
            
            if reply == QMessageBox.StandardButton.Save:
                if not self.control.save_csv_file():
                    event.ignore()  # Prevent closing if saving fails
                else:
                    event.accept()
            elif reply == QMessageBox.StandardButton.Discard:
                event.accept()
            else:
                event.ignore()
        else:
            event.accept()

# AQME RUN SETUP FUNCTIONS

    def validate_smiles_for_csearch(self):
        """Check if SMILES are valid for CSEARCH."""
        for smiles in self.csv_model["SMILES"]:
            if "." in smiles and self.control.gen_command_model["program"] == "rdkit":
                self.failure("RDKit doesn't support multi-molecule SMILES (containing '.'). "
                            "Please use separate entries or try CREST.")
                return False
        return True

    def select_output_directory(self):
        """THIS SHOULD BE IN CONTROLLER, CONNECTED TO MODEL. 
        
        Select the output directory for AQME results."""
        self.output_directory = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if self.output_directory:
            self.output_dir_input.setText(f"{self.output_directory}")
        else:
            return
        
    def toggle_csearch(self):
        """Toggle between starting and stopping CSEARCH."""
        if self.parent.worker.isRunning():
            # Stop the running thread
            self.stop_csearch()
        else:
            # Start a new thread
            self.start_csearch()

    def start_csearch(self):
        """Start CSEARCH thread."""
        if not self.validate_smiles_for_csearch():
            return
        self.run_button.setText("Stop CSEARCH")
        self.run_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_MediaStop))
        
        # Reset stop flag
        self.parent.worker._stop_requested = False
        self.parent.worker.start()
        self.shell_output.append("Starting CSEARCH...")
        # Begin progress loop to simulate progress bar updates
        self.progress_bar.setRange(0, 0)  # Indeterminate mode

    def stop_csearch(self):
        """Request CSEARCH to stop."""
        self.run_button.setEnabled(False)
        self.run_button.setText("Stopping...")
        self.parent.worker.request_stop()
        
        # Give it more time (10 seconds) to kill all spawned processes
        if not self.parent.worker.wait(10000):  # Wait 10 seconds
            self.parent.worker.terminate()
            self.parent.worker.wait()
            self.shell_output.append("CSEARCH thread forcefully terminated.")
        
        # Ensure killer loop has stopped
        self.parent.worker._killer_active = False

    def handle_qprep_confirm(self, message: str, destination: str) -> None:
        """Ask user and open QPREP if they say yes."""
        if self.yes_no_dialog(message):
            self.parent.open_qprep_after_csearch(destination)

    def on_thread_complete(self):
        """Called when thread finishes."""
        self.run_button.setEnabled(True)
        self.run_button.setText("Run CSEARCH")
        self.run_button.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay))
        self.progress_bar.setRange(0, 100)

# random 
    @Slot(str)
    def success(self, message:str) -> None:
        pixmap = QPixmap(Icons.green)
        msg = QMessageBox(self)
        msg.setWindowTitle("Success")
        msg.setText(message)
        if not pixmap.isNull():
            icon = QIcon(pixmap)
            msg.setWindowIcon(icon)
            msg.setIconPixmap(icon.pixmap(64, 64))
        else:
            msg.setIcon(QMessageBox.Icon.Information)
        msg.exec()

    @Slot(str)
    def failure(self, message:str) -> None:
        pixmap = QPixmap(Icons.red)
        msg = QMessageBox(self)
        msg.setWindowTitle("Failure")
        msg.setText(message)
        if not pixmap.isNull():
            icon = QIcon(pixmap)
            msg.setWindowIcon(icon)
            msg.setIconPixmap(icon.pixmap(64, 64))
        else:
            msg.setIcon(QMessageBox.Icon.Critical)
        msg.exec()

    @Slot(str)
    def yes_no_dialog(self, message:str) -> bool:
        pixmap = QPixmap(Icons.blue)
        msgBox = QMessageBox(self)
        msgBox.setWindowTitle("Input Required")
        msgBox.setText(message)
        if not pixmap.isNull():
            icon = QIcon(pixmap)
            msgBox.setWindowIcon(icon)
            msgBox.setIconPixmap(icon.pixmap(64, 64))
        else:
            msgBox.setIcon(QMessageBox.Icon.Question)
        msgBox.setStandardButtons(QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        reply = msgBox.exec()
        return reply == QMessageBox.StandardButton.Yes
    
    def _drag_enter(self,event):
        mime = event.mimeData()
        if mime and mime.hasUrls():
            event.acceptProposedAction()
            try:
                self.molecule_label.setStyleSheet(stylesheets.MoleculeLabelHover)
                self.molecule_label.setText("Drop in a ChemDraw, CSV or SDF file to import or input SMILES")
            except Exception as e:
                print(f"file hovered over widget, exception thrown: {e}")
        else:
            event.ignore()

    def _drag_leave(self,event):
        try:
            self.molecule_label.setStyleSheet(stylesheets.MoleculeLabel)
        except Exception as e:
            print(f"file unhovered, exception thrown: {e}")

    def _drop_event(self,event):
        mime = event.mimeData()
        if mime and mime.hasUrls():
            event.acceptProposedAction()
            try:
                self.molecule_label.setStyleSheet(stylesheets.MoleculeLabel)
            except Exception as e:
                print(f"file dropped over widget, exception thrown: {e}")
            urls = mime.urls()
            if urls:
                file_path = urls[0].toLocalFile()
                self.control.import_file(file_path)
        else:
            event.ignore()


class MessageBox(QMessageBox):
    def __init__(self, parent=None, title="", message="", type=""):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setText(message)
        if type == "success":
            icon = QIcon(QPixmap(Icons.green))
        elif type == "failure":
            icon = QIcon(QPixmap(Icons.red))
        else:
            icon = QIcon(QPixmap(Icons.blue))
        
        self.setWindowIcon(icon)
        self.setIconPixmap(icon.pixmap(64, 64))
        self.setStandardButtons(QMessageBox.StandardButton.Ok)