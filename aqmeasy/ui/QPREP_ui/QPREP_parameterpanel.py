from PySide6.QtWidgets import QWidget, QComboBox, QLabel, QGroupBox, QVBoxLayout, QHBoxLayout, QSpinBox, QSlider, QCheckBox, QTextBrowser
from PySide6.QtCore import Qt
from aqmeasy.models.QPREP_model import InputModel
from aqmeasy.ui.stylesheets import stylesheets
from aqmeasy.software_data import Orca, Gaussian, TranslationDict

class ParameterPanel(QWidget): 
    """Parameter selection panel, Deals with all calculation selections for GAUSSIAN and ORCA """
    def __init__(self, model=None):
        super().__init__()
        self.model = model or InputModel()
        self.setup_software_data()
        self.setup_ui()
    
    def setup_software_data(self):
        """Define software-specific functionals and basis sets"""

        self.ORCA_METHODS = Orca.ORCA_METHODS

        self.ORCA_FUNCTIONALS = Orca.ORCA_FUNCTIONALS
        
        self.ORCA_BASIS_SETS = Orca.ORCA_BASIS_SETS

        self.ORCA_DISPERSION_CORRECTIONS = Orca.ORCA_DISPERSION_CORRECTIONS
        
        self.ORCA_SOLVENT_MODELS = Orca.ORCA_SOLVENT_MODELS
        
        self.ORCA_SOLVENTS = Orca.ORCA_SOLVENTS

        self.GAUSSIAN_METHODS = Gaussian.GAUSSIAN_METHODS
        
        self.GAUSSIAN_FUNCTIONALS = Gaussian.GAUSSIAN_FUNCTIONALS
        
        self.GAUSSIAN_BASIS_SETS = Gaussian.GAUSSIAN_BASIS_SETS

        self.GAUSSIAN_SOLVENT_MODELS = Gaussian.GAUSSIAN_SOLVENT_MODELS
        
        self.GAUSSIAN_SOLVENTS = Gaussian.GAUSSIAN_SOLVENTS
    
    def setup_ui(self):
        layout = QVBoxLayout()

        info_group = QGroupBox("Information")
        info_layout = QHBoxLayout()
        self.info_line = QTextBrowser()
        self.info_line.setMaximumHeight(100)
        info_layout.addWidget(self.info_line)
        info_group.setLayout(info_layout)
        
        # Software selection
        software_group = QGroupBox("Software")
        software_layout = QHBoxLayout()
        
        self.software_combo = QComboBox()
        self.software_combo.addItems(["Orca", "Gaussian"]) 
        self.software_combo.currentTextChanged.connect(self.on_software_changed)
        self.software_combo.currentTextChanged.connect(self.model.set_software)
        software_layout.addWidget(QLabel("Software:",))
        software_layout.addWidget(self.software_combo)
        software_group.setLayout(software_layout)
        
        # Functional selection
        functional_group = QGroupBox("Level of Theory")
        functional_layout = QVBoxLayout()
        
        method_row = QHBoxLayout()
        method_row.addWidget(QLabel("Method"))
        self.method_combo = QComboBox()
        self.method_combo.currentTextChanged.connect(lambda text: self.setMethod(text))
        method_row.addWidget(self.method_combo)
        functional_layout.addLayout(method_row)

        func_row = QHBoxLayout()
        func_row.addWidget(QLabel("Functional"))
        self.functional_combo = QComboBox()
        self.functional_combo.setDisabled(True)
        self.functional_combo.currentTextChanged.connect(self.model.set_functional)
        func_row.addWidget(self.functional_combo)
        functional_layout.addLayout(func_row)

        disp_row = QHBoxLayout()
        self.dispersion_scheme = QLabel("Dispersion Scheme")
        disp_row.addWidget(self.dispersion_scheme)
        self.dispersion_combo = QComboBox()
        self.dispersion_combo.setDisabled(True)
        self.dispersion_combo.currentTextChanged.connect(self.model.set_dispersion)
        disp_row.addWidget(self.dispersion_combo)
        functional_layout.addLayout(disp_row)

        basis_row = QHBoxLayout()
        basis_row.addWidget(QLabel("Basis Set"))
        self.basis_combo = QComboBox()
        self.basis_combo.currentTextChanged.connect(self.model.set_basis_set)
        basis_row.addWidget(self.basis_combo)
        functional_layout.addLayout(basis_row)
        
        functional_group.setLayout(functional_layout)
        
        # Processor count selection
        processor_group = QGroupBox("Computational Resources")

        processor_layout = QVBoxLayout()
        
        proc_row = QHBoxLayout()
        proc_row.addWidget(QLabel("Number of Processors:"))
        self.nprocs_label = QLabel("8")
        proc_row.addWidget(self.nprocs_label)
        processor_layout.addLayout(proc_row)
        
        self.nprocs_slider = QSlider(Qt.Horizontal)
        self.nprocs_slider.setStyleSheet(stylesheets.QSlider)
        self.nprocs_slider.setRange(1, 60)
        self.nprocs_slider.setValue(8)
        self.nprocs_slider.valueChanged.connect(self.on_nprocs_changed)
        self.nprocs_slider.valueChanged.connect(self.model.set_nprocs)
        processor_layout.addWidget(self.nprocs_slider)
        
        processor_group.setLayout(processor_layout)

        # Memory in computational resources

        mem_row = QHBoxLayout()
        self.mem_title = QLabel('Memory:')
        mem_row.addWidget(self.mem_title)
        self.mem_label = QLabel("1")
        mem_row.addWidget(self.mem_label)
        processor_layout.addLayout(mem_row)

        self.mem_slider = QSlider(Qt.Horizontal)
        self.mem_slider.setRange(1, 16)
        self.mem_slider.setValue(1)
        self.mem_slider.valueChanged.connect(self.on_mem_changed)
        self.mem_slider.valueChanged.connect(self.model.set_mem)
        processor_layout.addWidget(self.mem_slider)

        processor_group.setLayout(processor_layout)
        
        # Solvation model
        solvation_group = QGroupBox("Solvation")
        solvation_group.setStyleSheet(stylesheets.QGroupBox)
        solvation_layout = QVBoxLayout()
        

        solv_row = QHBoxLayout()
        solv_row.addWidget(QLabel("Solvent Model"))
        self.solvation_combo = QComboBox()
        self.solvation_combo.addItems(self.ORCA_SOLVENT_MODELS)
        self.solvation_combo.currentTextChanged.connect(self.model.set_solvent_model)
        solv_row.addWidget(self.solvation_combo)
        solvation_layout.addLayout(solv_row)
        
        solvent_row = QHBoxLayout()
        solvent_row.addWidget(QLabel("Solvent"))
        self.solvent_combo = QComboBox()
        self.solvent_combo.addItems(self.ORCA_SOLVENTS)
        self.solvent_combo.currentTextChanged.connect(self.model.set_solvent)
        solvent_row.addWidget(self.solvent_combo)
        solvation_layout.addLayout(solvent_row)
        
        solvation_group.setLayout(solvation_layout)

        # Calculation type
        calc_group = QGroupBox("Calculation Type")
        calc_layout = QVBoxLayout()
        
        calc_row = QHBoxLayout()
        calc_row.addWidget(QLabel("Job Type:"))
        self.calc_combo = QComboBox()
        self.calc_combo.addItems(Orca.ORCA_RUNTYPES)
        self.calc_combo.currentTextChanged.connect(self.set_info)
        calc_row.addWidget(self.calc_combo)
        calc_layout.addLayout(calc_row)

        # Initialize with default software (Orca)
        self.on_software_changed("Orca")
        
        self.disable_auto_charge_mult = QCheckBox("Disable Auto Charge/Multiplicity Detection")
        self.disable_auto_charge_mult.stateChanged.connect(lambda state: self.setDisableAutoChargeMult(state))
        calc_layout.addWidget(self.disable_auto_charge_mult)

        charge_spin_row = QHBoxLayout()
        charge_spin_row.addWidget(QLabel("Charge:"))
        self.charge_spin = QSpinBox()
        self.charge_spin.setRange(-10, 10)
        self.charge_spin.setValue(0)
        self.charge_spin.setDisabled(True)
        charge_spin_row.addWidget(self.charge_spin)
        
        charge_spin_row.addWidget(QLabel("Multiplicity:"))
        self.multiplicity_spin = QSpinBox()
        self.multiplicity_spin.setRange(1, 10)
        self.multiplicity_spin.setValue(1)
        self.multiplicity_spin.setDisabled(True)
        charge_spin_row.addWidget(self.multiplicity_spin)
        calc_layout.addLayout(charge_spin_row)
        
        calc_group.setLayout(calc_layout)
        layout.addWidget(info_group)
        layout.addWidget(software_group)
        layout.addWidget(functional_group)
        layout.addWidget(processor_group)
        layout.addWidget(solvation_group)
        layout.addWidget(calc_group)
        layout.addStretch()
        
        self.setLayout(layout)

        # connect model signals to update UI
        self.model.software_Changed.connect(self.on_software_changed)
        self.model.functional_Changed.connect(self.set_functional)
        self.model.basis_set_Changed.connect(self.basis_combo.setCurrentText)
        self.model.basis_set_Changed.connect(self.set_info)
        self.model.nprocs_Changed.connect(self.nprocs_slider.setValue)
        self.model.mem_Changed.connect(self.mem_slider.setValue) 
        self.model.solvent_model_Changed.connect(self.solvation_combo.setCurrentText)
        self.model.solvent_Changed.connect(self.solvent_combo.setCurrentText)
        self.model.dispersion_Changed.connect(self.dispersion_combo.setCurrentText)
    
    def on_software_changed(self, software_name):
        self.method_combo.clear()
        self.functional_combo.clear()
        self.basis_combo.clear()
        self.solvation_combo.clear()
        self.solvent_combo.clear()
        self.calc_combo.clear()
        
        if software_name == "Orca":
            self.method_combo.addItems(self.ORCA_METHODS)
            self.functional_combo.addItems(self.ORCA_FUNCTIONALS)
            self.basis_combo.addItems(self.ORCA_BASIS_SETS)
            self.dispersion_combo.addItems(self.ORCA_DISPERSION_CORRECTIONS)
            self.solvation_combo.addItems(self.ORCA_SOLVENT_MODELS)
            self.solvent_combo.addItems(self.ORCA_SOLVENTS)
            self.software_combo.setCurrentText("Orca")
            self.mem_title.setText('Memory (per core in GB):')
            self.calc_combo.addItems(Orca.ORCA_RUNTYPES)
        elif software_name == "Gaussian":
            self.method_combo.addItems(self.GAUSSIAN_METHODS)
            self.functional_combo.addItems(self.GAUSSIAN_FUNCTIONALS)
            self.basis_combo.addItems(self.GAUSSIAN_BASIS_SETS)
            self.solvation_combo.addItems(self.GAUSSIAN_SOLVENT_MODELS)
            self.solvent_combo.addItems(self.GAUSSIAN_SOLVENTS)
            self.software_combo.setCurrentText("Gaussian")
            self.mem_title.setText('Memory (total in GB):')
            self.calc_combo.addItems(Gaussian.GAUSSIAN_RUNTYPES)
            
    def on_nprocs_changed(self, value):
        self.nprocs_label.setText(str(value))
        
    def on_mem_changed(self, value):
        self.mem_label.setText(str(value))

    def get_current_software(self):
        return self.software_combo.currentText()

    def get_current_functional(self):
        return self.functional_combo.currentText()

    def get_current_basis_set(self):
        return self.basis_combo.currentText()
    
    def get_current_nprocs(self):
        return self.nprocs_slider.value()
    
    def get_current_mem(self):
        return self.mem_slider.value()

    def get_current_job_type(self):
        return self.calc_combo.currentText()

    def get_current_solvent_model(self):
        return self.solvation_combo.currentText()

    def get_current_solvent(self):
        return self.solvent_combo.currentText()

    def get_current_charge(self):
        return self.charge_spin.value()

    def get_current_multiplicity(self):
        return self.multiplicity_spin.value()

    def setDisableAutoChargeMult(self, state):
        self.charge_spin.setDisabled(not state)
        self.multiplicity_spin.setDisabled(not state)

    def setMethod(self, method):
        if not method == "DFT":
            self.functional_combo.setDisabled(True)
            self.dispersion_combo.setDisabled(True)
            self.model.set_functional(method)
        else:
            self.functional_combo.setDisabled(False)
            self.dispersion_combo.setDisabled(False)
            self.model.set_functional(self.functional_combo.currentText())

    def set_functional(self, functional):
        if functional in self.ORCA_FUNCTIONALS or functional in self.GAUSSIAN_FUNCTIONALS:
            self.method_combo.setCurrentText("DFT")
            self.functional_combo.setCurrentText(functional)
        else:
            self.method_combo.setCurrentText(functional)
        self.set_info() 

    def set_info(self):
        calc_type = TranslationDict.get(self.calc_combo.currentText())
        method_expanded = TranslationDict.get(method:=self.model.functional(), "Density Functional Theory")
        basis = self.model.basis_set()
        info = f"Preparing:\n-> {calc_type} calculation\n\nLevel of theory:\n-> {method}/{basis}\n\nNote: {method} refers to {method_expanded}."
        self.info_line.setText(info) 