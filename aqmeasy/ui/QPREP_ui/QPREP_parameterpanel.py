from PySide6.QtWidgets import QWidget, QComboBox, QLabel, QGroupBox, QVBoxLayout, QHBoxLayout, QSpinBox, QSlider
from PySide6.QtCore import Qt
from aqmeasy.models.QPREP_model import InputModel
from aqmeasy.ui.stylesheets import stylesheets

class ParameterPanel(QWidget): 
    """Parameter selection panel, Deals with all calculation selections for GAUSSIAN and ORCA """
    def __init__(self, model=None):
        super().__init__()
        self.model = model or InputModel()
        self.setup_software_data()
        self.setup_ui()
    
    def setup_software_data(self):
        """Define software-specific functionals and basis sets"""
        self.ORCA_FUNCTIONALS = ["HF", "MP2", "CCSD", "CCSD(T)", "BLYP", "PBE",
            "revPBE", "B3LYP", "MO6L", "MO62X", "B97-3C"
        ]
        
        self.ORCA_BASIS_SETS = [
            '6-31G(d)', 'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'aug-cc-pVDZ',
            'aug-cc-pVTZ', 'aug-cc-pVQZ', 'def2-SVP', 'def2-TZVP', 'def2-QZVP',
            'def2-TVZPP', 'def2-QZVPP', 'def2-TZVPPD', 'def2-QZVPPD', 'ma-def2-SVP',
            'ma-def2-TZVP', 'ma-def2-QZVP',
        ]
        
        self.GAUSSIAN_FUNCTIONALS =['APFD', 'B3LYP', 'BPV86', 'B3PW91', 'CAM-B3LYP', 'HCTH', 'HSEH1PBE', 'LSDA', 'MPW1PW91', 'PBEPBE', 'TPSSTPSS', 'WB97XD']
        
        self.GAUSSIAN_BASIS_SETS = [
            'STO-3G', '3-21G', '6-31G(d)', '6-31G(d,p)', 'LANL2DZ',
            'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'cc-p5VZ', 'cc-p6VZ',
            'aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'aug-cc-p5VZ', 'aug-cc-p6VZ',
            'Def2SV', 'Def2TZV', 'Def2QZV', 'Def2SVP', 'Def2TZVP',
            'Def2QZVP', 'Def2SVPP', 'Def2TZVPP', 
        ]
    
    
    def setup_ui(self):
        layout = QVBoxLayout()
        
        # Software selection
        software_group = QGroupBox("Software")
        software_group.setStyleSheet(stylesheets.QGroupBox)
        software_layout = QHBoxLayout()
        
        self.software_combo = QComboBox()
        self.software_combo.setStyleSheet(stylesheets.QComboBox)
        self.software_combo.addItems(["Orca", "Gaussian"]) 
        self.software_combo.currentTextChanged.connect(self.on_software_changed)
        self.software_combo.currentTextChanged.connect(self.model.setSoftware)
        software_layout.addWidget(QLabel("Software:",))
        software_layout.addWidget(self.software_combo)
        software_group.setLayout(software_layout)
        
        # Functional selection
        functional_group = QGroupBox("Method")
        functional_group.setStyleSheet(stylesheets.QGroupBox)
        functional_layout = QVBoxLayout()
        
        func_row = QHBoxLayout()
        func_row.addWidget(QLabel("Functional"))
        self.functional_combo = QComboBox()
        self.functional_combo.setStyleSheet(stylesheets.QComboBox)
        self.functional_combo.currentTextChanged.connect(self.model.setFunctional)
        func_row.addWidget(self.functional_combo)
        functional_layout.addLayout(func_row)
        
        basis_row = QHBoxLayout()
        basis_row.addWidget(QLabel("Basis Set"))
        self.basis_combo = QComboBox()
        self.basis_combo.setStyleSheet(stylesheets.QComboBox)
        self.basis_combo.currentTextChanged.connect(self.model.setBasisSet)
        basis_row.addWidget(self.basis_combo)
        functional_layout.addLayout(basis_row)
        
        functional_group.setLayout(functional_layout)
        
        # Processor count selection
        processor_group = QGroupBox("Computational Resources")
        processor_group.setStyleSheet(stylesheets.QGroupBox)
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
        self.nprocs_slider.valueChanged.connect(self.model.setNprocs)
        processor_layout.addWidget(self.nprocs_slider)
        
        processor_group.setLayout(processor_layout)

        # Memory in computational resources

        mem_row = QHBoxLayout()
        mem_row.addWidget(QLabel('Memory:'))
        self.mem_label = QLabel("1")
        mem_row.addWidget(self.mem_label)
        processor_layout.addLayout(mem_row)

        self.mem_slider = QSlider(Qt.Horizontal)
        self.mem_slider.setStyleSheet(stylesheets.QSlider)
        self.mem_slider.setRange(1, 16)
        self.mem_slider.setValue(1)
        self.mem_slider.valueChanged.connect(self.on_mem_changed)
        self.mem_slider.valueChanged.connect(self.model.setMem)
        processor_layout.addWidget(self.mem_slider)

        processor_group.setLayout(processor_layout)
        
        # Solvation model
        solvation_group = QGroupBox("Solvation")
        solvation_group.setStyleSheet(stylesheets.QGroupBox)
        solvation_layout = QVBoxLayout()
        
        self.ORCA_SOLVENT_MODELS = ["None", "CPCM", "SMD", "COSMORS", "DRACO"]
        self.ORCA_SOLVENTS = [
            "None", "Water", "Methanol", "Ethanol", "Acetone", "DMSO", 'Toluene', 'Aniline', 'Benzene', 'Chloroform', 'Carbon Disulfide', 'DCM', 'diethyl ether', 'DMF', 'Ethyl Acetate', 'Nitromethane', 'THF',
        ]
        self.GAUSSIAN_SOLVENT_MODELS = ["None", "PCM", "SMD", "IEFPCM", "CPCM"]
        self.GAUSSIAN_SOLVENTS = [
            "None", "Water", "Methanol", "Ethanol", "Acetone", "DMSO", "Benzene", "Chloroform", "Toluene", "Acetonitrile", "Dichloromethane", "DiethylEther", "Hexane", "Heptane", "Octanol", "THF", "DMF", "Ethyl Acetate", "Nitromethane"
        ]

        solv_row = QHBoxLayout()
        solv_row.addWidget(QLabel("Solvent Model"))
        self.solvation_combo = QComboBox()
        self.solvation_combo.setStyleSheet(stylesheets.QComboBox)
        self.solvation_combo.addItems(self.ORCA_SOLVENT_MODELS)
        solv_row.addWidget(self.solvation_combo)
        solvation_layout.addLayout(solv_row)
        
        solvent_row = QHBoxLayout()
        solvent_row.addWidget(QLabel("Solvent"))
        self.solvent_combo = QComboBox()
        self.solvent_combo.setStyleSheet(stylesheets.QComboBox)
        self.solvent_combo.addItems(self.ORCA_SOLVENTS)
        solvent_row.addWidget(self.solvent_combo)
        solvation_layout.addLayout(solvent_row)
        
        solvation_group.setLayout(solvation_layout)

        # Initialize with default software (Orca)
        self.on_software_changed("Orca")
        
        # Calculation type
        calc_group = QGroupBox("Calculation Type")
        calc_group.setStyleSheet(stylesheets.QGroupBox)
        calc_layout = QVBoxLayout()
        
        calc_row = QHBoxLayout()
        calc_row.addWidget(QLabel("Job Type:"))
        self.calc_combo = QComboBox()
        self.calc_combo.setStyleSheet(stylesheets.QComboBox)
        self.calc_combo.addItems(["Single Point", "Geometry Optimization", "Frequency", "TD-DFT", "Opt+Freq"])
        calc_row.addWidget(self.calc_combo)
        calc_layout.addLayout(calc_row)
        
        charge_spin_row = QHBoxLayout()
        charge_spin_row.addWidget(QLabel("Charge:"))
        self.charge_spin = QSpinBox()
        self.charge_spin.setStyleSheet(stylesheets.QSpinBox)
        self.charge_spin.setRange(-10, 10)
        self.charge_spin.setValue(0)
        charge_spin_row.addWidget(self.charge_spin)
        
        charge_spin_row.addWidget(QLabel("Multiplicity:"))
        self.multiplicity_spin = QSpinBox()
        self.multiplicity_spin.setStyleSheet(stylesheets.QSpinBox)
        self.multiplicity_spin.setRange(1, 10)
        self.multiplicity_spin.setValue(1)
        charge_spin_row.addWidget(self.multiplicity_spin)
        calc_layout.addLayout(charge_spin_row)
        
        calc_group.setLayout(calc_layout)
        
        layout.addWidget(software_group)
        layout.addWidget(functional_group)
        layout.addWidget(processor_group)
        layout.addWidget(solvation_group)
        layout.addWidget(calc_group)
        layout.addStretch()
        
        self.setLayout(layout)
    
    def on_software_changed(self, software_name):
        self.functional_combo.clear()
        self.basis_combo.clear()
        
        if software_name == "Orca":
            self.functional_combo.addItems(self.ORCA_FUNCTIONALS)
            self.basis_combo.addItems(self.ORCA_BASIS_SETS)
            self.solvation_combo.clear()
            self.solvation_combo.addItems(self.ORCA_SOLVENT_MODELS)
            self.solvent_combo.clear()
            self.solvent_combo.addItems(self.ORCA_SOLVENTS)
        elif software_name == "Gaussian":
            self.functional_combo.addItems(self.GAUSSIAN_FUNCTIONALS)
            self.basis_combo.addItems(self.GAUSSIAN_BASIS_SETS)
            self.solvation_combo.clear()
            self.solvation_combo.addItems(self.GAUSSIAN_SOLVENT_MODELS)
            self.solvent_combo.clear()
            self.solvent_combo.addItems(self.GAUSSIAN_SOLVENTS)
            
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