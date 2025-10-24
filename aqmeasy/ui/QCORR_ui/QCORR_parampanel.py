import os
from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QPushButton,
    QLabel,
    QProgressBar,
    QCheckBox,
    QLineEdit,
    QSlider,
    QGroupBox,
    QDoubleSpinBox,
    QComboBox,
    QFileDialog,
)
from PySide6.QtCore import Qt
from aqmeasy.ui.QPREP_ui.QPREP_parameterpanel import ParameterPanel
from aqmeasy.models.QCORR_model.QCORR_parammodel import ParamModel
from aqmeasy.ui.stylesheets import stylesheets

class ParamPanel(QWidget):
    """Parameter panel for QCORR module"""

    def __init__(self):
        super().__init__()
        self.model = ParamModel()
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        param_group = QGroupBox("QCORR Parameters")
        param_group.setStyleSheet(stylesheets.QGroupBox)
        layout.addWidget(param_group)
        param_layout = QVBoxLayout()
        param_group.setLayout(param_layout)

        options_grid = QGridLayout()
        param_layout.addLayout(options_grid)

        
        ifreq_label = QLabel("Imaginary Frequency Cutoff (cm⁻¹): ")
        ifreq_label.setStyleSheet(stylesheets.QLabel)
        options_grid.addWidget(ifreq_label, 0, 0)
        self.ifreq_input = QDoubleSpinBox()
        self.ifreq_input.setRange(-100.0, 0.0)
        self.ifreq_input.setStyleSheet(stylesheets.QDoubleSpinBox)
        self.ifreq_input.setValue(self.model.ifreq_cutoff)
        self.ifreq_input.valueChanged.connect(lambda value: self.model.ifreqCutoffChanged.emit(float(value)))
        options_grid.addWidget(self.ifreq_input, 0, 1)

        ifreq_scale_label = QLabel(f"Amplitude Scale for Im. Frequencies: ")
        ifreq_scale_label.setStyleSheet(stylesheets.QLabel)
        options_grid.addWidget(ifreq_scale_label, 1, 0)
        slider_layout = QHBoxLayout()
        options_grid.addLayout(slider_layout, 1, 1)
        self.ifreq_scale_value = QLabel(f"{self.model.amplitude_ifreq:.2f}")
        self.ifreq_scale_value.setStyleSheet(stylesheets.QLabel)
        self.ifreq_scale_value.setFixedWidth(35)
        slider_layout.addWidget(self.ifreq_scale_value)

        self.ifreq_scale_slide = QSlider(Qt.Horizontal)
        self.ifreq_scale_slide.setStyleSheet(stylesheets.QSlider)
        self.ifreq_scale_slide.setRange(0,100)
        self.ifreq_scale_slide.setValue(int(self.model.amplitude_ifreq * 100))
        self.ifreq_scale_slide.setTickInterval(10)
        self.ifreq_scale_slide.valueChanged.connect(
            lambda value: self.model.amplitudeIfreqChanged.emit(float(value/100)))
        slider_layout.addWidget(self.ifreq_scale_slide)

        s2_threshold_label = QLabel("Spin Contamination Threshold (%): ")
        s2_threshold_label.setStyleSheet(stylesheets.QLabel)
        options_grid.addWidget(s2_threshold_label, 2, 0)
        self.s2_threshold_input = QDoubleSpinBox()
        self.s2_threshold_input.setStyleSheet(stylesheets.QDoubleSpinBox)
        self.s2_threshold_input.setRange(0.0, 100.0)
        self.s2_threshold_input.setValue(self.model.s2_threshold)
        options_grid.addWidget(self.s2_threshold_input, 2, 1)
        self.s2_threshold_input.valueChanged.connect(lambda value: self.model.s2ThresholdChanged.emit(float(value)))

        fullcheck_check = QCheckBox("Perform Full Check (level of theory)")
        fullcheck_check.setStyleSheet(stylesheets.QCheckBox)
        fullcheck_check.setChecked(self.model.fullcheck)
        options_grid.addWidget(fullcheck_check, 3, 0)
        self.fullcheck_check = fullcheck_check
        self.fullcheck_check.stateChanged.connect(lambda state: self.model.fullcheckChanged.emit(state == Qt.Checked))

        self.dup_check = QCheckBox("Disable Duplicate Filter")
        self.dup_check.setStyleSheet(stylesheets.QCheckBox)
        options_grid.addWidget(self.dup_check, 3, 1)
        self.dup_check.stateChanged.connect(lambda state: self.model.nodupCheckChanged.emit(state == Qt.Checked))
        self.dup_check.stateChanged.connect(lambda state: self.on_dup_check_state_changed(state))

        dup_group = QGroupBox("Duplicate Filter Parameters")
        dup_group.setStyleSheet(stylesheets.QGroupBox)
        param_layout.addWidget(dup_group)
        dup_layout = QGridLayout()
        dup_group.setLayout(dup_layout)

        dupThreshold_label = QLabel("Duplicate Energy Threshold (Hartree): ")
        dupThreshold_label.setStyleSheet(stylesheets.QLabel)
        dup_layout.addWidget(dupThreshold_label, 0, 0)
        self.dupThreshold_input = QDoubleSpinBox()
        self.dupThreshold_input.setStyleSheet(stylesheets.QDoubleSpinBox)
        self.dupThreshold_input.setRange(0.0, 1.0)
        self.dupThreshold_input.setSingleStep(0.0001)
        self.dupThreshold_input.setDecimals(4)
        self.dupThreshold_input.setValue(self.model.dup_threshold)
        dup_layout.addWidget(self.dupThreshold_input, 0, 1)
        self.dupThreshold_input.valueChanged.connect(lambda value: self.model.dupThresholdChanged.emit(float(value)))

        roThreshold_label = QLabel("Rotational Constant Threshold (cm⁻¹): ")
        roThreshold_label.setStyleSheet(stylesheets.QLabel)
        dup_layout.addWidget(roThreshold_label, 1, 0)
        self.roThreshold_input = QDoubleSpinBox()
        self.roThreshold_input.setStyleSheet(stylesheets.QDoubleSpinBox)
        self.roThreshold_input.setRange(0.0, 10.0)
        self.roThreshold_input.setSingleStep(0.1)
        self.roThreshold_input.setValue(self.model.ro_threshold)
        dup_layout.addWidget(self.roThreshold_input, 1, 1)
        self.roThreshold_input.valueChanged.connect(lambda value: self.model.roThresholdChanged.emit(float(value)))

        isom_group = QGroupBox("Isomerisation Filter Parameters")
        isom_group.setStyleSheet(stylesheets.QGroupBox)
        param_layout.addWidget(isom_group)
        isom_layout = QGridLayout()
        isom_group.setLayout(isom_layout)

        isom_type_label = QLabel("Isomerisation Input File Type: ")
        isom_type_label.setStyleSheet(stylesheets.QLabel)
        isom_layout.addWidget(isom_type_label, 0, 0)
        self.isom_type_input = QComboBox()
        self.isom_type_input.setStyleSheet(stylesheets.QComboBox)
        self.isom_type_input.addItems(["None", ".out", ".log", ".xyz"])
        isom_layout.addWidget(self.isom_type_input, 0, 1)
        self.isom_type_input.currentTextChanged.connect(lambda value: self.model.isomTypeChanged.emit(value))

        self.isom_inputs_label = QLabel("Isomerisation Input Files (optional): ")
        self.isom_inputs_label.setStyleSheet(stylesheets.QLabel)
        isom_layout.addWidget(self.isom_inputs_label, 1, 0)

        # input line for isomerisation input files
        isom_input_layout = QHBoxLayout()
        isom_layout.addLayout(isom_input_layout, 1, 1)
        self.isom_input_line = QLineEdit()
        self.isom_input_line.setStyleSheet(stylesheets.QLineEdit)
        isom_input_layout.addWidget(self.isom_input_line)
        self.isom_input_line.setReadOnly(True)
        # add a button to open file dialog
        self.isom_input_button = QPushButton("Browse")
        self.isom_input_button.setStyleSheet(stylesheets.QPushButton)
        isom_input_layout.addWidget(self.isom_input_button)
        # self.isom_input_button.clicked.connect(self.browse_isom_inputs)

        vdwfrac_label = QLabel("van der Waals Radii for Bond Detection: ")
        vdwfrac_label.setStyleSheet(stylesheets.QLabel)
        isom_layout.addWidget(vdwfrac_label, 2, 0)
        self.vdwfrac_input = QDoubleSpinBox()
        self.vdwfrac_input.setStyleSheet(stylesheets.QDoubleSpinBox)
        self.vdwfrac_input.setRange(0.0, 3.0)
        self.vdwfrac_input.setSingleStep(0.01)
        self.vdwfrac_input.setValue(self.model.vdwfrac)
        isom_layout.addWidget(self.vdwfrac_input, 2, 1)
        self.vdwfrac_input.valueChanged.connect(lambda value: self.model.vdwfracChanged.emit(float(value)))

        casevfrac_label = QLabel("Covalent Radii for Bond Detection: ")
        casevfrac_label.setStyleSheet(stylesheets.QLabel)
        isom_layout.addWidget(casevfrac_label, 3, 0)
        self.covfrac_input = QDoubleSpinBox()
        self.covfrac_input.setStyleSheet(stylesheets.QDoubleSpinBox)
        self.covfrac_input.setRange(0.0, 3.0)
        self.covfrac_input.setSingleStep(0.01)
        self.covfrac_input.setValue(self.model.covfrac)
        isom_layout.addWidget(self.covfrac_input, 3, 1)
        self.covfrac_input.valueChanged.connect(lambda value: self.model.covfracChanged.emit(float(value)))

        qm_input_group = QGroupBox("QM Input Parameters")
        qm_input_group.setStyleSheet(stylesheets.QGroupBox)
        param_layout.addWidget(qm_input_group)
        qm_input_layout = QHBoxLayout()
        qm_input_group.setLayout(qm_input_layout)

        qm_input = QLabel("QM Input Parameters (optional): ")
        qm_input.setStyleSheet(stylesheets.QLabel)
        qm_input_layout.addWidget(qm_input)

        qprep_button = QPushButton("Set QPREP Parameters")
        qprep_button.setStyleSheet(stylesheets.QPushButton)
        qm_input_layout.addWidget(qprep_button)
        qprep_button.clicked.connect(self.open_qprep_parameters)

        # Adjust view on model changes
        self.model.amplitudeIfreqChanged.connect(
            lambda value: self.ifreq_scale_value.setText(f"{value:.2f}")
        )
        self.model.s2ThresholdChanged.connect(
            lambda value: self.s2_threshold_input.setValue(value)
        )
        self.model.roThresholdChanged.connect(
            lambda value: self.roThreshold_input.setValue(value)
        )
        self.model.dupThresholdChanged.connect(
            lambda value: self.dupThreshold_input.setValue(value)
        )
        self.model.isomTypeChanged.connect(
            lambda value: self.isom_type_input.setCurrentText(value)
        )

    def open_qprep_parameters(self):
        self.qprep_dialog = qprep_parameters_dialog()
        self.qprep_dialog.show()
    
    def on_dup_check_state_changed(self, state):
        """Enables/disables duplicate threshold inputs based on checkbox state."""
        self.dupThreshold_input.setDisabled(state)
        self.roThreshold_input.setDisabled(state)


class qprep_parameters_dialog(QWidget):
    """Dialog to set QPREP parameters from QCORR"""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("QPREP Parameters")
        self.setStyleSheet(stylesheets.QWidget)
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)

        self.parameter_panel = ParameterPanel()
        layout.addWidget(self.parameter_panel)

        close_button = QPushButton("Close")
        close_button.setStyleSheet(stylesheets.QPushButton)
        close_button.clicked.connect(self.close)
        layout.addWidget(close_button)

        # Need to take data from qprep model in parampanel to use in qm_input in qcorr model
        # This can be done when the dialog is closed 
        # Might have to preset the QPREP model with whatever theory is being looked at but this can be done later