import os
from PySide6.QtWidgets import (
    QApplication,
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
    QSizePolicy,
    QStyle,
    QFrame
)
from PySide6.QtGui import QIcon
from PySide6.QtCore import Qt
from aqmeasy.ui.QPREP_ui.QPREP_parameterpanel import ParameterPanel
from aqmeasy.models.QCORR_model.QCORR_parammodel import ParamModel
from aqmeasy.ui.stylesheets import stylesheets
from aqmeasy.ui.icons import Icons


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

        checkboxes_layout = QHBoxLayout()
        param_layout.addLayout(checkboxes_layout)

        fullcheck_check = QCheckBox("Perform Full Check")
        fullcheck_check.setChecked(self.model.fullcheck)
        checkboxes_layout.addWidget(fullcheck_check)
        self.fullcheck_check = fullcheck_check
        self.fullcheck_check.stateChanged.connect(lambda state: self.model.fullcheckChanged.emit(state == Qt.Checked))

        self.dup_check = QCheckBox("Disable Duplicate Filter")
        checkboxes_layout.addWidget(self.dup_check)
        self.dup_check.stateChanged.connect(lambda state: self.model.nodupCheckChanged.emit(state == Qt.Checked))
        self.dup_check.stateChanged.connect(lambda state: self.on_dup_check_state_changed(state))

        reset_button = QPushButton("Reset Settings")
        checkboxes_layout.addWidget(reset_button)
        # reset_button.clicked.connect(self.model.resetParameters.emit)

        options_grid = QGridLayout()
        param_layout.addLayout(options_grid)
        
        ifreq_label = QLabel("Imaginary Frequency Cutoff (cm⁻¹): ")
        options_grid.addWidget(ifreq_label, 0, 0)
        self.ifreq_input = QDoubleSpinBox()
        self.ifreq_input.setRange(-100.0, 0.0)
        self.ifreq_input.setValue(self.model.ifreq_cutoff)
        self.ifreq_input.valueChanged.connect(lambda value: self.model.ifreqCutoffChanged.emit(float(value)))
        options_grid.addWidget(self.ifreq_input, 0, 1)

        ifreq_scale_label = QLabel(f"Amplitude Scale for Im. Frequencies: ")
        options_grid.addWidget(ifreq_scale_label, 1, 0)
        slider_layout = QHBoxLayout()
        options_grid.addLayout(slider_layout, 1, 1)
        self.ifreq_scale_value = QLabel(f"{self.model.amplitude_ifreq:.2f}")
        self.ifreq_scale_value.setFixedWidth(35)
        slider_layout.addWidget(self.ifreq_scale_value)

        self.ifreq_scale_slide = QSlider(Qt.Horizontal)
        self.ifreq_scale_slide.setRange(0,100)
        self.ifreq_scale_slide.setValue(int(self.model.amplitude_ifreq * 100))
        self.ifreq_scale_slide.setTickInterval(10)
        self.ifreq_scale_slide.valueChanged.connect(
            lambda value: self.model.amplitudeIfreqChanged.emit(float(value/100)))
        slider_layout.addWidget(self.ifreq_scale_slide)

        s2_threshold_label = QLabel("Spin Contamination Threshold (%): ")
        options_grid.addWidget(s2_threshold_label, 2, 0)
        self.s2_threshold_input = QDoubleSpinBox()
        self.s2_threshold_input.setRange(0.0, 100.0)
        self.s2_threshold_input.setValue(self.model.s2_threshold)
        options_grid.addWidget(self.s2_threshold_input, 2, 1)
        self.s2_threshold_input.valueChanged.connect(lambda value: self.model.s2ThresholdChanged.emit(float(value)))

        dup_group = QGroupBox("Duplicate Filter Parameters")
        param_layout.addWidget(dup_group)
        dup_layout = QGridLayout()
        dup_group.setLayout(dup_layout)

        dupThreshold_label = QLabel("Duplicate Energy Threshold (Hartree): ")
        dup_layout.addWidget(dupThreshold_label, 0, 0)
        self.dupThreshold_input = QDoubleSpinBox()
        self.dupThreshold_input.setRange(0.0, 1.0)
        self.dupThreshold_input.setSingleStep(0.0001)
        self.dupThreshold_input.setDecimals(4)
        self.dupThreshold_input.setValue(self.model.dup_threshold)
        dup_layout.addWidget(self.dupThreshold_input, 0, 1)
        self.dupThreshold_input.valueChanged.connect(lambda value: self.model.dupThresholdChanged.emit(float(value)))

        roThreshold_label = QLabel("Rotational Constant Threshold (cm⁻¹): ")
        dup_layout.addWidget(roThreshold_label, 1, 0)
        self.roThreshold_input = QDoubleSpinBox()
        self.roThreshold_input.setRange(0.0, 10.0)
        self.roThreshold_input.setSingleStep(0.1)
        self.roThreshold_input.setValue(self.model.ro_threshold)
        dup_layout.addWidget(self.roThreshold_input, 1, 1)
        self.roThreshold_input.valueChanged.connect(lambda value: self.model.roThresholdChanged.emit(float(value)))

        isom_group = QGroupBox("Isomerisation Filter Parameters")
        param_layout.addWidget(isom_group)
        isom_layout = QGridLayout()
        isom_group.setLayout(isom_layout)

        isom_type_label = QLabel("Isomerisation Input File Type: ")
        isom_layout.addWidget(isom_type_label, 0, 0)
        self.isom_type_input = QComboBox()
        self.isom_type_input.addItems(["None", ".out", ".log", ".xyz", ".sdf", ".mol"])
        isom_layout.addWidget(self.isom_type_input, 0, 1)
        self.isom_type_input.currentTextChanged.connect(lambda value: self.model.isomTypeChanged.emit(value))

        self.isom_inputs_label = QLabel("Isomerisation Input Files (optional): ")
        isom_layout.addWidget(self.isom_inputs_label, 1, 0)

        # input line for isomerisation input files
        isom_input_layout = QHBoxLayout()
        isom_layout.addLayout(isom_input_layout, 1, 1)
        self.isom_input_line = QLineEdit()
        isom_input_layout.addWidget(self.isom_input_line)
        self.isom_input_line.setReadOnly(True)
        # add a button to open file dialog
        self.isom_input_button = QPushButton()
        self.isom_input_button.setIcon(QIcon(Icons.file_open))
        isom_input_layout.addWidget(self.isom_input_button)
        # self.isom_input_button.clicked.connect(self.browse_isom_inputs)

        vdwfrac_label = QLabel("van der Waals Radii for Bond Detection: ")
        isom_layout.addWidget(vdwfrac_label, 2, 0)
        self.vdwfrac_input = QDoubleSpinBox()
        self.vdwfrac_input.setRange(0.0, 3.0)
        self.vdwfrac_input.setSingleStep(0.01)
        self.vdwfrac_input.setValue(self.model.vdwfrac)
        isom_layout.addWidget(self.vdwfrac_input, 2, 1)
        self.vdwfrac_input.valueChanged.connect(lambda value: self.model.vdwfracChanged.emit(float(value)))

        casevfrac_label = QLabel("Covalent Radii for Bond Detection: ")
        isom_layout.addWidget(casevfrac_label, 3, 0)
        self.covfrac_input = QDoubleSpinBox()
        self.covfrac_input.setRange(0.0, 3.0)
        self.covfrac_input.setSingleStep(0.01)
        self.covfrac_input.setValue(self.model.covfrac)
        isom_layout.addWidget(self.covfrac_input, 3, 1)
        self.covfrac_input.valueChanged.connect(lambda value: self.model.covfracChanged.emit(float(value)))

        qm_input_group = QGroupBox("QM Input Parameters")
        param_layout.addWidget(qm_input_group)
        qm_input_layout = QHBoxLayout()
        qm_input_group.setLayout(qm_input_layout)

        qm_input = QLabel("QM Input Parameters (optional): ")
        qm_input_layout.addWidget(qm_input)

        qprep_button = QPushButton("Set QPREP Parameters")
        qm_input_layout.addWidget(qprep_button)
        qprep_button.clicked.connect(self.open_qprep_parameters)

        advanced_settings_button = QPushButton("QPREP Settings")
        layout.addWidget(advanced_settings_button)
        advanced_settings_button.setIcon(QIcon(Icons.eye))
        advanced_settings_button.setCheckable(True)
        advanced_settings_button.clicked.connect(lambda: self.toggle_panel(self.height(), self.width()))
        self.advanced_settings_button = advanced_settings_button

        self.advanced_panel = QFrame()
        self.advanced_panel.setLayout(QVBoxLayout())
        self.advanced_panel.layout().addWidget(param_group)
        self.advanced_panel.setFixedHeight(0)
        layout.addWidget(self.advanced_panel)

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

    def toggle_panel(self, height, width):
        expanded_height = 500
        if self.advanced_settings_button.isChecked():
            self.advanced_panel.setFixedHeight(expanded_height)
            self.resize(width, height + expanded_height)
            self.advanced_settings_button.setIcon(QIcon(Icons.eye_crossed))
        else:
            self.advanced_panel.setFixedHeight(0)
            self.resize(width, height- expanded_height)
            self.advanced_settings_button.setIcon(QIcon(Icons.eye))


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