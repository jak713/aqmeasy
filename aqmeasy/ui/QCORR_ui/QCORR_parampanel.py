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
    QSlider
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
        
        options_grid = QGridLayout()
        layout.addLayout(options_grid)

        ifreq_label = QLabel("Imaginary Frequency Cutoff (cm⁻¹): ")
        ifreq_label.setStyleSheet(stylesheets.QLabel)
        options_grid.addWidget(ifreq_label, 0, 0)
        self.ifreq_input = QLineEdit(f"{self.model.ifreq_cutoff}")
        self.ifreq_input.setStyleSheet(stylesheets.QLineEdit)
        options_grid.addWidget(self.ifreq_input, 0, 1)

        ifreq_scale_label = QLabel(f"Amplitude Scale for Im. Frequencies ({self.model.amplitude_ifreq}) : ")
        ifreq_scale_label.setStyleSheet(stylesheets.QLabel)
        options_grid.addWidget(ifreq_scale_label, 1, 0)
        self.ifreq_scale_slide = QSlider(Qt.Horizontal)
        self.ifreq_scale_slide.setStyleSheet(stylesheets.QSlider)
        self.ifreq_scale_slide.setRange(0,100)
        self.ifreq_scale_slide.setValue(int(self.model.amplitude_ifreq * 100))
        self.ifreq_scale_slide.setTickInterval(10)
        self.ifreq_scale_slide.setTickPosition(QSlider.TicksBelow)        
        options_grid.addWidget(self.ifreq_scale_slide, 1, 1)

        s2_threshold_label = QLabel("Spin Contamination Threshold (%): ")
        s2_threshold_label.setStyleSheet(stylesheets.QLabel)
        options_grid.addWidget(s2_threshold_label, 2, 0)
        self.s2_threshold_input = QLineEdit(f"{self.model.s2_threshold}")
        self.s2_threshold_input.setStyleSheet(stylesheets.QLineEdit)
        options_grid.addWidget(self.s2_threshold_input, 2, 1)

        self.dup_check = QCheckBox("Enable Duplicate Filter")
        self.dup_check.setStyleSheet(stylesheets.QCheckBox)
        options_grid.addWidget(self.dup_check, 3, 0)
        self.dup_check.stateChanged.connect(lambda state: self.on_dup_check_state_changed(state))

        qprep_button = QPushButton("Set QPREP Parameters")
        qprep_button.setStyleSheet(stylesheets.QPushButton)
        layout.addWidget(qprep_button)
        qprep_button.clicked.connect(self.open_qprep_parameters)

    def open_qprep_parameters(self):
        self.qprep_dialog = qprep_parameters_dialog()
        self.qprep_dialog.show()
    
    def on_dup_check_state_changed(self, state):
        pass


class qprep_parameters_dialog(QWidget):
    """Dialog to set QPREP parameters from QCORR"""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("QPREP Parameters")
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)

        self.parameter_panel = ParameterPanel()
        layout.addWidget(self.parameter_panel)