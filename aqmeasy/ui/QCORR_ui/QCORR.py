"""
JUST FOR ME TO LOOK OVER PHILOSOPHY:
QCORR is a cclib-based module that detects issues and errors
in QM output files, structures all output data, and creates ready-to-submit input files to correct those issues. User-
specified criteria (i.e., spin contamination, isomerization, etc.) can be defined to filter output data. 

Typically, a tedious manual search and correction for error terminations, convergence issues, and extra imaginary
frequencies is necessary after running QM calculations. Based on our experience with structure optimizations and
frequency calculations for large databases (i.e., many thousands) of organic compounds, such occurrences are relatively
common. QCORR structures output data and automatically detects issues or errors, creating new input files that try to
correct those issues, a cycle that can be repeated several times 
"""


import os
import sys
import subprocess
import tempfile

from PySide6.QtWidgets import  QLabel, QVBoxLayout, QWidget, QPushButton, QHBoxLayout, QTextBrowser,QLineEdit, QTextEdit, QCheckBox, QMessageBox, QSizePolicy, QFileDialog, QTableWidget, QTableWidgetItem, QHeaderView, QApplication, QComboBox, QSpinBox, QStyle, QTableWidgetItem, QFrame, QGridLayout, QDoubleSpinBox
from PySide6.QtCore import Qt, QProcess
from PySide6.QtGui import QPixmap, QKeySequence, QShortcut, QMouseEvent, QIcon, QDoubleValidator, QTextCursor, QIntValidator
from aqmeasy.ui.QCORR_ui.QCORR_filepanel import FilePanel
from aqmeasy.ui.QCORR_ui.QCORR_viewpanel import ViewPanel
from aqmeasy.ui.QCORR_ui.QCORR_parampanel import ParamPanel
from aqmeasy.ui.QCORR_ui.QCORR_jsonpanel import JSONpanel

from aqmeasy.controllers.QCORR_controller import FileController
from aqmeasy.models.QCORR_model.QCORR_filemodel import FileModel

from aqmeasy.ui.stylesheets import stylesheets



class QCORR(QWidget):
    def __init__(self):
        super().__init__()
        self.file_model = FileModel()
        self.view_panel = ViewPanel()
        self.file_panel = FilePanel(self.file_model, self.view_panel)


        self.setStyleSheet(stylesheets.QWidget)
        # control.set_parent(self)
        self.setWindowTitle("QCORR")

        # main layout = horizontal box
        main_layout = QHBoxLayout()
        self.setLayout(main_layout)
        panel1 = QVBoxLayout()
        panel2 = QVBoxLayout()
        panel3 = QVBoxLayout()

        panel1.addWidget(self.file_panel)
        self.param_panel = ParamPanel()
        panel1.addWidget(self.param_panel)

        panel2.addWidget(self.view_panel)

        self.json_panel = JSONpanel()
        panel3.addWidget(self.json_panel)

        main_layout.addLayout(panel1)
        main_layout.addLayout(panel2)
        main_layout.addLayout(panel3)

