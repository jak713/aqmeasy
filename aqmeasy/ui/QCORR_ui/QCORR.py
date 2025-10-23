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
from aqmeasy.ui.stylesheets import stylesheets



class QCORR(QWidget):
    def __init__(self):
        super().__init__()
        self.setStyleSheet(stylesheets.QWidget)
        # control.set_parent(self)
        self.setWindowTitle("QCORR")

        # main layout = horizontal box
        main_layout = QHBoxLayout()
        self.setLayout(main_layout)
        # left panel, right panel = vertical box
        left_panel = QVBoxLayout()

        self.file_panel = FilePanel()
        left_panel.addWidget(self.file_panel)
        self.param_panel = ParamPanel()
        left_panel.addWidget(self.param_panel)
        
        right_panel = QVBoxLayout()
        
        main_layout.addLayout(left_panel)
        main_layout.addLayout(right_panel)

