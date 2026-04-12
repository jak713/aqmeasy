import logging
from pathlib import Path
from PySide6.QtWidgets import (
    QLabel, QVBoxLayout, QWidget, QPushButton, QHBoxLayout, QTextBrowser,
    QLineEdit, QCheckBox, QMessageBox, QSizePolicy, QFileDialog,
    QComboBox,
)
from PySide6.QtCore import Qt, Slot, QThread, Signal
from PySide6.QtGui import QPixmap, QKeySequence, QShortcut, QIcon, QDoubleValidator
from aqme.qdescp import qdescp
from aqmeasy.ui.stylesheets import stylesheets
from aqmeasy.ui.icons import Icons
from aqmeasy.ui.CSEARCH_ui.CSEARCH import CSEARCH
from aqmeasy.ui.QPREP_ui.QPREP_molecularviewer import MoleculeViewer
class QDESCP(QWidget):
    def __init__(self, parent=None):
        super().__init__()
        if parent is not None:
            self.parent = parent
        self.worker = None
        self.file_paths = []
        self.setStyleSheet(stylesheets.QWidget)
        self.setWindowTitle("QDESCP")
        self.set_up_ui()

    def set_up_ui(self):
        """Requirements:
        - Drop in/import for files 
        - CSEARCH button for opening CSEARCH window in a new tab 
        - visualiser of the SDF/xyz/pdb files from embedded py3dmol
        - Program option (xTB/alternative)
        - Boltzmann averaging (True/False checkbox)"""
    
        layout = QVBoxLayout()
        self.setLayout(layout)

        # File input section
        file_input_layout = QHBoxLayout()
        self.file_input = QLineEdit()
        self.file_input.setPlaceholderText("Drag and drop files here or click to browse")
        self.file_input.setReadOnly(True)
        self.file_input.setAcceptDrops(True)
        self.file_input.dragEnterEvent = self.drag_enter_event
        self.file_input.dropEvent = self.drop_event
        browse_button = QPushButton("Browse")
        browse_button.clicked.connect(self.browse_files)
        file_input_layout.addWidget(self.file_input)
        file_input_layout.addWidget(browse_button)
        layout.addLayout(file_input_layout)

        # CSEARCH button
        csearch_button = QPushButton("Work from SMILES via CSEARCH")
        csearch_button.clicked.connect(self.open_csearch)
        layout.addWidget(csearch_button)

        # Visualiser
        self.visualiser = MoleculeViewer()
        layout.addWidget(self.visualiser)

        
        # Program options
        program_options_layout = QHBoxLayout()
        self.program_combo = QComboBox()
        self.program_combo.addItems(["xTB", "alternative"])
        program_options_layout.addWidget(QLabel("Select Program:"))
        program_options_layout.addWidget(self.program_combo)
        layout.addLayout(program_options_layout)

        # Boltzmann averaging checkbox
        self.boltzmann_checkbox = QCheckBox("Boltzmann Averaging")
        layout.addWidget(self.boltzmann_checkbox)  

        # Run button
        run_button = QPushButton("Run QDESCP")
        run_button.clicked.connect(self.run_qdescp)
    
        run_button.setStyleSheet(stylesheets.RunButton)
        layout.addWidget(run_button)

    def drag_enter_event(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
        else:
            event.ignore()

    def drop_event(self, event):
        urls = event.mimeData().urls()
        if urls:
            self.set_input_files([url.toLocalFile() for url in urls])

    def browse_files(self):
        file_dialog = QFileDialog(self)
        file_dialog.setFileMode(QFileDialog.FileMode.ExistingFiles)
        file_dialog.setNameFilters(["SDF files (*.sdf)", "XYZ files (*.xyz)", "PDB files (*.pdb)"])
        if file_dialog.exec():
            self.set_input_files(file_dialog.selectedFiles())

    def set_input_files(self, file_paths):
        """Set QDESCP input files and update UI and viewer."""
        cleaned = [str(Path(path).resolve()) for path in (file_paths or []) if path]
        self.file_paths = cleaned
        self.file_input.setText(", ".join(self.file_paths))
        self.visualiser.load_molecules_from_files(self.file_paths)

    def open_csearch(self):
        """Open CSEARCH new CSEARCH window"""
        self.csearch_window = CSEARCH_window()
        self.csearch_window.show()

    def run_qdescp(self):
        """Run QDESCP in a separate thread to avoid blocking the UI"""
        if not self.file_paths:
            self.failure("Please provide at least one input file.")
            return
        
        program = self.program_combo.currentText()
        boltzmann_averaging = self.boltzmann_checkbox.isChecked()

        ...
    
    
    @Slot(str)
    def success(self, message: str):
        """Show success message"""
        msg = QMessageBox(self)
        msg.setWindowTitle("Success")
        msg.setText(message)
        pixmap = QPixmap(Icons.green)
        if not pixmap.isNull():
            icon = QIcon(pixmap)
            msg.setWindowIcon(icon)
            msg.setIconPixmap(icon.pixmap(64, 64))
        else:
            msg.setIcon(QMessageBox.Icon.Information)
        msg.exec()
    
    @Slot(str)
    def failure(self, message: str):
        """Show failure message"""
        msg = QMessageBox(self)
        msg.setWindowTitle("Error")
        msg.setText(message)
        pixmap = QPixmap(Icons.red)
        if not pixmap.isNull():
            icon = QIcon(pixmap)
            msg.setWindowIcon(icon)
            msg.setIconPixmap(icon.pixmap(64, 64))
        else:
            msg.setIcon(QMessageBox.Icon.Critical)
        msg.exec()

    def closeEvent(self, event):
        if self.worker and self.worker.isRunning():
            self.worker.request_stop()
            self.worker.quit()
            self.worker.wait()
        if hasattr(self, 'parent') and self.parent is not None:
            self.parent.button_for_qdescp.setEnabled(True)
        return super().closeEvent(event)
    

class CSEARCH_window(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CSEARCH for QDESCP")
        self.setStyleSheet(stylesheets.QWidget)
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)

        self.csearch = CSEARCH()
        layout.addWidget(self.csearch)

        close_button = QPushButton("Close")
        close_button.setStyleSheet(stylesheets.QPushButton)
        close_button.clicked.connect(self.close)
        layout.addWidget(close_button)