import logging
from pathlib import Path
from PySide6.QtWidgets import (
    QLabel, QVBoxLayout, QWidget, QPushButton, QHBoxLayout, QTextBrowser,
    QLineEdit, QCheckBox, QMessageBox, QFileDialog, QListWidgetItem
)
from PySide6.QtCore import Qt, Slot, QThread, Signal
from PySide6.QtGui import QCloseEvent, QPixmap, QKeySequence, QShortcut, QIcon, QDoubleValidator

from aqmeasy.ui.stylesheets import stylesheets
from aqmeasy.ui.icons import Icons
from aqmeasy.ui.QPREP_ui.QPREP_molecularviewer import MoleculeViewer
from aqmeasy.ui.CMIN_ui.CMIN_filepanel import FilePanel
from aqmeasy.models.CMIN_model import FileModel
from aqme.cmin import cmin


class CMIN(QWidget):
    def __init__(self,  parent=None):
        super().__init__()
        if parent:
            self.parent = parent
        self.worker = None
        self.file_paths = []
        self.setWindowTitle("CMIN")
        self.setGeometry(150,150,1200, 600)
        self.setStyleSheet(stylesheets.QWidget)

        self.layout = QVBoxLayout()
        view_and_import_layout = QHBoxLayout()

        self.molecule_viewer = MoleculeViewer()
        # self.molecule_viewer.file_group.hide()
        view_and_import_layout.addWidget(self.molecule_viewer)

        self.file_panel = FilePanel(self, model=FileModel())
        self.file_panel.setMaximumWidth(400)
        view_and_import_layout.addWidget(self.file_panel)
        self.layout.addLayout(view_and_import_layout)
        self.setLayout(self.layout)


    
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

    def closeEvent(self, event: QCloseEvent) -> None:
        if hasattr(self, 'parent') and self.parent is not None:
            self.parent.button_for_cmin.setEnabled(True)
        return super().closeEvent(event)