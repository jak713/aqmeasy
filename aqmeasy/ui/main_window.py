import csv
from PySide6.QtWidgets import QMainWindow, QLabel, QVBoxLayout, QWidget, QPushButton, QHBoxLayout, QMessageBox, QApplication
from PySide6.QtGui import QAction, QPixmap, QDesktopServices, QIcon
from aqmeasy.ui.CSEARCH_ui.CSEARCH import CSEARCH
from aqmeasy.ui.QPREP_ui.QPREP import QPREP
from aqmeasy.ui.QCORR_ui.QCORR import QCORR
from aqmeasy.ui.stylesheets import stylesheets
from PySide6.QtCore import Qt, QUrl
import os
import requests

class MainWindow(QMainWindow):
#FE6253 coral
#84A1FD light blue
    def __init__(self):
        super().__init__()
        self.setStyleSheet(stylesheets.QMainWindow)
        self.setWindowTitle("AQMEasy 0.1")
        self.create_menu()
        self.move(0, 0)

        #on close destroy all child windows too
        self.setAttribute(Qt.WA_DeleteOnClose)

        central_widget = QWidget(self)
        layout = QVBoxLayout(central_widget)

        logo = QLabel(self)
        logo_path = os.path.join(os.path.dirname(__file__), "resources", "aqme-logo-blue.svg")

        logo_icon = QPixmap(logo_path)
        if not logo_icon.isNull():
            logo.setPixmap(logo_icon)
        else:
            logo.setText("Logo not found")
        logo.setOpenExternalLinks(True)
        logo.mousePressEvent = lambda event: QDesktopServices.openUrl(QUrl("https://github.com/jvalegre/aqme"))
        logo.setCursor(Qt.CursorShape.PointingHandCursor)


        readthedocs_label = QLabel(self)

        readthedocs_label.setOpenExternalLinks(True)
        readthedocs_label.setCursor(Qt.CursorShape.PointingHandCursor)
        readthedocs_path = os.path.join(os.path.dirname(__file__), "resources", "readthedocs_logo.png")

        readthedocs_icon = QPixmap(readthedocs_path)
        if not readthedocs_icon.isNull():
            img = readthedocs_icon.toImage()
            img.invertPixels()
            readthedocs_icon = QPixmap.fromImage(img)
            readthedocs_icon = readthedocs_icon.scaled(logo.pixmap().size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            readthedocs_label.setPixmap(readthedocs_icon)
        else:
            readthedocs_label.setText('<a href="https://aqme.readthedocs.io/en/latest/">Docs</a>')
        
        readthedocs_label.mousePressEvent = lambda event: QDesktopServices.openUrl(QUrl("https://aqme.readthedocs.io/en/latest/"))
        

        layout.addWidget(logo,2, alignment=Qt.AlignCenter)
        layout.addWidget(readthedocs_label,1, alignment=Qt.AlignCenter)
        self.setCentralWidget(central_widget)

        label_layout = QHBoxLayout()
        button_layout = QHBoxLayout()

        label_for_csearch = QLabel("Conformational\nSearch:")
        label_layout.addWidget(label_for_csearch)
        self.button_for_csearch = QPushButton("CSEARCH")
        self.button_for_csearch.setFixedHeight(50)
        self.button_for_csearch.clicked.connect(self.new_csearch_widget)
        button_layout.addWidget(self.button_for_csearch)
        
        label_for_qprep = QLabel("Quantum Chem\nPreprocessing:")
        label_layout.addWidget(label_for_qprep)
        button_for_qprep = QPushButton("QPREP")
        button_for_qprep.setFixedHeight(50)
        button_for_qprep.clicked.connect(self.new_qprep_widget)
        button_layout.addWidget(button_for_qprep)
        
        label_for_qcorr = QLabel("Quantum Chem\nCorrections:")
        label_layout.addWidget(label_for_qcorr)
        button_for_qcorr = QPushButton("QCORR")
        button_for_qcorr.setFixedHeight(50)
        button_for_qcorr.clicked.connect(self.new_qcorr_widget)
        button_layout.addWidget(button_for_qcorr)

        label_for_qdescp = QLabel("Quantum Chem\nDescriptors:")
        label_layout.addWidget(label_for_qdescp)
        button_for_qdescp = QPushButton("QDESCP")
        button_for_qdescp.setFixedHeight(50)
        button_for_qdescp.clicked.connect(self.new_qdescp_widget)
        button_layout.addWidget(button_for_qdescp)
        
        labels = [label_for_csearch, label_for_qprep, label_for_qcorr, label_for_qdescp]
        label_font = labels[0].font()
        label_font.setBold(True)

        for label in labels:
            label.setFont(label_font)
            label.setAlignment(Qt.AlignCenter)

        layout.addStretch()  
        layout.addLayout(label_layout)  
        layout.addLayout(button_layout)


    def create_menu(self):
        menu = self.menuBar()
        menu.setNativeMenuBar(False)

        help_menu = menu.addMenu("Help")
        about_action = QAction("About AQME", self)
        about_action.triggered.connect(self.showAboutDialog)
        help_menu.addAction(about_action)

        about_action = QAction("About Qt", self)
        about_action.triggered.connect(self.aboutQt)
        help_menu.addAction(about_action)

        about_action = QAction("About RDKit", self)
        about_action.triggered.connect(self.aboutRDKit)
        help_menu.addAction(about_action)

        about_action = QAction("About xTB", self)
        about_action.triggered.connect(self.aboutXTB)
        help_menu.addAction(about_action)
        
        about_action = QAction("About CREST", self)
        about_action.triggered.connect(self.aboutCREST)
        help_menu.addAction(about_action)

        about_action = QAction("About cclib", self)
        about_action.triggered.connect(self.aboutcclib)
        help_menu.addAction(about_action)

        about_action = QAction("About py3Dmol", self)
        about_action.triggered.connect(self.aboutpy3Dmol)
        help_menu.addAction(about_action)

        # should one include check for updates?

    def showAboutDialog(self):
        about_text = """
        <div>
            <h3>Automated Quantum Mechanical Environments</h3>
            <p>Code <a href="https://github.com/jvalegre/aqme">here</a> (Github)</p>
            <p>ReadTheDocs <a href="https://aqme.readthedocs.io/en/latest/">here</a></p>
            <p>Citation: Alegre-Requena JV, Sowndarya S. V. S, Pérez-Soto R, Alturaifi TM, Paton RS. AQME: Automated quantum mechanical environments for researchers and educators. WIREs Qut Mol Sci. 2023; 13(5):e1663. https://doi.org/10.1002/wcms.1663</p>
            <p>GUI developed by <a href="https://github.com/jak713">Julia Kaczmarek</a></p>
        </div>
        """
        QMessageBox.about(self, "About AQMEasy", about_text)

    def aboutQt(self):
        QMessageBox.aboutQt(self, QApplication.aboutQt())

    def aboutRDKit(self):
        rdkit_text = """
        <div >
            <h3>RDKit</h3>
            <p>RDKit is an open-source cheminformatics software.</p>
            <p>Find out more at: <a href="https://www.rdkit.org/">https://www.rdkit.org/</a></p>
        </div>
        """
        QMessageBox.about(self, "About RDKit", rdkit_text)

    def aboutXTB(self):
        xtb_text = """
        <div>
            <h3>xTB</h3>
            <p>xTB is a semiempirical tight-binding quantum chemistry software.</p>
            <p>Find out more at: <a href="https://xtb-docs.readthedocs.io/en/latest/">https://xtb-docs.readthedocs.io/en/latest/</a></p>
        </div>
        """
        QMessageBox.about(self, "About xTB", xtb_text)

    def aboutCREST(self):
        crest_text = """
        <div>
            <h3>CREST</h3>
            <p>CREST is a conformer-rotamer ensemble sampling tool based on xTB.</p>
            <p>Find out more at: <a href="https://xtb-docs.readthedocs.io/en/latest/crest.html">https://xtb-docs.readthedocs.io/en/latest/crest.html</a></p>
        </div>
        """
        QMessageBox.about(self, "About CREST", crest_text)

    def aboutcclib(self):
        pass

    def aboutpy3Dmol(self):
        pass

    def new_csearch_widget(self):
        """Opens a new window for converting SMILES to CSV"""
        self.new_csearch = CSEARCH()
        self.new_csearch.show()

    def new_qprep_widget(self):
        self.new_qprep = QPREP()
        self.new_qprep.show()

    def new_qcorr_widget(self):
        self.new_qcorr = QCORR()
        self.new_qcorr.show()

    def new_qdescp_widget(self):
        """Opens a new window for secret stuff"""
        pass

    def closeEvent(self, event):
        event.accept()