import csv
from PySide6.QtWidgets import QMainWindow, QLabel, QVBoxLayout, QWidget, QPushButton, QHBoxLayout, QMessageBox, QApplication
from PySide6.QtGui import QAction, QPixmap, QDesktopServices
from aqmeasy.ui.CSEARCH import CSEARCH
from aqmeasy.ui.QPREP import QPREP
from PySide6.QtCore import Qt, QUrl
import os

class MainWindow(QMainWindow):
#FE6253 coral
#84A1FD light blue
    def __init__(self):
        super().__init__()
        self.setWindowTitle("aqmeasy 1.0")
        self.resize(400, 500)
        self.create_menu()
        self.move(0, 0)

        central_widget = QWidget(self)
        layout = QVBoxLayout(central_widget)
        self.logo = QLabel(self)
        logo_path = os.path.join(os.path.dirname(__file__), "resources", "aqme-logo-grey-transparent.svg")
        logo = QPixmap(logo_path)
        if not logo.isNull():
            self.logo.setPixmap(logo)
        else:
            self.logo.setText("Logo not found")
        layout.addWidget(self.logo,2, alignment=Qt.AlignCenter)

        readthedocs_label = QLabel(self)
        readthedocs_icon = QPixmap(os.path.join(os.path.dirname(__file__), "resources", "readthedocs_logo.png"))
        if not readthedocs_icon.isNull():
            readthedocs_icon = readthedocs_icon.scaled(self.logo.pixmap().size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            readthedocs_label.setPixmap(readthedocs_icon)
        else:
            readthedocs_label.setText('<a href="https://aqme.readthedocs.io/en/latest/">Docs</a>')
            readthedocs_label.setTextFormat(Qt.RichText)
            readthedocs_label.setTextInteractionFlags(Qt.TextBrowserInteraction)
            readthedocs_label.setOpenExternalLinks(True)
        readthedocs_label.mousePressEvent = lambda event: QDesktopServices.openUrl(QUrl("https://aqme.readthedocs.io/en/latest/"))
        layout.addWidget(readthedocs_label,1, alignment=Qt.AlignCenter)

        self.setCentralWidget(central_widget)

        label_layout = QHBoxLayout()
        button_layout = QHBoxLayout()

        label_for_csearch = QLabel("Conformational\nSearch:")
        label_layout.addWidget(label_for_csearch)
        self.button_for_csearch = QPushButton("CSEARCH")
        self.button_for_csearch.setFixedHeight(50)
        self.button_for_csearch.clicked.connect(self.new_smiles2csv_window)
        button_layout.addWidget(self.button_for_csearch)
        
        label_for_qprep = QLabel("Quantum\nPreprocessing:")
        label_layout.addWidget(label_for_qprep)
        button_for_qprep = QPushButton("QPREP")
        button_for_qprep.setFixedHeight(50)
        button_for_qprep.clicked.connect(self.new_qprep_widget)
        button_layout.addWidget(button_for_qprep)
        
        label_for_qcorr = QLabel("Quantum\nCorrections:")
        label_layout.addWidget(label_for_qcorr)
        button_for_qcorr = QPushButton("QCORR")
        button_for_qcorr.setFixedHeight(50)
        button_for_qcorr.clicked.connect(self.new_qcorr_widget)
        button_layout.addWidget(button_for_qcorr)

        label_for_qdescp = QLabel("Quantum\nDescriptors:")
        label_layout.addWidget(label_for_qdescp)
        button_for_qdescp = QPushButton("QDESCP")
        button_for_qdescp.setFixedHeight(50)
        button_for_qdescp.clicked.connect(self.new_qdescp_widget)
        button_layout.addWidget(button_for_qdescp)
        
        labels = [label_for_csearch, label_for_qprep, label_for_qcorr, label_for_qdescp]
        label_font = labels[0].font()
        label_font.setPointSize(10)
        label_font.setBold(True)

        for label in labels:
            label.setFont(label_font)
            label.setAlignment(Qt.AlignCenter)
            label.setStyleSheet(f"color: {self.palette().color(self.foregroundRole()).name()};")

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

    def showAboutDialog(self):
        about_text = """
        <div style="font-size:10pt;">
            <h3>Automated Quantum Mechanical Environments</h3>
            <p>Code <a href="https://github.com/jvalegre/aqme">here</a> (Github)</p>
            <p>ReadTheDocs <a href="https://aqme.readthedocs.io/en/latest/">here</a></p>
            <p>Citation: Alegre-Requena JV, Sowndarya S. V. S, Pérez-Soto R, Alturaifi TM, Paton RS. AQME: Automated quantum mechanical environments for researchers and educators. WIREs Comput Mol Sci. 2023; 13(5):e1663. https://doi.org/10.1002/wcms.1663</p>
            <p>GUI developed by <a href="https://github.com/jak713">Julia Kaczmarek</a></p>
        </div>
        """
        QMessageBox.about(self, "About AQMEasy", about_text)

        

    def aboutQt(self):
        QMessageBox.aboutQt(self, QApplication.aboutQt())

    def new_smiles2csv_window(self):
        """Opens a new window for converting SMILES to CSV"""
        self.new_smiles2csv = CSEARCH()
        self.new_smiles2csv.show()

    def new_qprep_widget(self):
        self.new_qprep = QPREP()
        print("QPREP created:", self.new_qprep)
        self.new_qprep.show()

    def new_qcorr_widget(self):
        """Opens a new window for qcorr"""
        pass

    def new_qdescp_widget(self):
        """Opens a new window for secret stuff"""
        pass

    def closeEvent(self, event):
        event.accept()