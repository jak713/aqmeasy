import csv
from PySide6.QtWidgets import QMainWindow, QLabel, QVBoxLayout, QWidget, QPushButton, QHBoxLayout, QMessageBox, QApplication
from PySide6.QtGui import QAction, QPixmap, QDesktopServices
from ui.smiles2csv import smiles_to_csv
from PySide6.QtCore import Qt, QUrl

class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("AQMEasy v1.0")
        self.resize(400, 400)
        self.create_menu()
        self.move(0, 0)


        central_widget = QWidget(self)
        layout = QVBoxLayout(central_widget)
        self.logo = QLabel(self)
        logo = QPixmap("/Users/user/Documents/aqme/aqmeasy/ui/resources/aqmeasy_logo.png")
        if not logo.isNull():
            self.logo.setPixmap(logo)
        else:
            self.logo.setText("Logo not found")
        layout.addWidget(self.logo, alignment=Qt.AlignCenter)

        readthedocs_label = QLabel(self)
        readthedocs_icon = QPixmap("/Users/user/Documents/aqme/aqmeasy/ui/resources/readthedocs_logo.png")
        if not readthedocs_icon.isNull():
            readthedocs_icon = readthedocs_icon.scaled(self.logo.pixmap().size() * 0.5, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            readthedocs_label.setPixmap(readthedocs_icon)
        else:
            readthedocs_label.setText('<a href="https://aqme.readthedocs.io/en/latest/">Docs</a>')
            readthedocs_label.setTextFormat(Qt.RichText)
            readthedocs_label.setTextInteractionFlags(Qt.TextBrowserInteraction)
            readthedocs_label.setOpenExternalLinks(True)
        readthedocs_label.mousePressEvent = lambda event: QDesktopServices.openUrl(QUrl("https://aqme.readthedocs.io/en/latest/"))
        layout.addWidget(readthedocs_label, alignment=Qt.AlignCenter)

        self.setCentralWidget(central_widget)

        self.smiles2csv_list = []
        button_layout = QHBoxLayout()
        button_for_csv = QPushButton("New csv")
        button_for_csv.setFixedHeight(50)
        button_for_csv.clicked.connect(self.new_smiles2csv_window)
        button_layout.addWidget(button_for_csv)
        
        button_for_descriptors = QPushButton("Descriptors")
        button_for_descriptors.setFixedHeight(50)
        button_for_descriptors.clicked.connect(self.new_descriptors_widget)
        button_layout.addWidget(button_for_descriptors)
        
        button_for_analysis = QPushButton("Analysis")
        button_for_analysis.setFixedHeight(50)
        button_for_analysis.clicked.connect(self.new_analysis_widget)
        button_layout.addWidget(button_for_analysis)

        button_for_secret = QPushButton("GoodVibes")
        button_for_secret.setFixedHeight(50)
        button_for_secret.clicked.connect(self.new_secret_widget)
        button_layout.addWidget(button_for_secret)

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
        <h3>Automated Quantum Mechanical Environments</h3>
        <p>Code <a href="https://github.com/jvalegre/aqme">here</a> (Github)</p>
        <p>ReadTheDocs <a href="https://aqme.readthedocs.io/en/latest/">here</a></p>
        <p>Citation: Alegre-Requena JV, Sowndarya S. V. S, Pérez-Soto R, Alturaifi TM, Paton RS. AQME: Automated quantum mechanical environments for researchers and educators. WIREs Comput Mol Sci. 2023; 13(5):e1663. https://doi.org/10.1002/wcms.1663</p>
        """
        QMessageBox.about(self, "About AQMEasy", about_text)

        # <p>GUI developed by <a href="https://github.com/jak713">Julia Kaczmarek</a></p>

    def aboutQt(self):
        QMessageBox.aboutQt(self, QApplication.aboutQt())

    def new_smiles2csv_window(self):
        """Opens a new window for converting SMILES to CSV"""
        new_smiles2csv = smiles_to_csv()
        self.smiles2csv_list.append(new_smiles2csv)
        new_smiles2csv.destroyed.connect(lambda obj: self.smiles2csv_list.remove(new_smiles2csv))
        print(f"Number of smiles2csv windows running: {len(self.smiles2csv_list)}")
        new_smiles2csv.show()

    def new_descriptors_widget(self):
        """Opens a new window for selecting descriptors"""
        pass

    def new_analysis_widget(self):
        """Opens a new window for analysing results"""
        pass

    def new_secret_widget(self):
        """Opens a new window for secret stuff"""
        pass

    def closeEvent(self, event):
        for smiles2csv in self.smiles2csv_list:
            smiles2csv.close()
        event.accept()