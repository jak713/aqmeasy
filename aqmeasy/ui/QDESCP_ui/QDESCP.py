from PySide6.QtWidgets import QWidget
from aqmeasy.ui.stylesheets import stylesheets

class QDESCP(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("QDESCP")
        self.setGeometry(200, 200, 800, 600)
        self.setStyleSheet(stylesheets.QWidget)
