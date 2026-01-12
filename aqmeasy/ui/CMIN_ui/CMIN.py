from PySide6.QtWidgets import QWidget
from aqmeasy.ui.stylesheets import stylesheets

class CMIN(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("CMIN")
        self.setGeometry(150, 150, 800, 600)
        self.setStyleSheet(stylesheets.QWidget)
