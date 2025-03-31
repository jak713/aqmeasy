import sys
from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QIcon
from ui.main_window import MainWindow


if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon('/Users/user/Documents/aqme/aqmeasy/ui/resources/AQME_icon.png'))
    app.setStyle("Fusion")
    window = MainWindow()
    window.show()
    window.destroyed.connect(app.quit)
    sys.exit(app.exec())