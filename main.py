import sys
from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QIcon
from aqmeasy.ui.main_window import MainWindow
from pathlib import Path
from aqmeasy.ui.stylesheets import stylesheets

if __name__ == "__main__":
    app = QApplication([]) # won't need any args if ran as executable
    base_path = Path(__file__).parent
    icon_path = base_path / "ui" / "resources" / "AQME_icon.png"
    app.setWindowIcon(QIcon(str(icon_path)))
    app.setStyle("Fusion")
    window = MainWindow()
    window.setStyleSheet(stylesheets.QMainWindow)
    window.show()
    window.destroyed.connect(app.quit)
    sys.exit(app.exec())