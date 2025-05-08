import sys
from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QIcon
from ui.main_window import MainWindow
from pathlib import Path


if __name__ == "__main__":
    app = QApplication(sys.argv)
    base_path = Path(__file__).parent
    icon_path = base_path / "ui" / "resources" / "AQME_icon.png"
    app.setWindowIcon(QIcon(str(icon_path)))
    app.setStyle("Fusion")
    window = MainWindow()
    window.show()
    window.destroyed.connect(app.quit)
    sys.exit(app.exec())