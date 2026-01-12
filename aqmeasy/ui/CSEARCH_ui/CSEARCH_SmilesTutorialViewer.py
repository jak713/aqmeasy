from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWebEngineCore import QWebEnginePage
from PySide6.QtCore import QUrl
from PySide6.QtWidgets import QPushButton
from aqmeasy.ui.stylesheets import stylesheets

class QuietWebPage(QWebEnginePage):
    """Custom page that suppresses mixed content warnings."""
    def javaScriptConsoleMessage(self, level, message, lineNumber, sourceID):
        if "Mixed Content" in message:
            return
        # Pass through other messages
        super().javaScriptConsoleMessage(level, message, lineNumber, sourceID)

class SmilesTutorialViewer(QWebEngineView):
    """A simple web viewer to display the SMILES tutorial."""
    def __init__(self):
        super().__init__()
        self.setPage(QuietWebPage(self))
        self.load(QUrl("https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html"))
        self.setWindowTitle("SMILES Tutorial")
        self.resize(800, 800)

        back_button = QPushButton("<", self)
        back_button.setStyleSheet(stylesheets.ToggleButton)
        back_button.clicked.connect(self.go_back)
        back_button.setFixedSize(30, 30)
        
        forward_button = QPushButton(">", self)
        forward_button.setStyleSheet(stylesheets.ToggleButton)
        forward_button.clicked.connect(self.go_forward)
        forward_button.setFixedSize(30, 30)
        
        back_button.move(10, 10)
        forward_button.move(45, 10)
        
        self.show()

    def go_back(self):
        """Go back to the previous page."""
        if self.history().canGoBack():
            self.back()
        
    def go_forward(self):
        """Go forward to the next page."""
        if self.history().canGoForward():
            self.forward()