from PySide6.QtWidgets import QWidget, QVBoxLayout, QPushButton, QFileDialog
from PySide6.QtWebEngineWidgets import QWebEngineView
from aqmeasy.ui.QCORR_ui.QCORR_analysiswidget import QCORRAnalysisWidget
import json


class JSONpanel(QWidget):
    def __init__(self):
        super().__init__()
        self.analysis_widget = None
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # Button to load JSON file
        load_button = QPushButton("Load QCORR JSON")
        load_button.clicked.connect(self.load_json_file)
        layout.addWidget(load_button)
        
        # Analysis widget
        self.analysis_widget = QCORRAnalysisWidget()
        layout.addWidget(self.analysis_widget)
    
    def load_json_file(self):
        """Open file dialog and load JSON data."""
        filepath, _ = QFileDialog.getOpenFileName(
            self,
            "Select QCORR JSON File",
            "",
            "JSON Files (*.json)"
        )
        
        if filepath:
            try:
                with open(filepath, 'r') as f:
                    data = json.load(f)
                self.analysis_widget.load_json_data(data)
            except Exception as e:
                print(f"Error loading JSON: {e}")
    
    def load_json_data(self, data: dict):
        """Programmatically load JSON data (for integration with other parts)."""
        self.analysis_widget.load_json_data(data)