import sys

from PySide6.QtWidgets import (QApplication, QWidget, QVBoxLayout, 
                            QHBoxLayout)
from PySide6.QtCore import Qt
from PySide6.QtGui import QCloseEvent
from aqmeasy.ui.QPREP_ui.QPREP_molecularviewer import MoleculeViewer
from aqmeasy.ui.QPREP_ui.QPREP_parameterpanel import ParameterPanel 
from aqmeasy.ui.QPREP_ui.QPREP_filepanel import FilePanel   
from aqmeasy.models.QPREP_model import InputModel
from aqmeasy.ui.stylesheets import stylesheets

class QPREP(QWidget):
    """QPREP Widget opens a window"""
    def __init__(self, parent=None):
        super().__init__()
        if parent is not None:
            self.parent = parent
        self.setStyleSheet(stylesheets.QWidget)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self.setWindowTitle("QPREP")
        # self.resize(1200, 800) 
        # Main vertical layout to hold everything
        main_v_layout = QVBoxLayout()
        
        # Top horizontal layout for the three panels
        top_h_layout = QHBoxLayout()
        
        # Left Panel - Parameters
        left_panel = QVBoxLayout()
        self.input_model = InputModel()
        self.parameter_panel = ParameterPanel(model=self.input_model)
        left_panel.addWidget(self.parameter_panel)
        
        # Middle Panel - File operations
        middle_panel = QVBoxLayout()
        self.molecular_viewer = MoleculeViewer()
        self.file_panel = FilePanel(
            model=self.input_model, 
            parameter_panel=self.parameter_panel,
            molecular_viewer=self.molecular_viewer
        )
        middle_panel.addWidget(self.file_panel)
        
        # Right panel - Molecular Analysis (options part)
        right_panel = QVBoxLayout()
        right_panel.addWidget(self.molecular_viewer.options_container, 1)
        right_panel.addWidget(self.molecular_viewer.viewer_group, 2)
        # Add the three top panels to the top horizontal layout
        top_h_layout.addLayout(left_panel, 1)  
        top_h_layout.addLayout(middle_panel, 1)
        top_h_layout.addLayout(right_panel, 1) 
        
        # Add the top panels and the 3D viewer to the main vertical layout
        main_v_layout.addLayout(top_h_layout, 1)
        
        self.setLayout(main_v_layout)

    def closeEvent(self, event: QCloseEvent) -> None:
        if hasattr(self, 'parent') and self.parent is not None:
            self.parent.button_for_qprep.setEnabled(True)
        return super().closeEvent(event)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = QPREP()
    window.show()
    sys.exit(app.exec())