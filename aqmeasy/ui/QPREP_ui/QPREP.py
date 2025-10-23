import sys
import os

from PySide6.QtWidgets import (QApplication, QWidget, QVBoxLayout, 
                            QHBoxLayout, QLabel, QFrame, QGroupBox)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from aqmeasy.ui.QPREP_ui.QPREP_molecularviewer import MoleculeViewer
from aqmeasy.ui.QPREP_ui.QPREP_parameterpanel import ParameterPanel 
from aqmeasy.ui.QPREP_ui.QPREP_filepanel import FilePanel   
from aqmeasy.models.QPREP_model import InputModel
from aqmeasy.ui.stylesheets import stylesheets

class QPREP(QWidget):
    """QPREP Widget opens a window"""
    def __init__(self):
        super().__init__()
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
        param_label = QLabel("Calculation Parameters")
        param_label.setStyleSheet(stylesheets.QLabel)
        param_label.setAlignment(Qt.AlignCenter)
        param_label.setObjectName("panelTitle")
        
        self.input_model = InputModel()
        self.parameter_panel = ParameterPanel(model=self.input_model)
        
        left_panel.addWidget(param_label)
        left_panel.addWidget(self.parameter_panel)
        
        # Middle Panel - File operations
        middle_panel = QVBoxLayout()
        file_label = QLabel("File Operations")
        file_label.setStyleSheet(stylesheets.QLabel)
        file_label.setAlignment(Qt.AlignCenter)
        file_label.setObjectName("panelTitle")
        
        self.molecular_viewer = MoleculeViewer()
        self.file_panel = FilePanel(
            model=self.input_model, 
            parameter_panel=self.parameter_panel,
            molecular_viewer=self.molecular_viewer
        )
        
        middle_panel.addWidget(file_label)
        middle_panel.addWidget(self.file_panel)
        
        # Right panel - Molecular Analysis (options part)
        right_panel = QVBoxLayout()
        viewer_label = QLabel("Molecular Analysis")
        viewer_label.setStyleSheet(stylesheets.QLabel)
        viewer_label.setAlignment(Qt.AlignCenter)
        viewer_label.setObjectName("panelTitle")
        
        right_panel.addWidget(viewer_label)
        right_panel.addWidget(self.molecular_viewer.options_container, 1)

        # Add the three top panels to the top horizontal layout
        top_h_layout.addLayout(left_panel, 1)  
        top_h_layout.addLayout(middle_panel, 1)
        top_h_layout.addLayout(right_panel, 1) 
        
        # Add the top panels and the 3D viewer to the main vertical layout
        main_v_layout.addLayout(top_h_layout, 1)
        right_panel.addWidget(self.molecular_viewer.viewer_group, 2)
        
        self.setLayout(main_v_layout)
        

# class EasyChemApp(QApplication):
#     """Main application class that inherits from QApplication"""
    
#     def __init__(self, argv):
#         """Constructor that initializes the application"""
#         super().__init__(argv) 
#         self.setApplicationName("EasyChem")  
#         self.setStyle('Fusion')
        
#         self.setStyleSheet("""
#             QWidget {
#                 background-color: #2b2d42;
#                 color: #edf2f4;
#                 font-family: Arial, sans-serif;
#             }
#             QMainWindow {
#                 background-color: #2b2d42;
#             }
#             QGroupBox {
#                 border: 2px solid #8d99ae;
#                 border-radius: 5px;
#                 margin-top: 2ex;
#                 background-color: #2b2d42;
#             }
#             QGroupBox::title {
#                 subcontrol-origin: margin;
#                 subcontrol-position: top center;
#                 padding: 0 3px;
#                 background-color: #2b2d42;
#                 color: #edf2f4;
#                 font-weight: bold;
#                 font-size: 14px;
#             }
#             QLabel {
#                 color: #edf2f4;
#             }
#             QLabel#panelTitle {
#                 font-size: 20px;
#                 font-weight: bold;
#                 margin: 5px;
#             }
            
#             QPushButton {
#                 background-color: #ef233c;
#                 color: #edf2f4;
#                 border: none;
#                 padding: 6px;
#                 border-radius: 3px;
#                 font-weight: bold;
#             }
#             QPushButton:hover {
#                 background-color: #d90429;
#             }
#             QSlider::groove:horizontal {
#                 border: 1px solid #8d99ae;
#                 height: 8px;
#                 background: #8d99ae;
#                 margin: 2px 0;
#                 border-radius: 4px;
#             }
#             QSlider::handle:horizontal {
#                 background: #ef233c;
#                 border: 1px solid #d90429;
#                 width: 18px;
#                 margin: -2px 0;
#                 border-radius: 9px;
#             }
            
#             QGroupBox > QWidget {
#                 background-color: transparent;
#             }
#         """)