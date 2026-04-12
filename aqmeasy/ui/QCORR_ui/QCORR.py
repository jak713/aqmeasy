from PySide6.QtGui import QCloseEvent
from PySide6.QtWidgets import  QWidget, QVBoxLayout, QHBoxLayout, QPushButton
from PySide6.QtCore import Qt

from aqmeasy.ui.QCORR_ui.QCORR_filepanel import FilePanel
from aqmeasy.ui.QCORR_ui.QCORR_viewpanel import ViewPanel
from aqmeasy.ui.QCORR_ui.QCORR_parampanel import ParamPanel
from aqmeasy.ui.QCORR_ui.QCORR_jsonpanel import JSONpanel

from aqmeasy.models.QCORR_model.QCORR_filemodel import FileModel
from aqmeasy.models.QCORR_model.QCORR_parammodel import ParamModel

from aqmeasy.controllers.QCORR_controllers.QCORR_worker import QCORRWorker

from aqmeasy.ui.stylesheets import stylesheets

class QCORR(QWidget):
    def __init__(self, parent=None):
        super().__init__()
        if parent is not None:
            self.parent = parent
        self.file_model = FileModel()
        self.parameter_model = ParamModel()
        self.worker = QCORRWorker(self)
        self.worker.file_model = self.file_model
        self.worker.result.connect(self.on_qcorr_result)
        
        self.view_panel = ViewPanel(self.file_model)
        self.file_panel = FilePanel(self,self.file_model)

        self.setStyleSheet(stylesheets.QWidget)
        self.setWindowTitle("QCORR")

        # main layout = horizontal box
        main_layout = QHBoxLayout()
        self.setLayout(main_layout)
        panel1 = QVBoxLayout()
        panel2 = QVBoxLayout()
        panel3 = QVBoxLayout()

        panel1.addWidget(self.file_panel)
        self.param_panel = ParamPanel()
        panel1.addWidget(self.param_panel)

        self.toggle_button1 = QPushButton("<")
        self.toggle_button1.setStyleSheet(stylesheets.ToggleButton)
        self.toggle_button1.setCheckable(True)
        self.toggle_button1.toggled.connect(self.toggle_panel1)

        self.toggle_button2 = QPushButton("<")
        self.toggle_button2.setStyleSheet(stylesheets.ToggleButton)
        self.toggle_button2.setCheckable(True)
        self.toggle_button2.toggled.connect(self.toggle_panel2)

        panel2.addWidget(self.view_panel)

        self.json_panel = JSONpanel()
        panel3.addWidget(self.json_panel)

        main_layout.addLayout(panel1)
        main_layout.addWidget(self.toggle_button1, alignment=Qt.AlignmentFlag.AlignTop)
        main_layout.addLayout(panel2)
        main_layout.addWidget(self.toggle_button2, alignment=Qt.AlignmentFlag.AlignTop)
        main_layout.addLayout(panel3)


    def toggle_panel1(self, checked):
        if checked:
            self.file_panel.hide()
            self.param_panel.hide()
            self.toggle_button1.setText(">")
        else:
            self.file_panel.show()
            self.param_panel.show()
            self.toggle_button1.setText("<")

    def toggle_panel2(self, checked):
        if checked:
            self.view_panel.hide()
            self.toggle_button2.setText(">")
        else:
            self.view_panel.show()
            self.toggle_button2.setText("<")

    def run_qcorr(self):
        self.worker.run_qcorr()

    def on_qcorr_result(self, result):
        self.view_panel.file_viewer.setText(result)
        output_dir = self.file_panel.output_dir_label.text()
        if output_dir:
            self.file_model.update_files_after_qcorr(output_dir)

    def closeEvent(self, event: QCloseEvent) -> None:
        self.parent.button_for_qcorr.setEnabled(True)
        return super().closeEvent(event)
    