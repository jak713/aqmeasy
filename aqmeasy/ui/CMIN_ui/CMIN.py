from PySide6.QtWidgets import (
    QVBoxLayout, QWidget, QHBoxLayout, QMessageBox, QStatusBar, QDialog, QPushButton,
    QLabel, QCheckBox, QDialogButtonBox
)
from PySide6.QtCore import Slot, QThread
from PySide6.QtGui import QCloseEvent, QPixmap, QIcon

from aqmeasy.ui.stylesheets import stylesheets
from aqmeasy.ui.icons import Icons
from aqmeasy.ui.QPREP_ui.QPREP_molecularviewer import MoleculeViewer
from aqmeasy.ui.CMIN_ui.CMIN_filepanel import FilePanel
from aqmeasy.ui.CMIN_ui.CMIN_parameters import ParameterPanel
from aqmeasy.ui.CMIN_ui.CMIN_results import ResultsPanel
from aqmeasy.controllers.CMIN_worker import CMINWorker
from aqmeasy.models.CMIN_model import FileModel
from aqmeasy.utils import discover_aqme_result_files


class CMIN(QWidget):
    def __init__(self,  parent=None):
        super().__init__()
        self.parent_window = parent if parent else None
        self.worker = None
        self.worker_thread = None
        self.results_dialog = None
        self.popup_results_panel = None
        self.latest_results = None
        self.file_paths = []
        self.setWindowTitle("CMIN")
        self.setGeometry(150,150,1200, 800)
        self.setStyleSheet(stylesheets.QWidget)

        self.main_layout = QVBoxLayout()
        view_and_import_layout = QHBoxLayout()

        self.molecule_viewer = MoleculeViewer()
        view_and_import_layout.addWidget(self.molecule_viewer)

        self.file_panel = FilePanel(self, model=FileModel())
        self.file_panel.setMaximumWidth(400)
        view_and_import_layout.addWidget(self.file_panel)
        
        self.main_layout.addLayout(view_and_import_layout)

        # Parameters section
        self.parameters_panel = ParameterPanel(self)
        self.main_layout.addWidget(self.parameters_panel)
        
        self.open_results_button = QPushButton("Open CMIN Results")
        self.open_results_button.setEnabled(False)
        self.open_results_button.clicked.connect(self._reopen_results_popup)
        self.main_layout.addWidget(self.open_results_button)
        
        # Status bar
        self.status_bar = QStatusBar()
        self.status_bar.showMessage("Ready")
        self.main_layout.addWidget(self.status_bar)
        
        self.setLayout(self.main_layout)

    def run_cmin(self):
        """Start CMIN execution in background thread"""
        # Validate inputs
        if not self.file_panel.model.files:
            self.failure("No files selected. Please add input files first.")
            return
        
        if not self.file_panel.model.w_dir_main:
            self.failure("No output directory selected. Please select an output directory.")
            return
        
        # Get parameters
        try:
            parameters = self.parameters_panel.get_parameters()
        except ValueError as exc:
            self.failure(str(exc))
            return
        parameters['w_dir_main'] = self.file_panel.model.w_dir_main
        
        # Clear previous results
        self.latest_results = None
        self.open_results_button.setEnabled(False)
        if self.popup_results_panel is not None:
            self.popup_results_panel.clear()
        
        # Create worker and thread
        self.worker = CMINWorker(self.file_panel.model.files, parameters)
        self.worker_thread = QThread()
        self.worker.moveToThread(self.worker_thread)
        
        # Connect signals
        self.worker_thread.started.connect(self.worker.run)
        self.worker.started.connect(self.on_cmin_started)
        self.worker.progress.connect(self.on_cmin_progress)
        self.worker.finished.connect(self.on_cmin_finished)
        self.worker.error.connect(self.on_cmin_error)
        self.worker.finished.connect(self.worker_thread.quit)
        self.worker.error.connect(self.worker_thread.quit)
        self.worker_thread.finished.connect(self.cleanup_thread)
        
        # Start execution
        self.worker_thread.start()
    
    def stop_cmin(self):
        """Stop CMIN execution"""
        if self.worker:
            self.worker.stop()
        if self.worker_thread and self.worker_thread.isRunning():
            self.worker_thread.quit()
            self.worker_thread.wait()
        self.status_bar.showMessage("Stopped by user")
        self.file_panel.set_running(False)
    
    @Slot()
    def on_cmin_started(self):
        """Handle CMIN start"""
        self.file_panel.set_running(True)
        self.status_bar.showMessage("CMIN running...")
    
    @Slot(str)
    def on_cmin_progress(self, message):
        """Handle progress updates"""
        self.status_bar.showMessage(message)
    
    @Slot(dict)
    def on_cmin_finished(self, results):
        """Handle successful completion"""
        self.file_panel.set_running(False)
        self.status_bar.showMessage("CMIN completed successfully")

        self.latest_results = results
        self.open_results_button.setEnabled(True)
        self._show_results_popup(results)
        
        # Color input files by status and add generated output files
        for item in results.get('per_file_data', []):
            input_file = item.get('input_file')
            status = item.get('status', 'success')
            if input_file:
                self.file_panel.set_file_status(input_file, status)

        if results.get('output_files'):
            for output_file in results['output_files']:
                self.file_panel.add_output_file(output_file, 'success')

        warnings = results.get('warnings', [])
        warning_text = ""
        if warnings:
            warning_text = "\nWarnings:\n" + "\n".join(f"- {warning}" for warning in warnings)

        per_file = results.get('per_file_data', [])
        execution_success = sum(1 for item in per_file if item.get('execution_status') == 'success')
        filtered_files = sum(1 for item in per_file if item.get('filter_outcome') == 'filtered')
        retained_all = sum(1 for item in per_file if item.get('filter_outcome') == 'retained_all')
        eliminated = sum(1 for item in per_file if item.get('filter_outcome') == 'eliminated')
        failed_files = results.get('failed_file_count', 0)
        
        # Show success message
        self.success(f"CMIN completed!\n\n"
                    f"Input files: {results.get('input_count', 0)}\n"
                    f"Output files: {results.get('output_count', 0)}\n"
                    f"Execution successful: {execution_success}/{results.get('input_count', 0)}\n"
                    f"Execution failed: {failed_files}\n"
                    f"Retained all conformers: {retained_all}\n"
                    f"Filtered conformers: {filtered_files}\n"
                    f"Eliminated after filtering: {eliminated}\n"
                    f"Filtered to zero files: {results.get('file_level_eliminated', results.get('eliminated_count', 0))}\n"
                    f"Conformer-level eliminated: {results.get('conformer_level_eliminated', 0)}\n"
                    f"Input conformers: {results.get('input_conformer_count', 0)}\n"
                    f"Output conformers: {results.get('output_conformer_count', 0)}"
                    f"{warning_text}")

        self._prompt_follow_up_modules(results)

    def _prompt_follow_up_modules(self, results):
        """Ask whether CMIN results should be opened in QPREP and/or QDESCP."""
        if self.parent_window is None:
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Open Follow-up Modules")
        layout = QVBoxLayout(dialog)
        layout.addWidget(QLabel("Open follow-up modules with CMIN result files:"))

        qprep_box = QCheckBox("Open QPREP")
        qdescp_box = QCheckBox("Open QDESCP")
        layout.addWidget(qprep_box)
        layout.addWidget(qdescp_box)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)

        if dialog.exec() != QDialog.DialogCode.Accepted:
            return

        if not qprep_box.isChecked() and not qdescp_box.isChecked():
            return

        output_dir = results.get('output_dir') or self.file_panel.model.w_dir_main
        cmin_sdf_files = discover_aqme_result_files(
            output_dir,
            source="cmin",
            extensions=(".sdf",),
            recursive=True,
        )

        if qprep_box.isChecked():
            if cmin_sdf_files:
                qprep_widget = self.parent_window.new_qprep_widget()
                qprep_widget.file_panel.get_files_from_csearch(cmin_sdf_files)
            else:
                self.failure("No CMIN SDF result files were found to open in QPREP.")

        if qdescp_box.isChecked():
            if cmin_sdf_files:
                qdescp_widget = self.parent_window.new_qdescp_widget()
                qdescp_widget.set_input_files(cmin_sdf_files)
            else:
                self.failure("No CMIN SDF result files were found to open in QDESCP.")

    def _show_results_popup(self, results):
        """Show results in a dedicated dialog to keep main CMIN window compact."""
        if self.results_dialog is None:
            self.results_dialog = QDialog(self)
            self.results_dialog.setWindowTitle("CMIN Results")
            self.results_dialog.setMinimumSize(1000, 700)
            dialog_layout = QVBoxLayout(self.results_dialog)
            self.popup_results_panel = ResultsPanel(self.results_dialog)
            dialog_layout.addWidget(self.popup_results_panel)

        if self.popup_results_panel is None:
            return

        self.popup_results_panel.clear()
        self.popup_results_panel.update_results(results)
        self.results_dialog.show()
        self.results_dialog.raise_()
        self.results_dialog.activateWindow()

    def _reopen_results_popup(self):
        """Reopen the most recent results popup without rerunning CMIN."""
        if self.latest_results is None:
            return
        self._show_results_popup(self.latest_results)
    
    @Slot(str)
    def on_cmin_error(self, error_message):
        """Handle execution errors"""
        self.file_panel.set_running(False)
        self.status_bar.showMessage("CMIN failed")
        self.failure(error_message)
    
    def cleanup_thread(self):
        """Clean up thread resources"""
        if self.worker_thread:
            self.worker_thread.deleteLater()
            self.worker_thread = None
        if self.worker:
            self.worker.deleteLater()
            self.worker = None

    @Slot(str)
    def success(self, message: str):
        """Show success message"""
        msg = QMessageBox(self)
        msg.setWindowTitle("Success")
        msg.setText(message)
        pixmap = QPixmap(Icons.green)
        if not pixmap.isNull():
            icon = QIcon(pixmap)
            msg.setWindowIcon(icon)
            msg.setIconPixmap(icon.pixmap(64, 64))
        else:
            msg.setIcon(QMessageBox.Icon.Information)
        msg.exec()
    
    @Slot(str)
    def failure(self, message: str):
        """Show failure message"""
        msg = QMessageBox(self)
        msg.setWindowTitle("Error")
        msg.setText(message)
        pixmap = QPixmap(Icons.red)
        if not pixmap.isNull():
            icon = QIcon(pixmap)
            msg.setWindowIcon(icon)
            msg.setIconPixmap(icon.pixmap(64, 64))
        else:
            msg.setIcon(QMessageBox.Icon.Critical)
        msg.exec()

    def closeEvent(self, event: QCloseEvent) -> None:
        # Stop any running operations
        if self.worker_thread and self.worker_thread.isRunning():
            self.stop_cmin()
        
        if self.parent_window is not None and hasattr(self.parent_window, 'button_for_cmin'):
            self.parent_window.button_for_cmin.setEnabled(True)
        return super().closeEvent(event)
