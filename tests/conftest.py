"""Pytest configuration and shared fixtures for GUI testing."""
import pytest
from unittest.mock import MagicMock, Mock, patch
from PySide6.QtWidgets import QMessageBox, QFileDialog, QInputDialog


@pytest.fixture(scope="session", autouse=True)
def qapp_args():
    """Configure QApplication for headless testing.
    
    This prevents GUI windows from actually rendering on screen during tests.
    """
    return ['--platform', 'offscreen']


@pytest.fixture(autouse=True)
def mock_dialogs(monkeypatch):
    """Mock all Qt dialogs to prevent them from blocking tests.
    
    This fixture automatically applies to all tests and prevents:
    - QMessageBox (warning, information, critical, question dialogs)
    - QFileDialog (file open/save dialogs)
    - QInputDialog (text input dialogs)
    """
    # Create a mock QMessageBox instance that doesn't show
    mock_msgbox_instance = Mock()
    mock_msgbox_instance.exec.return_value = QMessageBox.StandardButton.Ok
    mock_msgbox_instance.exec_.return_value = QMessageBox.StandardButton.Ok
    mock_msgbox_instance.show.return_value = None
    mock_msgbox_instance.setWindowTitle.return_value = None
    mock_msgbox_instance.setText.return_value = None
    mock_msgbox_instance.setWindowIcon.return_value = None
    mock_msgbox_instance.setIconPixmap.return_value = None
    mock_msgbox_instance.setStandardButtons.return_value = None
    mock_msgbox_instance.setIcon.return_value = None
    
    # Mock QMessageBox constructor to return our mock instance
    mock_msgbox_class = MagicMock(return_value=mock_msgbox_instance)
    
    # Mock QMessageBox static methods
    mock_msgbox_class.warning = MagicMock(return_value=QMessageBox.StandardButton.Ok)
    mock_msgbox_class.information = MagicMock(return_value=QMessageBox.StandardButton.Ok)
    mock_msgbox_class.critical = MagicMock(return_value=QMessageBox.StandardButton.Ok)
    mock_msgbox_class.question = MagicMock(return_value=QMessageBox.StandardButton.Yes)
    mock_msgbox_class.about = MagicMock(return_value=None)
    mock_msgbox_class.aboutQt = MagicMock(return_value=None)
    
    # Preserve StandardButton enum
    mock_msgbox_class.StandardButton = QMessageBox.StandardButton
    mock_msgbox_class.Ok = QMessageBox.StandardButton.Ok
    mock_msgbox_class.Save = QMessageBox.StandardButton.Save
    mock_msgbox_class.Discard = QMessageBox.StandardButton.Discard
    mock_msgbox_class.Cancel = QMessageBox.StandardButton.Cancel
    mock_msgbox_class.Yes = QMessageBox.StandardButton.Yes
    mock_msgbox_class.No = QMessageBox.StandardButton.No
    mock_msgbox_class.Information = QMessageBox.Icon.Information
    mock_msgbox_class.Critical = QMessageBox.Icon.Critical
    mock_msgbox_class.Icon = QMessageBox.Icon
    
    # Patch QMessageBox in all modules that import it
    modules_to_patch = [
        "PySide6.QtWidgets.QMessageBox",
        "aqmeasy.controllers.CSEARCH_controller.QMessageBox",
        "aqmeasy.ui.CSEARCH_ui.CSEARCH_widget.QMessageBox",
        "aqmeasy.ui.QPREP_ui.QPREP_molecularviewer.QMessageBox",
        "aqmeasy.ui.main_window.QMessageBox",
    ]
    
    for module in modules_to_patch:
        monkeypatch.setattr(module, mock_msgbox_class)
    
    # Mock QFileDialog methods
    monkeypatch.setattr(QFileDialog, "getSaveFileName", MagicMock(return_value=("", "")))
    monkeypatch.setattr(QFileDialog, "getOpenFileName", MagicMock(return_value=("", "")))
    monkeypatch.setattr(QFileDialog, "getExistingDirectory", MagicMock(return_value=""))
    
    # Also patch in specific modules
    monkeypatch.setattr("aqmeasy.controllers.CSEARCH_controller.QFileDialog.getSaveFileName", MagicMock(return_value=("", "")))
    monkeypatch.setattr("aqmeasy.controllers.CSEARCH_controller.QFileDialog.getOpenFileName", MagicMock(return_value=("", "")))
    
    # Mock QInputDialog methods
    monkeypatch.setattr(QInputDialog, "getText", MagicMock(return_value=("", False)))
