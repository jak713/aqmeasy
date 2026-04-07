import os
import sys
import pytest
from pytestqt import qtbot
from PySide6.QtCore import Qt

from aqmeasy.ui.QPREP_ui.QPREP import QPREP

@pytest.fixture
def qprep_widget(qtbot):
    widget = QPREP()
    widget.hide()  # Prevent the widget from showing during tests
    qtbot.addWidget(widget)
    return widget

class TestOutputResults:
    ...