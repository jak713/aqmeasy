##########################################
#  Stylesheets for AQMEasy UI components #
##########################################
class stylesheets:

    QMainWindow = """
QMainWindow {
    /* background-color:#818589; */
}"""

    QWidget = """
QWidget {
}"""

    QMenuBar = """
QMenuBar {
    /* background-color:#818589; */
    color: #ffffff;
    font-size: 14px;
}
QMenuBar::item {
    color: #000000;
    padding: 5px 15px;
}

QMenuBar::item:selected {
    color: #ffffff;
    background-color: transparent;
}
QMenu {
    color: #000000;
    /* background-color:#818589; */
}

QMenu::item:selected {
    color: #ffffff;
    background-color: #3c3c3c;
}"""


    QComboBox = """
QComboBox {
    background-color: #edf2f4;
    color: #2b2d42;
    border: 1px solid #8d99ae;
    border-radius: 3px;
    padding: 2px;
    combobox-popup: 0;
}"""

    QSpinBox = """
QSpinBox {
    background-color: #edf2f4;
    color: #2b2d42;
    border: 1px solid #8d99ae;
    border-radius: 3px;
    padding: 2px;
}"""

    QDoubleSpinBox = """
QDoubleSpinBox {
    background-color: #edf2f4;
    color: #2b2d42;
    border: 1px solid #8d99ae;
    border-radius: 3px;
    padding: 2px;
}"""

    QPushButton = """
QPushButton {
    background-color: #3c3c3c;
    color: #ffffff;
    border: 1px solid #5c5c5c;
    border-radius: 5px;
    padding: 5px;
}
QPushButton:hover {
    background-color: #5c5c5c;
}"""

    QLineEdit = """
QLineEdit {
    background-color: #ffffff;
    color: #000000;
    border: 1px solid #000000;
    border-radius: 3px;
    padding: 2px;
}"""

    QTextEdit = """
QTextEdit {
    background-color: #ffffff;
    color: #000000;
    border: 1px solid #000000;
    border-radius: 3px;
    padding: 2px;
}"""

    QTextBrowser = """
QTextBrowser {
    background-color: #ffffff;
    color: #000000;
    border: 1px solid #000000;
    border-radius: 3px;
    padding: 2px;
}"""

    QLabel = """
QLabel {
    padding: 2px;
    color: solid #000000;
    font-size: 11px; 

}"""

    QCheckBox = """
QCheckBox {
    color: #000000;
    padding: 2px;
    font-size: 11px; 
}"""

    QSlider = """
QSlider {
    height: 10px;
    margin: 1px 0;
}

QSlider::groove:horizontal {
    border: 1px solid #999999;
    height: 3px;
    background: #b0b0b0;
    margin: 2px 0;
    border-radius: 4px;
}
QSlider::handle:horizontal {
    background: #5c5c5c;
    border: 1px solid #3c3c3c;
    width: 10px;
    margin: -3px 0;
    border-radius: 2px;
}
QSlider::handle:horizontal:hover {
    background: #7c7c7c;
}

QSlider::sub-page:horizontal {
    background: #3c3c3c;
    margin: 2px 0;
    border: 1px solid #3c3c3c;
    height: 3px;
    border-radius: 4px;
}"""

    QGroupBox = """
QGroupBox {
    border: 1px solid rgba(0, 0, 0, 0.5);
    margin-top: 10px;
    border-radius: 5px;
    padding: 5px;
}
QGroupBox:title {
    subcontrol-origin: margin;
    subcontrol-position: top center;
    padding: 0 3px;
    color: rgba(0, 0, 0, 0.5);
}
"""

    QListWidget = """
QListWidget { 
    background-color: #edf2f4; 
    color: #2b2d42; 
    border: 1px solid #8d99ae; 
    border-radius: 3px; 
    padding: 5px; 
}"""

    QTreeView = """
QTreeView {
    background-color: #2b2d42;
    color: #edf2f4;
    border: 1px solid #8d99ae;
}"""

    QTableWidget = """
QTableWidget {
    border: 1px solid #dcdcdc;
    gridline-color: #dcdcdc;
    font-family: Arial, sans-serif;
    font-size: 12px;
}

QTableWidget::item {
    padding: 5px;
}

QTableWidget::item:selected {
    background-color: #cce7ff;
    color: #000;
}

QTabWidget::pane {
    border: 1px solid #3c3c3c;
    background-color: #2b2d42;
}
QTabBar::tab {
    background-color: #3c3c3c;
    color: #ffffff;
    padding: 5px;
    border: 1px solid #5c5c5c;
    border-bottom: none;
}
QTabBar::tab:selected {
    background-color: #5c5c5c;
}
QTabBar::tab:hover {
    background-color: #4c4c4c;
}"""

    QWebEngineView = """
QWebEngineView {
    border: 1px solid #3c3c3c;
}"""

    ShellOutput = """
QTextBrowser {
    background-color: #1e1e1e;
    color: #dcdcdc;
    border: 1px solid #3c3c3c;
    padding: 10px;
    border-radius: 5px;
    font-family: menlo,"Lucida Console", "Courier New", monospace;
    font-size: 12px;
}"""

    MoleculeLabel = """
QLabel {
    background-color: #f8f8f8;
    border: 1px solid black; 
    border-radius: 5px;
    color: black; 
    font-size: 12px;
}"""

##########################################

# Custom title bar?