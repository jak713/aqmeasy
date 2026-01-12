##########################################
#  Stylesheets for AQMEasy UI components #
##########################################

class stylesheets:

    QPushButton = """
QPushButton {
    background-color: #0077B6;
    color: #CAF0F8;
    border-radius: 5px;
    padding: 5px;
}
QPushButton:hover {
    background-color: #00B4D8;
}"""

    QLabelMain = """
QLabel {
    color: #90E0EF;
    padding: 2px;
    font-size: 10pt;
}"""

    QLabel = """
QLabel {
    color: #90E0EF;
    padding: 2px;
    font-size: 12pt; 
}"""

    QMessageBox = """
QMessageBox {
    background-color:#020338;
}

QMessageBox QLabel {
    font-size: 12pt;
    color: #CAF0EF;
}

QMessageBox QComboBox {
    background-color: #020338;
    color: #CAF0F8;
    border: 1px solid #00B4D8;
    border-radius: 3px;
    padding: 2px;
    combobox-popup: 0;
}

QMessageBox QComboBox:hover {
    background-color: #0077B6;
}

QMessageBox QComboBox QAbstractItemView {
    background-color: #020338;
    color: #CAF0F8;
    selection-background-color: #0077B6;
}

QMessageBox QPushButton {
    background-color: #0077B6;
    color: #CAF0F8;
    border: 1px solid #0077B6;
    border-radius: 5px;
    padding: 5px;
    padding-left: 10px;
    padding-right: 10px;
}
QMessageBox QPushButton:hover {
    background-color: #00B4D8;
}"""

    QScrollBar = """
QScrollBar:vertical {
    background: #020338;
    width: 10px;
    margin: 0px 0px 0px 0px;
    border: 1px solid #00B4D8;
    border-radius: 5px;
}

QScrollBar::handle:vertical {
    background: #00B4D8;
    min-height: 10px;
    border-radius: 5px;
}       
QScrollBar::handle:vertical:hover {
    background: #0077B6;
}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0px;
}
QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical {
    background: none;
}
QScrollBar:horizontal {
    background: #020338;
    height: 10px;
    margin: 0px 0px 0px 0px;
    border: 1px solid #00B4D8;
    border-radius: 5px;
}   
QScrollBar::handle:horizontal {
    background: #00B4D8;
    min-width: 10px;
    border-radius: 5px;
}
QScrollBar::handle:horizontal:hover {
    background: #0077B6;
}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
    width: 0px;
}
QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal {
    background: none;
}"""
    
    QMenuBar = """
QMenuBar {
    background-color:#020338;
    color: #90E0EF;
    font-size: 14px;
}
QMenuBar::item {
    color: #CAF0F8;
    padding: 5px 15px;
}

QMenuBar::item:selected {
    color: #FE6253;
    background-color: transparent;
}
QMenu {
    color: #CAF0F8;
    background-color:#0077B6;
}

QMenu::item:selected {
    color: #FE6253;
    background-color: #020338;
}"""

    QComboBox = """
QComboBox {
    background-color: #020338;
    color: #CAF0F8;
    border: 1px solid #00B4D8;
    border-radius: 3px;
    padding: 2px;
    combobox-popup: 0;
}

QComboBox:hover {
    background-color: #0077B6;
}

QComboBox QAbstractItemView {
    background-color: #020338;
    color: #CAF0F8;
    selection-background-color: #0077B6;
}

QComboBox:drop-down {
    subcontrol-origin: padding;
    subcontrol-position: top right;
    width: 15px;
    border-left: 1px solid #00B4D8;
}

QComboBox:disabled {
    color: #000000;
}
"""

    QSpinBox = """
QSpinBox {
    background-color: #020338;
    color: #CAF0F8;
    border: 1px solid #00B4D8;
    border-radius: 3px;
    padding: 2px;
}

QSpinBox:hover {
    background-color: #0077B6;
}

QSpinBox:disabled {
    color: #000000;
}"""

    QDoubleSpinBox = """
QDoubleSpinBox {
    background-color:  #020338;
    color: #CAF0F8;
    border: 1px solid #00B4D8;
    border-radius: 3px;
    padding: 2px;
}

QDoubleSpinBox:hover {
    background-color:#0077B6;
}

QAbstractSpinBox:disabled {
    color: #000000;
}
"""


    QLineEdit = """
QLineEdit {
    background-color: #020338;
    color: #CAF0F8;
    border: 1px solid #00B4D8;
    border-radius: 3px;
    padding: 2px;
}"""

    QTextEdit = """
QTextEdit {
    background-color: #020338;
    color: #CAF0F8;
    border: 1px solid #00B4D8;
    border-radius: 3px;
    padding: 2px;
}

QTextEdit::selection {
    background-color: #00B4D8;
}"""

    QPlainTextEdit = """
QPlainTextEdit {
    background-color: #020338;
    color: #CAF0F8;
    border: 1px solid #00B4D8;
    border-radius: 3px;
    padding: 2px;
}

QPlainTextEdit::selection {
    background-color: #00B4D8;
}"""


    QTextBrowser = """
QTextBrowser {
    background-color: #0077B6;
    color: #CAF0F8;
    border: 1px solid #00B4D8;
    border-radius: 3px;
    padding: 2px;
}
QTextBrowser::selection {
    background: #00B4D8;
}

QScrollBar:vertical {
    background: #0077B6;
}"""

    QCheckBox = """
QCheckBox {
    color: #CAF0EF;
    padding: 2px;
    font-size: 12pt; 
}
QCheckBox::indicator {
    border: 1px solid #CAF0F8;
    border-radius: 3px;
}

QCheckBox::indicator:checked {
    background-color: #41dc8e;
}"""

    QSlider = """
QSlider {
    height: 10px;
    margin: 1px 0;
}

QSlider::groove:horizontal {
    border: 1px solid #00B4D8;
    height: 3px;
    background: #CAF0F8;
    margin: 2px 0;
    border-radius: 4px;
}
QSlider::handle:horizontal {
    background: #00B4D8;
    border: 1px solid #CAF0F8;
    width: 10px;
    margin: -3px 0;
    border-radius: 2px;
}
QSlider::handle:horizontal:hover {
    background: #0077B6;
}

QSlider::sub-page:horizontal {
    background: #0077B6;
    margin: 2px 0;
    border: 1px solid #00B4D8;
    height: 3px;
    border-radius: 4px;
}"""

    QGroupBox = """
QGroupBox {
    border: 1px solid #00B4D8;
    margin-top: 10px;
    border-radius: 5px;
    padding: 5px;
}
QGroupBox:title {
    subcontrol-origin: margin;
    subcontrol-position: top center;
    padding: 0 3px;
    color: #00B4D8;
}
"""

    QListWidget = """
QListWidget { 
    background-color: #020338; 
    color: #CAF0F8; 
    border: 1px solid #00B4D8; 
    border-radius: 3px; 
    padding: 5px; 
}
QListWidget::item { 
    padding: 5px;
}
QListWidget::item:selected { 
    background-color: #90E0EF; 
    color: #020338; 
}"""

    QTreeView = """
QTreeView {
    background-color: #020338;
    color: #CAF0F8;
    border: 1px solid #00B4D8;
}

QTreeView::item {
    padding: 5px;
}

QTreeView::item:selected {
    background-color: #90E0EF;
    color: #020338;
}"""

    QTableWidget = """
QTableWidget {
    background-color: #020338;
    border: 1px solid #00B4D8;
    gridline-color: #00B4D8;
    color: #CAF0F8;
}

QTableWidget::header {
    background-color: #020338;
    color: #CAF0F8;
    padding: 5px;
    border: 1px solid #00B4D8;
}
QTableWidget::header:selected {
    background-color: #90E0EF;
    color: #020338;
}

QTableWidget::item {
    color: #CAF0F8;
    background-color: #020338;
}

QTableWidget::item:selected {
    background-color: #90E0EF;
    color: #020338;
}

QHeaderView::section {
    background-color: #020338;
    color: #CAF0F8;
    padding: 5px;
    border: 1px solid #00B4D8;
}

QHeaderView {
    background-color: #020338;
    border: 1px solid #00B4D8;
}

QTableCornerButton::section {
    background-color: #020338;
    border: 1px solid #00B4D8;
}
""" + QScrollBar

    QWebEngineView = """
QWebEngineView {
    background-color: #020338;
    color: #CAF0EF;
}"""

    QDialog = """
QDialog {
    background-color: #020338;
    color: #CAF0EF;
}"""

    QInputDialog = """
    QInputDialog {
    background-color: #020338;
    color: #CAF0EF;
}
QInputDialog QLabel {
    color: #CAF0EF;
    font-size: 12pt;
}
QInputDialog QComboBox {
    background-color: #020338;
    color: #CAF0F8;
    border: 1px solid #00B4D8;
    border-radius: 3px;
    padding: 2px;
    combobox-popup: 0;
}
QInputDialog QComboBox:hover {
    background-color: #0077B6;
}
QInputDialog QComboBox QAbstractItemView {
    background-color: #020338;
    color: #CAF0F8;
    selection-background-color: #0077B6;
}
QInputDialog QPushButton {
    background-color: #0077B6;
    color: #CAF0F8;
    border: 1px solid #0077B6;
    border-radius: 5px;
    padding: 5px;
    padding-left: 10px;
    padding-right: 10px;
}
QInputDialog QPushButton:hover {
    background-color: #00B4D8;
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
}
QTextBrowser::selection {
    background: #00B4D8;
}"""

    MoleculeLabel = """
QLabel {
    background-color: #f8f8f8;
    border: 1px solid #0077B6; 
    border-radius: 5px;
    color: black; 
    font-size: 12px;
}"""

    MoleculeLabelHover = """
QLabel {
    background-color: #41dc8e; 
    font-size: 12px;
    border-radius: 5px;
    border: 1px solid black; 
    color: black;
}"""

    ToggleButton = """
QPushButton {
    background-color: #020338;
    color: #CAF0F8;
    font-size: 15px;
    font-weight: bold;
}
QPushButton:hover {
    color: #FE6253;
}
"""

    ProgressBar = """
QProgressBar {
    border: 1px solid #00B4D8;
    border-radius: 5px;
    text-align: center;
    background-color: #020338;
    color: #CAF0F8;
}
QProgressBar::chunk {
    background-color: #00B4D8;
    border-radius: 5px;
}"""
##########################################

    QWidget = """
QWidget {
    background-color: #020338;
    color: #CAF0EF;
}

QToolTip {
    background-color: #020338;
    color: #CAF0EF;
    border: 1px solid #0077B6;
    padding: 3px;
    font-size: 12px;
}   
""" + QScrollBar + QLabel + QPushButton + QComboBox + QLineEdit + QTextEdit + QCheckBox + QSlider + QGroupBox + QTableWidget + QListWidget + QTreeView + QSpinBox + QDoubleSpinBox + QDialog + QInputDialog + QWebEngineView + QPlainTextEdit + ProgressBar

    QMainWindow = """
QMainWindow {
    background-color:#020338; 
    color: #90E0EF;
}
""" + QMessageBox + QMenuBar + QPushButton + QLabelMain

    RunButton = """
QPushButton {
    font-size: 14px; 
    color: black; 
    background-color: lightblue;
    padding: 5px;
    height: 35px;
    width: 110px;
    border-radius: 10px;
}

QPushButton:hover {
    background-color: azure;
}"""

    StopButton = """
QPushButton {
    font-size: 14px; 
    color: black; 
    background-color: lightcoral;
    padding: 5px;
    height: 35px;
    max-width: 110px;
    border-radius: 10px;
    alignment: center;
}

QPushButton:hover {
    background-color: darkred;
}
"""
