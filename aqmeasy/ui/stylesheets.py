##########################################
#  Stylesheets for AQMEasy UI components #
##########################################
class stylesheets:

    QMainWindow = """
QMainWindow {
    background-color:#03045E; 
    color: #90E0EF;
}
"""

    QWidget = """
QWidget {
    background-color: #03045E;
    color: #CAF0EF;
}

QInputDialog {
    background-color: #03045E;
    color: #CAF0EF;
}

QToolTip {
    background-color: #03045E;
    color: #CAF0EF;
    border: 1px solid #0077B6;
    padding: 3px;
    font-size: 12px;
}   

QScrollBar:vertical {
    background: #03045E;
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
"""

# Placeholder hover style for drag and drop file widgets
#     QWidget_file_hover = """
# QWidget {
#     background-color: #41dc8e; 
# }"""

    QMenuBar = """
QMenuBar {
    background-color:#03045E;
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
    background-color: #03045E;
}"""


    QComboBox = """
QComboBox {
    background-color: #03045E;
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
    background-color: #03045E;
    color: #CAF0F8;
    selection-background-color: #0077B6;
}
"""

    QSpinBox = """
QSpinBox {
    background-color: #03045E;
    color: #CAF0F8;
    border: 1px solid #00B4D8;
    border-radius: 3px;
    padding: 2px;
}

QSpinBox:hover {
    background-color: #0077B6;
}"""

    QDoubleSpinBox = """
QDoubleSpinBox {
    background-color:  #03045E;
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
    QPushButton = """
QPushButton {
    background-color: #0077B6;
    color: #CAF0F8;
    border: 1px solid #0077B6;
    border-radius: 5px;
    padding: 5px;
}
QPushButton:hover {
    background-color: #00B4D8;
}"""

    QLineEdit = """
QLineEdit {
    background-color: #03045E;
    color: #CAF0F8;
    border: 1px solid #00B4D8;
    border-radius: 3px;
    padding: 2px;
}"""

    QTextEdit = """
QTextEdit {
    background-color: #03045E;
    color: #CAF0F8;
    border: 1px solid #00B4D8;
    border-radius: 3px;
    padding: 2px;

}"""

    QTextBrowser = """
QTextBrowser {
    background-color: #0077B6;
    color: #CAF0F8;
    border: 1px solid #00B4D8;
    border-radius: 3px;
    padding: 2px;
}

QScrollBar:vertical {
    background: #0077B6;

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
    border: 1px solid #CAF0F8;
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
    border: 1px solid #CAF0F8;
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
    background-color: #03045E; 
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
    color: #03045E; 
}"""

    QTreeView = """
QTreeView {
    background-color: #03045E;
    color: #CAF0F8;
    border: 1px solid #00B4D8;
}"""

    QTableWidget = """
QTableWidget {
    background-color: #03045E;
    border: 1px solid #00B4D8;
    gridline-color: #00B4D8;
    color: #CAF0F8;
}

QTableWidget::header {
    background-color: #03045E;
    color: #CAF0F8;
    padding: 5px;
    border: 1px solid #00B4D8;
}
QTableWidget::header:selected {
    background-color: #90E0EF;
    color: #03045E;
}

QTableWidget::item {
    color: #CAF0F8;
    background-color: #03045E;
}

QTableWidget::item:selected {
    background-color: #90E0EF;
    color: #03045E;
}

QHeaderView::section {
    background-color: #03045E;
    color: #CAF0F8;
    padding: 5px;
    border: 1px solid #00B4D8;
}

QHeaderView {
    background-color: #03045E;
    border: 1px solid #00B4D8;
}

QTableCornerButton::section {
    background-color: #03045E;
    border: 1px solid #00B4D8;
}

QScrollBar:vertical {
    background: #03045E;
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
    background: #03045E;
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

"""

    QWebEngineView = """
QWebEngineView {
    border: 1px solid #3c3c3c;
}"""

    QDialog = """
QDialog {
    background-color: #03045E;
    color: #CAF0EF;
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

##########################################

# Custom title bar?