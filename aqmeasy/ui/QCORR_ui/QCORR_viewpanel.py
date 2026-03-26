import os
import json
from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QGroupBox, 
    QFileSystemModel,
    QTreeView,
    QTextBrowser,
    QLineEdit,
    QPushButton,
    QProgressBar,
    QHBoxLayout,
    
)
from PySide6.QtCore import Qt, Signal, Slot
from PySide6.QtGui import QTextOption, QTextCharFormat, QColor, QKeySequence, QShortcut, QTextCursor

from aqmeasy.ui.icons import Icons
from aqmeasy.controllers.QCORR_controllers.QCORR_ViewController import ViewController

class ViewPanel(QWidget):

    def __init__(self,model):
        super().__init__()
        self.model = model
        self.controller = ViewController(model, self)
        self.init_ui()
        self.update_ui()
        self.setMinimumWidth(400)

        self.model.currentlySelectedFileChanged.connect(self.on_file_selected)


    def init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)

        view_group = QGroupBox("File Contents")
        view_layout = QVBoxLayout()
        view_group.setLayout(view_layout)


        search_dialog = SearchDialog()
        view_layout.addWidget(search_dialog)
        search_dialog.setVisible(False)
        search_shortcut = QShortcut(QKeySequence("Ctrl+F"), self)
        search_shortcut.activated.connect(lambda: search_dialog.setVisible(not search_dialog.isVisible()))
        search_dialog.search_query_signal.connect(self.search)


        self.file_viewer = QTextBrowser()
        self.file_viewer.setFontFamily("Menlo")
        self.file_viewer.setWordWrapMode(QTextOption.WrapMode.NoWrap)
        view_layout.addWidget(self.file_viewer)

        self.progress_bar = QProgressBar()
        self.progress_bar.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.progress_bar.setMaximumHeight(10)

        results_group = QGroupBox("QCORR Results")
        results_layout = QVBoxLayout()
        results_group.setLayout(results_layout)

        self.processed_files = QFileSystemModel()
        self.processed_files.setRootPath('')

        self.results_view = QTreeView()
        self.results_view.setModel(self.processed_files)
        self.results_view.setRootIsDecorated(True)
        self.results_view.setColumnWidth(0, 200)

        self.results_view.doubleClicked.connect(
            lambda index: self.display_file_content(
                self.processed_files.filePath(index)
            )
        )
        results_layout.addWidget(self.results_view)
        layout.addWidget(view_group)
        layout.addWidget(self.progress_bar)
        layout.addWidget(results_group)

    def update_ui(self):
        """Update the UI elements based on the model's state."""
        w_dir = self.model.__get__w_dir_main__()
        if w_dir and os.path.exists(w_dir):
            # self.processed_files.setRootPath(w_dir)
            self.results_view.setRootIndex(self.processed_files.index(w_dir))

    @Slot(str)
    def on_file_selected(self, filepath):
        """
        Handle file selection from the file panel.
        Load both the file content and associated JSON.
        """
        if not filepath or not os.path.isfile(filepath):
            return
        
        self.display_file_content(filepath)

    def display_file_content(self, file_path):
        """Display the content of the selected file in the text viewer."""
        try:
            with open(file_path, 'r', encoding='utf-8', errors='replace') as file:
                content = file.read()
                self.file_viewer.setPlainText(content)
        except Exception as e:
            self.file_viewer.setPlainText(f"Error reading file: {e}")
        
        # Try to load associated JSON
        self._try_load_associated_json(file_path)

    def _try_load_associated_json(self, file_path):
        """
        Load JSON metadata for the selected file.
        
        JSON files are always located in a 'json_files' subdirectory
        within the same directory as the output file.
        
        Structure:
            QCORR_X/success/file.log          -> QCORR_X/success/json_files/file.json
            QCORR_X/failed/error/file.log     -> QCORR_X/failed/error/json_files/file.json
        """
        if not file_path or not os.path.isfile(file_path):
            return
        
        # Get directory containing the file
        file_dir = os.path.dirname(file_path)
        
        # JSON is in json_files subdirectory
        json_dir = os.path.join(file_dir, 'json_files')
        
        # Match basename: calc.log -> calc.json
        filename = os.path.basename(file_path)
        base_name = os.path.splitext(filename)[0]
        json_path = os.path.join(json_dir, base_name + '.json')
        
        if os.path.isfile(json_path):
            self._load_json_data(json_path)

    def _load_json_data(self, json_path):
        """Load JSON and pass to the analysis widget."""
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)
            
            # Find the parent QCORR widget and load into json_panel
            parent = self.parent()
            while parent and not hasattr(parent, 'json_panel'):
                parent = parent.parent()
            
            if parent and hasattr(parent, 'json_panel'):
                parent.json_panel.load_json_data(data)
        except Exception as e:
            print(f"Error loading JSON {json_path}: {e}")

        

    def _display_selected_files(self, filenames):
        """Display the list of selected files in the file viewer."""
        display_text = "\n".join(filenames)
        self.file_viewer.setPlainText(display_text)

    def clear_file_viewer(self):
        """Clear the file viewer content."""
        self.file_viewer.clear()

    def refresh_results_tree(self):
        """Force refresh of the file system model to show latest changes."""
        w_dir = self.model.__get__w_dir_main__()
        if w_dir and os.path.exists(w_dir):
            # Force model refresh by resetting root path
            self.processed_files.setRootPath('')
            self.processed_files.setRootPath(w_dir)
            self.results_view.setRootIndex(self.processed_files.index(w_dir))

    @Slot(str)
    def search(self, query):
        ...
        
class SearchDialog(QWidget):
    search_query_signal = Signal(str)

    def __init__(self):
        super().__init__()
        
        layout = QHBoxLayout()
        search_line = QLineEdit()
        search_line.setPlaceholderText("Search...")
        layout.addWidget(search_line)
        self.setLayout(layout)
        # search_line.textChanged.connect(lambda text: self.search_query_signal.emit(text))
        search_line.textChanged.connect(self.search)
        self.search_line = search_line

    def search(self):
        query = self.search_line.text()
        self.search_query_signal.emit(query)