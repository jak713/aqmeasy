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


class FileLoadWorker(QObject):
    """Background worker to load file text and associated JSON metadata."""

    finished = Signal(int, str, str, object)
    failed = Signal(int, str, str)
    cancelled = Signal(int, str)

    MAX_DISPLAY_BYTES = 5 * 1024 * 1024
    CHUNK_SIZE = 256 * 1024

    def __init__(self, request_id, file_path):
        super().__init__()
        self.request_id = request_id
        self.file_path = file_path

    @Slot()
    def run(self):
        thread = self.thread()

        try:
            if not self.file_path or not os.path.isfile(self.file_path):
                self.failed.emit(self.request_id, self.file_path, "Selected path is not a valid file.")
                return

            file_size = os.path.getsize(self.file_path)
            bytes_to_read = min(file_size, self.MAX_DISPLAY_BYTES)
            chunks = []
            bytes_read = 0

            with open(self.file_path, 'rb') as f:
                while bytes_read < bytes_to_read:
                    if thread and thread.isInterruptionRequested():
                        self.cancelled.emit(self.request_id, self.file_path)
                        return

                    remaining = bytes_to_read - bytes_read
                    chunk = f.read(min(self.CHUNK_SIZE, remaining))
                    if not chunk:
                        break
                    chunks.append(chunk)
                    bytes_read += len(chunk)

            content = b"".join(chunks).decode('utf-8', errors='replace')

            if file_size > self.MAX_DISPLAY_BYTES:
                content += (
                    "\n\n--- OUTPUT TRUNCATED ---\n"
                    f"Showing first {self.MAX_DISPLAY_BYTES:,} bytes of {file_size:,} bytes.\n"
                    "Open the file externally if you need the full content."
                )

            json_data = self._load_associated_json()

            if thread and thread.isInterruptionRequested():
                self.cancelled.emit(self.request_id, self.file_path)
                return

            self.finished.emit(self.request_id, self.file_path, content, json_data)
        except Exception as e:
            self.failed.emit(self.request_id, self.file_path, str(e))

    def _load_associated_json(self):
        """Load optional JSON metadata associated with an output file."""
        file_dir = os.path.dirname(self.file_path)
        json_dir = os.path.join(file_dir, 'json_files')
        base_name = os.path.splitext(os.path.basename(self.file_path))[0]
        json_path = os.path.join(json_dir, base_name + '.json')

        if not os.path.isfile(json_path):
            return None

        try:
            with open(json_path, 'r', encoding='utf-8', errors='replace') as f:
                return json.load(f)
        except Exception:
            # Do not fail full file loading because metadata JSON is malformed.
            return None

class ViewPanel(QWidget):

    MAX_SEARCH_MATCHES = 1500
    MAX_HIGHLIGHTED_MATCHES = 300

    def __init__(self,model):
        super().__init__()
        self.model = model
        self.controller = ViewController(model, self)
        self.init_ui()
        self.update_ui()
        self.setMinimumWidth(400)

        self.model.currentlySelectedFileChanged.connect(self.on_file_selected)

        self.current_search_query = ""
        self.search_matches = []
        self.search_ranges = []
        self.current_match_index = -1
        self._load_request_id = 0
        self._load_thread = None
        self._load_worker = None
        self._pending_file_path = None
        self._pending_json_data = None

        self._search_timer = QTimer(self)
        self._search_timer.setSingleShot(True)
        self._search_timer.setInterval(180)
        self._search_timer.timeout.connect(self._run_pending_search)
        self._pending_search_query = ""

        self._json_timer = QTimer(self)
        self._json_timer.setSingleShot(True)
        self._json_timer.setInterval(250)
        self._json_timer.timeout.connect(self._apply_pending_json)



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

        search_dialog.search_query_signal.connect(self.schedule_search)
        search_dialog.next_match_signal.connect(self.search_next)
        search_dialog.previous_match_signal.connect(self.search_previous)

        self.search_dialog = search_dialog


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

    @Slot()
    def on_qcorr_started(self):
        """Show an indeterminate progress state while QCORR is running."""
        self.progress_bar.setRange(0, 0)

    @Slot(int, int)
    def on_qcorr_progress(self, value, maximum):
        """Update determinate progress when coarse worker milestones are available."""
        if maximum <= 0:
            self.progress_bar.setRange(0, 0)
            return

        self.progress_bar.setRange(0, maximum)
        self.progress_bar.setValue(min(max(value, 0), maximum))

    @Slot(str)
    def on_qcorr_error(self, message):
        """Keep the last error visible and restore progress bar to idle state."""
        if message:
            self.file_viewer.setPlainText(message)
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(0)

    @Slot()
    def on_qcorr_finished(self):
        """Return the progress bar to an idle completed state."""
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(1)

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
        """Display selected file content asynchronously to keep UI responsive."""
        if not file_path or not os.path.isfile(file_path):
            return

        self._reset_search_state()

        self._load_request_id += 1
        self._pending_file_path = file_path

        if self._load_thread and self._load_thread.isRunning():
            self._load_thread.requestInterruption()
            self.progress_bar.setRange(0, 0)
            return

        self._start_file_load(self._load_request_id, self._pending_file_path)

    def _start_file_load(self, request_id, file_path):
        """Start a background file load for the provided request id and path."""
        if not file_path or not os.path.isfile(file_path):
            return

        self.progress_bar.setRange(0, 0)

        thread = QThread(self)
        worker = FileLoadWorker(request_id, file_path)
        worker.moveToThread(thread)

        thread.started.connect(worker.run)
        worker.finished.connect(self._on_file_loaded)
        worker.failed.connect(self._on_file_load_failed)
        worker.cancelled.connect(self._on_file_load_cancelled)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        worker.cancelled.connect(thread.quit)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(self._on_loader_thread_finished)

        self._load_thread = thread
        self._load_worker = worker
        self._pending_file_path = None
        thread.start()

    @Slot(int, str, str, object)
    def _on_file_loaded(self, request_id, file_path, content, json_data):
        if request_id != self._load_request_id:
            return

        self.file_viewer.setPlainText(content)

        # If user is still switching files, defer expensive JSON panel updates.
        if self._pending_file_path:
            self.progress_bar.setRange(0, 1)
            return

        if json_data is not None:
            self._pending_json_data = json_data
            self._json_timer.start()

        self.progress_bar.setRange(0, 1)

    @Slot(int, str, str)
    def _on_file_load_failed(self, request_id, file_path, error):
        if request_id != self._load_request_id:
            return

        if error == "__cancelled__":
            return

        self._reset_search_state()
        self.file_viewer.setPlainText(f"Error reading file: {error}")
        self.progress_bar.setRange(0, 1)

    @Slot(int, str)
    def _on_file_load_cancelled(self, request_id, file_path):
        if request_id != self._load_request_id:
            return

    @Slot()
    def _on_loader_thread_finished(self):
        if self._load_thread and not self._load_thread.isRunning():
            self._load_thread = None
            self._load_worker = None

        if self._pending_file_path and os.path.isfile(self._pending_file_path):
            self._start_file_load(self._load_request_id, self._pending_file_path)

    @Slot()
    def _apply_pending_json(self):
        data = self._pending_json_data
        self._pending_json_data = None
        if data is None:
            return

        parent = self.parent()
        while parent and not hasattr(parent, 'json_panel'):
            parent = parent.parent()

        if parent and hasattr(parent, 'json_panel'):
            json_panel = getattr(parent, 'json_panel', None)
            if json_panel is not None:
                json_panel.load_json_data(data)

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
                json_panel = getattr(parent, 'json_panel', None)
                if json_panel is not None:
                    json_panel.load_json_data(data)
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
    def schedule_search(self, query):
        """Debounce search to avoid expensive full-document scans on every keystroke."""
        self._pending_search_query = query
        self._search_timer.start()

    @Slot()
    def _run_pending_search(self):
        self.search(self._pending_search_query)

    @Slot(str)
    def search(self, query):
        """
        Search for text in the file viewer and highlight all matches.
        """
        if not query:
            self._reset_search_state()
            if hasattr(self, 'search_dialog'):
                self.search_dialog.update_match_count(0, 0)
            return
        
        self.current_search_query = query
        document = self.file_viewer.document()
        
        self._clear_search_highlighting()
        
        # Find all matches
        self.search_matches = []
        self.search_ranges = []
        cursor = QTextCursor(document)
        
        while True:
            cursor = document.find(query, cursor)
            if cursor.isNull():
                break
            
            self.search_matches.append(cursor.selectionEnd())
            self.search_ranges.append((cursor.selectionStart(), cursor.selectionEnd()))

            if len(self.search_matches) >= self.MAX_SEARCH_MATCHES:
                break
        
        # Highlight first match differently and navigate to it
        if self.search_matches:
            self.current_match_index = 0
            self._highlight_current_match()

        if hasattr(self, 'search_dialog'):
            if self.search_matches:
                self.search_dialog.update_match_count(1, len(self.search_matches))
            else:
                self.search_dialog.update_match_count(0, 0)

    def search_next(self):
        """Navigate to next search match."""
        if not self.search_matches:
            return
        
        self.current_match_index = (self.current_match_index + 1) % len(self.search_matches)
        self._highlight_current_match()

    def search_previous(self):
        """Navigate to previous search match."""
        if not self.search_matches:
            return
        
        self.current_match_index = (self.current_match_index - 1) % len(self.search_matches)
        self._highlight_current_match()

    def _highlight_current_match(self):
        """Highlight the current match with a different color."""
        if not self.search_matches or self.current_match_index < 0:
            return

        self._apply_search_highlights()

        current_start, current_end = self.search_ranges[self.current_match_index]
        cursor = self.file_viewer.textCursor()
        cursor.setPosition(current_start)
        cursor.setPosition(current_end, QTextCursor.MoveMode.KeepAnchor)
        self.file_viewer.setTextCursor(cursor)
        self.file_viewer.ensureCursorVisible()
        
        if hasattr(self, 'search_dialog'):
            self.search_dialog.update_match_count(
            self.current_match_index + 1, 
            len(self.search_matches)
        )

    def _clear_search_highlighting(self):
        """Remove all search highlighting from the text viewer."""
        self.file_viewer.setExtraSelections([])

    def _apply_search_highlights(self):
        """Apply capped visual highlights without mutating document text formatting."""
        if not self.search_ranges:
            self.file_viewer.setExtraSelections([])
            return

        selections = []

        base_format = QTextCharFormat()
        base_format.setBackground(QColor(255, 255, 0, 100))

        current_format = QTextCharFormat()
        current_format.setBackground(QColor(255, 140, 0, 180))

        limit = min(len(self.search_ranges), self.MAX_HIGHLIGHTED_MATCHES)
        for idx in range(limit):
            start, end = self.search_ranges[idx]
            sel = QTextEdit.ExtraSelection()
            sel_cursor = self.file_viewer.textCursor()
            sel_cursor.setPosition(start)
            sel_cursor.setPosition(end, QTextCursor.MoveMode.KeepAnchor)
            sel.cursor = sel_cursor  # type: ignore[attr-defined]
            sel.format = current_format if idx == self.current_match_index else base_format  # type: ignore[attr-defined]
            selections.append(sel)

        self.file_viewer.setExtraSelections(selections)

    def _reset_search_state(self):
        self.current_search_query = ""
        self.search_matches = []
        self.search_ranges = []
        self.current_match_index = -1
        self._clear_search_highlighting()
        
class SearchDialog(QWidget):
    search_query_signal = Signal(str)
    next_match_signal = Signal()
    previous_match_signal = Signal()

    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        layout = QHBoxLayout()
        self.setLayout(layout)
        
        # Search input
        self.search_line = QLineEdit()
        self.search_line.setPlaceholderText("Search...")
        self.search_line.textChanged.connect(self.on_search_text_changed)
        self.search_line.returnPressed.connect(self.on_return_pressed)  # Enter = next match
        layout.addWidget(self.search_line)
        
        # Previous button
        prev_button = QPushButton("↑")
        prev_button.setToolTip("Previous match")
        prev_button.setMaximumWidth(30)
        prev_button.clicked.connect(self.previous_match_signal.emit)
        layout.addWidget(prev_button)
        
        # Next button
        next_button = QPushButton("↓")
        next_button.setToolTip("Next match (Enter)")
        next_button.setMaximumWidth(30)
        next_button.clicked.connect(self.next_match_signal.emit)
        layout.addWidget(next_button)
        
        # Match counter label
        self.match_label = QLabel("")
        self.match_label.setMinimumWidth(60)
        layout.addWidget(self.match_label)

    def on_search_text_changed(self, text):
        """Emit search signal when text changes."""
        self.search_query_signal.emit(text)
    
    def on_return_pressed(self):
        """Enter key = next match."""
        if self.search_line.text():
            self.next_match_signal.emit()
    
    def update_match_count(self, current, total):
        """Update the match counter display."""
        if total > 0:
            self.match_label.setText(f"{current}/{total}")
        else:
            self.match_label.setText("0 matches")
