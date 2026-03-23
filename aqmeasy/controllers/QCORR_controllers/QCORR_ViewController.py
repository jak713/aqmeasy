from PySide6.QtCore import QObject
from PySide6.QtCore import QThread, Signal

class FileDisplayWorker(QThread):
    """Thread for displaying file contents"""
    contentReady = Signal(str)
    
    def __init__(self, file_path):
        super().__init__()
        self.file_path = file_path
    
    def run(self):
        """Display the contents of the selected file in the text viewer."""
        print(f"Displaying contents of file: {self.file_path}")
        try:
            with open(self.file_path, 'r') as file:
                content = file.read()
                self.contentReady.emit(content)
        except Exception as e:
            self.contentReady.emit(f"Error reading {self.file_path}: {e}")


class ViewController(QObject):
    """Controller for the QCORR module, connecting the model and the view panel."""
    def __init__(self, model, view) :
        super().__init__()
        self.model = model
        self.view = view
        self.worker_thread = None

        self.model.currentlySelectedFileChanged.connect(self.thread_display_file)
        self.model.filesChanged.connect(self.check_clear_file_viewer)
        self.model.wdirChanged.connect(self.view.update_ui)

    def check_clear_file_viewer(self):
        # if model.files == [], clear the file viewer
        if not self.model.__get__files__():
            self.view.file_viewer.clear()
    
    def thread_display_file(self):
        if self.model.isSelectable == False:
            return
        
        file_path = self.model.__get__currently_selected_file__()   
        if file_path != "":
            print(f"Starting thread to display file: {file_path}")
            
            # Clean up previous thread if it exists
            if self.worker_thread is not None:
                self.worker_thread.quit()
                self.worker_thread.wait()
            
            self.worker_thread = FileDisplayWorker(file_path)
            self.worker_thread.contentReady.connect(self.view.file_viewer.setText)
            self.worker_thread.start()