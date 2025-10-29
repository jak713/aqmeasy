from PySide6.QtCore import QObject, Slot, QRunnable, QThreadPool

class ViewController(QObject):
    """Controller for the QCORR module, connecting the model and the view panel."""
    def __init__(self, model, view):
        super().__init__()
        self.model = model
        self.view = view
        self.threadpool = QThreadPool()

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
            worker = Worker()
            self.view.file_viewer.setText(worker.display_file_contents(file_path))

class Worker(QRunnable):
    """Thread for displaying text"""

    @Slot()
    def display_file_contents(self, file_path):
        """Display the contents of the selected file in the text viewer."""
        print(f"Displaying contents of file: {file_path}")
        try:
            with open(file_path, 'r') as file:
                content = file.read()
                return content
                    
        except Exception as e:
            return f"Error reading {file_path}: {e}"