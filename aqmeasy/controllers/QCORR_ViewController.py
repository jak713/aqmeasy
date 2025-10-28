from PySide6.QtCore import QObject, Slot, QRunnable

class ViewController(QObject):
    """Controller for the QCORR module, connecting the model and the view."""
    def __init__(self, model, view):
        super().__init__()
        self.model = model
        self.view = view
        model.currentlySelectedFileChanged.connect(self.thread_display_file)

    def thread_display_file(self):
        thread = Worker()
        file_path = self.model.__get__currently_selected_file__()
        self.view.file_viewer.setText(thread.display_file_contents(file_path))

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