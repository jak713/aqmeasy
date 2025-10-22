from PySide6.QtCore import QObject, Signal

class FileModel(QObject):
    """Signals for the QCORR file model"""
    
    filesChanged = Signal(list)
    wdirChanged = Signal(str)

    def __init__(self):
        super().__init__()
        self.files = []
        self.w_dir_main = ""


    def update_from_dict(self, params):
        """Updates the model parameters from a given dictionary."""
        for key, value in params.items():
            if hasattr(self, key):
                setattr(self, key, value)

    def as_dict(self):
        """Returns the current state of choices as a dict, and is once again updated upon changing, so that correct user selection is created in the file."""
        return {
            'files': self.files,
            'w_dir_main': self.w_dir_main,
        }