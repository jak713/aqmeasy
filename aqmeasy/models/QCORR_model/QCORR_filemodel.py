from PySide6.QtCore import QObject, Signal

class FileModel(QObject):
    """Signals for the QCORR file model"""
    
    filesChanged = Signal(list)
    wdirChanged = Signal(str)
    currentlySelectedFileChanged = Signal(str)

    def __init__(self):
        super().__init__()
        self.files = []
        self.w_dir_main = ""
        self.currently_selected_file = ""

    def __get__files__(self):
        return self.files

    def __get__w_dir_main__(self):
        return self.w_dir_main


    def __get__currently_selected_file__(self):
        return self.currently_selected_file


    # def as_dict(self):
    #     """Returns the current state of choices as a dict, and is once again updated upon changing, so that correct user selection is created in the file."""
    #     return {
    #         'files': self.files,
    #         'w_dir_main': self.w_dir_main,
    #     }