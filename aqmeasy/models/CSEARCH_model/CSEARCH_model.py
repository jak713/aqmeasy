from collections import UserDict
from PySide6.QtCore import Signal, QObject

class csv_model_signals(QObject):
    updated = Signal()

class csv_model(UserDict):
    """Replacement for the current dictionary model present in smiles2csv"""

    def __init__(self, *args, **kwargs):
        self.signals = csv_model_signals()
        super().__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        prev_value = self.get(key)
        super().__setitem__(key, value)
        if value != prev_value:
            self.signals.updated.emit()

csv_dictionary = csv_model(
    SMILES = [""],
    code_name = [""],
    charge = [""],
    multiplicity = [""],
    constraints_atoms = [""],
    constraints_dist = [""],
    constraints_angle = [""],
    constraints_dihedral = [""],
    complex_type = [""],
    geom = [""],
)