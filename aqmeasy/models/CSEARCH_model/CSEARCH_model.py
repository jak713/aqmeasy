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

    def __getitem__(self, key):
        return super().__getitem__(key)
    
    def __len__(self) -> int:
        return super().__len__()
    
    def append(self, key, value):
        """Append a value to the list at the given key and emit updated signal."""
        if key in self.data:
            self.data[key].append(value)
            self.signals.updated.emit()
        else:
            raise KeyError(f"Key {key} not found in csv_model.")
    
    def add_row(self, row_data: dict):
        """Add a new row to the model. row_data should be a dictionary with keys matching the model's keys."""
        for key in self.data.keys():
            if key in row_data:
                self.data[key].append(row_data[key])
            else:
                self.data[key].append("")  # Append empty string if key not in row_data
        self.signals.updated.emit()

    def get_row(self, index: int) -> dict:
        """Get a row of data as a dictionary."""
        return {key: self.data[key][index] for key in self.data.keys()}
    
    def get_row_as_list_of_tuples(self, index: int) -> list:
        """Get a row of data as a list of (key, value) tuples."""
        return [(key, self.data[key][index]) for key in self.data.keys()]

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