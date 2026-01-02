from collections import UserDict
from PySide6.QtCore import Signal, QObject

class command_model_signals(QObject):
    """Signals for the command model"""
    updated = Signal()

class command_model(UserDict):
    """Replacement for the current dictionary model present in smiles2csv"""

    def __init__(self, *args, **kwargs):
        self.signals = command_model_signals()
        super().__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        prev_value = self.get(key)
        super().__setitem__(key, value)
        if value != prev_value:
            self.signals.updated.emit()

    def __getitem__(self, key):
        return super().__getitem__(key)


general_command_default = {
    "input": None,
    "destination": "",
    "program": "",
    "stacksize": "1GB",
    "sample": 25,
    "auto_sample": "mid",
    "ewin_csearch": 5.0,
    "initial_energy_threshold": 0.0001,
    "energy_threshold": 0.25,
    "rms_threshold": 0.25,
    "opt_steps_rdkit": 1000,
    "heavyonly": True,
    "max_matches_rmsd": 1000,
    "max_mol_wt": 0.1,
    "max_torsions": 0,
    "seed": 62609,
    "geom": [],
    "bond_thres": 0.2,
    "angle_thres": 30,
    "dihedral_thres": 30,
    "auto_metal_atoms": True,
    "complex_type": ""}

general_command_model = command_model(**general_command_default) # ** unpacks the dictionary into keyword arguments

summ_command_default = {
    "degree": 120.0
}

summ_command_model = command_model(**summ_command_default)

fullmonte_command_default = {
    "ewin_fullmonte": 5.0,
    "ewin_sample_fullmonte": 2.0,
    "nsteps_fullmonte": 100,
    "nrot_fullmonte": 3,
    "ang_fullmonte": 30
}

fullmonte_command_model = command_model(**fullmonte_command_default)


crest_command_default = {
    "nprocs": 8,
    "constraints_atoms": [],
    "constraints_dist": [],
    "constraints_angle":[],
    "constraints_dihedral":[],
    "crest_force": 0.5,
    "crest_keywords": None,
    "cregen": True,
    "cregen_keywords": None,
    "xtb_keywords": None,
    "crest_runs": 1
}

crest_command_model = command_model(**crest_command_default)