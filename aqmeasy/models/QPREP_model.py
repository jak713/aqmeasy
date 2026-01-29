from PySide6.QtCore import QObject, Signal

# Our MVC. The buttons are all in the default state of how they are when the app is opened, and then if they are changed
# The new values are returned and updated in a dictionary fashion (bottom of code)

class InputModel(QObject):
    """Able to read when button is changed, i.e going from GAUSSIAN -> ORCA ETC"""
    software_Changed = Signal(str)
    functional_Changed = Signal(str)
    basis_set_Changed = Signal(str)
    nprocs_Changed = Signal(int)
    mem_Changed = Signal(int)
    solvent_model_Changed = Signal(str)
    solvent_Changed = Signal(str)
    dispersion_Changed = Signal(str)

    def __init__(self):
        super().__init__()
        self._software = 'Orca'
        self._functional = 'B3LYP'
        self._basis_set = '6-31G(d)'
        self._nprocs = 8
        self._mem = 1
        self._solvent_model = ""
        self._solvent = ""
        self._dispersion = ""

    def check_functional_for_dispersion(self, functional):
        f = (functional or "").strip().lower()
        return ("-d3" in f) or ("-d4" in f) or f.endswith("-v")

    def software(self):
        return self._software
    def set_software(self, value):
        if self._software != value:
            self._software = value
            self.software_Changed.emit(value)

    def functional(self):
        return self._functional
    def set_functional(self, value):
        if self._functional != value:
            self._functional = value
            self.functional_Changed.emit(value)
            if self.check_functional_for_dispersion(value):
                self.set_dispersion("")

    def basis_set(self):
        return self._basis_set
    def set_basis_set(self, value):
        if self._basis_set != value:
            self._basis_set = value
            self.basis_set_Changed.emit(value)

    def nprocs(self):
        return self._nprocs
    def set_nprocs(self, value):
        if self._nprocs != value:
            self._nprocs = value
            self.nprocs_Changed.emit(value)

    def mem(self):
        return self._mem
    def set_mem(self, value):
        if self._mem != value:
            self._mem = value
            self.mem_Changed.emit(value)

    def solvent_model(self):
        return self._solvent_model
    def set_solvent_model(self, value):
        if self._solvent_model != value:
            self._solvent_model = value
            self.solvent_model_Changed.emit(value)

    def solvent(self):
        return self._solvent
    def set_solvent(self, value):
        if self._solvent != value:
            self._solvent = value
            self.solvent_Changed.emit(value)

    def dispersion(self):
        return self._dispersion
    def set_dispersion(self, value):
        if self._dispersion != value:
            self._dispersion = value
            self.dispersion_Changed.emit(value)

    def as_dict(self):
        """Returns the current state of choices as a dict, and is once again updated upon changing, so that correct user selection is created in the file."""
        return {
            'software': self._software,
            'functional': self._functional,
            'basis_set': self._basis_set,
            'nprocs': self._nprocs,
            'mem': self._mem,
            'solvent_model': self._solvent_model,
            'solvent': self._solvent,
            'dispersion': self._dispersion,
        }
    
    def input_as_dict(self, qm_input:dict) -> None:
        """Updates the model attributes from a given qm_input dictionary."""
        for key, value in qm_input.items():
            if hasattr(self, f"set_{key}"):
                setter = getattr(self, f"set_{key}")
                setter(value)
                # Emit the corresponding signal
                signal_name = f"{key}_Changed"
                signal = getattr(self, signal_name, None)
                if signal:
                    print(f"Emitting signal for {key} with value {value}")
                    signal.emit(value)