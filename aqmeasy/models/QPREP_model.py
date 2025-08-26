from PySide6.QtCore import QObject, Signal

# Our MVC. The buttons are all in the default state of how they are when the app is opened, and then if they are changed
# The new values are returned and updated in a dictionary fashion (bottom of code)

class InputModel(QObject):
    """Able to read when buttom is changed, i.e going from GAUSSIAN -> ORCA ETC"""
    softwareChanged = Signal(str)
    functionalChanged = Signal(str)
    basisSetChanged = Signal(str)
    nprocsChanged = Signal(int)
    memChanged = Signal(int)

    def __init__(self):
        super().__init__()
        self._software = 'Orca'
        self._functional = 'B3LYP'
        self._basis_set = '6-31G(d)'
        self._nprocs = 8 # just because the default is usually 8
        self._mem = 1

    def software(self):
        return self._software
    def setSoftware(self, value):
        if self._software != value:
            self._software = value
            self.softwareChanged.emit(value)

    def functional(self):
        return self._functional
    def setFunctional(self, value):
        if self._functional != value:
            self._functional = value
            self.functionalChanged.emit(value)

    def basisSet(self):
        return self._basis_set
    def setBasisSet(self, value):
        if self._basis_set != value:
            self._basis_set = value
            self.basisSetChanged.emit(value)

    def nprocs(self):
        return self._nprocs
    def setNprocs(self, value):
        if self._nprocs != value:
            self._nprocs = value
            self.nprocsChanged.emit(value)

    def mem(self):
        return self._mem
    def setMem(self, value):
        if self._mem != value:
            self._mem = value
            self.memChanged.emit(value)

    def as_dict(self):
        """Returns the current state of choices as a dict, and is once again updated upon changing, so that correct user selection is created in the file."""
        return {
            'software': self._software,
            'functional': self._functional,
            'basis_set': self._basis_set,
            'nprocs': self._nprocs,
            'mem': self._mem
        }