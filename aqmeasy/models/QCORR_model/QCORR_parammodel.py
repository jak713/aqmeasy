from typing import Any
from PySide6.QtCore import QObject, Signal

class ParamModel(QObject):
    """Signals for the QCORR parameter model"""

    fullcheckChanged = Signal(bool)
    varfileChanged = Signal(str)

    ifreq_cutoffChanged = Signal(float)
    amplitude_ifreqChanged = Signal(float)
    freq_convChanged = Signal(str)

    s2_thresholdChanged = Signal(float)

    nodup_checkChanged = Signal(bool)
    dup_thresholdChanged = Signal(float)
    ro_thresholdChanged = Signal(float)
    isom_typeChanged = Signal(str)
    isom_inputsChanged = Signal(str)
    vdwfracChanged = Signal(float)
    covfracChanged = Signal(float)

    qmInputChanged = Signal(str)

    def __init__(self):
        super().__init__()
        self.fullcheck = True
        self.varfile = None # YAML file input
        
        self.ifreq_cutoff = 0.0 # Ignore -ve frequencies above this
        self.amplitude_ifreq = 0.2 # Scale for imaginary normal mode displacement
        
        self.freq_conv = None # If a string is defined, it will remove calculations that converged during optimization but did not convergence in the subsequent frequency calculation. <<< ask juanvi what this means

        self.s2_threshold = 10.0 # Cutoff for spin contamination in terms of % i.e. ± 10%
        
        self.nodup_check = False # Duplicate filter
        self.dup_threshold = 0.0001 # Hartree, difference in E, H, G for duplicates
        self.ro_threshold = 0.1 # Rotational constant for detecting duplicates

        self.isom_type = None # File extension of initial input file
        self.isom_inputs = None # If none use qcorr defaults

        self.vdwfrac = 0.50 # Fraction of summed VDW radii for detecting bonds in isomerization checks
        self.covfrac = 1.10 # Fraction of summed covalent radii for detecting bonds in isomerization checks

        self.qm_input = ""

    def __setattr__(self, name: str, value: Any) -> None:
        return super().__setattr__(name, value)
    
    def __getattribute__(self, name: str) -> Any:
        return super().__getattribute__(name)

    def update_from_dict(self, params):
        """Updates the model parameters from a given dictionary."""
        for key, value in params.items():
            if hasattr(self, key):
                setattr(self, key, value)
                # Emit the corresponding signal
                signal_name = f"{key}Changed"
                signal = getattr(self, signal_name, None)
                if signal:
                    print(f"Emitting signal for {key} with value {value}")
                    signal.emit(value)

    def as_dict(self):
        """Returns the current state of choices as a dict, and is once again updated upon changing, so that correct user selection is created in the file."""
        return {
            'fullcheck': self.fullcheck,
            'varfile': self.varfile,
            'ifreq_cutoff': self.ifreq_cutoff,
            'amplitude_ifreq': self.amplitude_ifreq,
            'freq_conv': self.freq_conv,
            's2_threshold': self.s2_threshold,
            'nodup_check': self.nodup_check,
            'dup_threshold': self.dup_threshold,
            'ro_threshold': self.ro_threshold,
            'isom_type': self.isom_type,
            'isom_inputs': self.isom_inputs,
            'vdwfrac': self.vdwfrac,
            'covfrac': self.covfrac,
            'qm_input': self.qm_input
        }
    
    def reset_to_defaults(self):
        """Resets all parameters to their default values."""
        defaults = {
            'fullcheck': True,
            'varfile': None,
            'ifreq_cutoff': 0.0,
            'amplitude_ifreq': 0.2,
            'freq_conv': None,
            's2_threshold': 10.0,
            'nodup_check': False,
            'dup_threshold': 0.0001,
            'ro_threshold': 0.1,
            'isom_type': None,
            'isom_inputs': None,
            'vdwfrac': 0.50,
            'covfrac': 1.10,
            'qm_input': ""
        }
        self.update_from_dict(defaults)
    
default_values = {
    'fullcheck': True,
    'varfile': None,
    'ifreq_cutoff': 0.0,
    'amplitude_ifreq': 0.2,
    'freq_conv': None,
    's2_threshold': 10.0,
    'nodup_check': False,
    'dup_threshold': 0.0001,
    'ro_threshold': 0.1,
    'isom_type': None,
    'isom_inputs': None,
    'vdwfrac': 0.50,
    'covfrac': 1.10,
    'qm_input': ""
}