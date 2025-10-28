from PySide6.QtCore import QObject, Signal

class ParamModel(QObject):
    """Signals for the QCORR parameter model"""

    fullcheckChanged = Signal(bool)
    varfileChanged = Signal(str)

    ifreqCutoffChanged = Signal(float)
    amplitudeIfreqChanged = Signal(float)
    freqConvChanged = Signal(str)

    s2ThresholdChanged = Signal(float)

    nodupCheckChanged = Signal(bool)
    dupThresholdChanged = Signal(float)
    roThresholdChanged = Signal(float)

    isomTypeChanged = Signal(str)
    isomInputsChanged = Signal(str)
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

    def update_from_dict(self, params):
        """Updates the model parameters from a given dictionary."""
        for key, value in params.items():
            if hasattr(self, key):
                setattr(self, key, value)

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