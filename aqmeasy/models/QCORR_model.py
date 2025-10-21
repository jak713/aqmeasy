from PySide6.QtCore import QObject, Signal

# I haven't fully thought out how this is going to work
class QCORRModelSignals(QObject):
    """Signals for the QCORR model"""
    updated = Signal()