from __future__ import annotations

import math

from PySide6.QtCore import QThread, Signal, Slot

from aqme.qdescp import qdescp

from aqmeasy.models.QDESCP_model.aqmetab_model import find_common_smarts_details


class MCSProcessWorker(QThread):
    smartsReady = Signal(str, int, int)
    error = Signal(str)
    timeout = Signal()

    def __init__(self, smiles_list, timeout_ms: int = 30000):
        super().__init__()
        self.smiles_list = list(smiles_list)
        self.timeout_ms = timeout_ms

    @Slot()
    def run(self):
        try:
            timeout_seconds = max(1, math.ceil(self.timeout_ms / 1000))
            result = find_common_smarts_details(self.smiles_list, timeout_seconds=timeout_seconds)
            self.smartsReady.emit(result.smarts, result.atom_count, result.bond_count)
        except TimeoutError:
            self.timeout.emit()
        except Exception as exc:
            self.error.emit(str(exc))


class QDESCPWorker(QThread):
    success_signal = Signal(str)
    failure_signal = Signal(str)
    finished_signal = Signal()

    def __init__(self, params: dict):
        super().__init__()
        self.params = params
        self._stop_requested = False

    def request_stop(self):
        self._stop_requested = True

    @Slot()
    def run(self):
        try:
            if self._stop_requested:
                self.finished_signal.emit()
                return

            qdescp(**self.params)
            self.success_signal.emit("QDESCP finished successfully.")
        except Exception as exc:
            self.failure_signal.emit(str(exc))
        finally:
            self.finished_signal.emit()
