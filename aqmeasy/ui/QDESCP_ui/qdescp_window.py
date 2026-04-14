from __future__ import annotations

import csv
import tempfile
from io import BytesIO
from pathlib import Path
from time import time
from typing import Iterable

from PySide6.QtCore import Qt, Signal, QUrl
from PySide6.QtGui import QDesktopServices, QIcon, QImage, QMouseEvent, QPixmap
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QStatusBar,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdDepictor, rdMolDraw2D

from aqmeasy.controllers.QDESCP_controller import MCSProcessWorker, QDESCPWorker
from aqmeasy.models.QDESCP_model.aqmetab_model import (
    detect_smiles_column,
    generate_mapped_smiles,
    load_molecules_from_path,
    smart_read_csv,
)
from aqmeasy.ui.CSEARCH_ui.CSEARCH import CSEARCH
from aqmeasy.ui.icons import Icons
from aqmeasy.ui.stylesheets import stylesheets


class FileListWidget(QListWidget):
    filesDropped = Signal(list)

    def __init__(self, parent=None, allowed_suffixes: tuple[str, ...] = (".sdf", ".xyz", ".pdb", ".csv", ".cdxml", ".cdx", ".mol")):
        super().__init__(parent)
        self.allowed_suffixes = allowed_suffixes
        self.setAcceptDrops(True)
        self.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
        else:
            event.ignore()

    def dropEvent(self, event):
        paths = [url.toLocalFile() for url in event.mimeData().urls() if url.toLocalFile()]
        accepted = []
        for path in paths:
            path_obj = Path(path)
            if path_obj.is_dir():
                for child in path_obj.rglob("*"):
                    if child.is_file() and child.suffix.lower() in self.allowed_suffixes:
                        accepted.append(str(child))
                continue
            if path_obj.suffix.lower() in self.allowed_suffixes:
                accepted.append(str(path_obj))
        if not accepted:
            event.ignore()
            return
        # Keep order stable while removing duplicates.
        deduped = list(dict.fromkeys(accepted))
        self.filesDropped.emit(deduped)
        event.acceptProposedAction()


class ClickableMoleculeLabel(QLabel):
    clicked = Signal(float, float)

    def mousePressEvent(self, event: QMouseEvent):
        if event.button() == Qt.MouseButton.LeftButton:
            pos = event.position()
            self.clicked.emit(pos.x(), pos.y())
        super().mousePressEvent(event)


class MoleculeTableDialog(QDialog):
    def __init__(self, parent=None, mols: list[Chem.Mol] | None = None, csv_df=None, smiles_column: str | None = None):
        super().__init__(parent)
        self.mols = mols or []
        self.csv_df = csv_df
        self.smiles_column = smiles_column
        self.saved_csv_path: str | None = None
        self.setAcceptDrops(True)
        self.setWindowTitle("Descriptor Input Table")
        self.resize(1000, 700)
        self._build_ui()
        self._populate_rows()

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
        else:
            event.ignore()

    def dropEvent(self, event):
        paths = [url.toLocalFile() for url in event.mimeData().urls() if url.toLocalFile()]
        if not paths:
            event.ignore()
            return
        self._load_source(paths[0])
        event.acceptProposedAction()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["Image", "SMILES", "code_name"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)

        options = QHBoxLayout()
        self.charge_checkbox = QCheckBox("charge")
        self.charge_checkbox.stateChanged.connect(lambda state: self._toggle_column("charge", state))
        self.mult_checkbox = QCheckBox("mult")
        self.mult_checkbox.stateChanged.connect(lambda state: self._toggle_column("mult", state))
        options.addWidget(self.charge_checkbox)
        options.addWidget(self.mult_checkbox)
        options.addStretch(1)
        layout.addLayout(options)

        button_row = QHBoxLayout()
        save_button = QPushButton("Save as CSV")
        save_button.clicked.connect(self._save_to_csv)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_row.addStretch(1)
        button_row.addWidget(save_button)
        button_row.addWidget(cancel_button)
        layout.addLayout(button_row)

    def _populate_rows(self):
        rows: list[tuple[str, str]] = []
        if self.csv_df is not None and self.smiles_column:
            for row_num, (_index, row) in enumerate(self.csv_df.iterrows(), start=1):
                smiles = str(row.get(self.smiles_column, "")).strip()
                if not smiles:
                    continue
                code_name = str(row.get("code_name", f"mol_{row_num}")).strip()
                rows.append((smiles, code_name))
        else:
            for index, mol in enumerate(self.mols, start=1):
                rows.append((Chem.MolToSmiles(mol, canonical=False), f"mol_{index}"))

        self.table.setRowCount(len(rows))
        for row_index, (smiles, code_name) in enumerate(rows):
            self._set_row(row_index, smiles, code_name)

    def _load_source(self, source_path: str):
        source_suffix = Path(source_path).suffix.lower()
        if source_suffix == ".csv":
            csv_df = smart_read_csv(source_path)
            smiles_column = detect_smiles_column(csv_df)
            if smiles_column is None:
                QMessageBox.warning(self, "CSV Error", "No smiles column found in this CSV.")
                return
            self.csv_df = csv_df
            self.smiles_column = smiles_column
            self.mols = []
        else:
            self.mols = load_molecules_from_path(source_path)
            self.csv_df = None
            self.smiles_column = None
            if not self.mols:
                QMessageBox.warning(self, "Input Error", "No valid molecules found in the dropped file.")
                return

        self.table.setRowCount(0)
        self._populate_rows()

    def _set_row(self, row: int, smiles: str, code_name: str):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            image = Draw.MolToImage(mol, size=(100, 100))
            buffer = BytesIO()
            image.save(buffer, format="PNG")
            qimg = QImage.fromData(buffer.getvalue())
            image_label = QLabel()
            image_label.setPixmap(
                QPixmap.fromImage(qimg).scaled(
                    100,
                    100,
                    Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation,
                )
            )
            holder = QWidget()
            holder_layout = QHBoxLayout(holder)
            holder_layout.setContentsMargins(0, 0, 0, 0)
            holder_layout.addWidget(image_label)
            holder_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
            self.table.setCellWidget(row, 0, holder)

        self.table.setItem(row, 1, QTableWidgetItem(smiles))
        self.table.setItem(row, 2, QTableWidgetItem(code_name))

    def _table_headers(self) -> list[str]:
        headers: list[str] = []
        for index in range(self.table.columnCount()):
            item = self.table.horizontalHeaderItem(index)
            headers.append(item.text() if item is not None else "")
        return headers

    def _toggle_column(self, name: str, state: int):
        headers = self._table_headers()
        if state:
            if name in headers:
                return
            new_index = self.table.columnCount()
            self.table.insertColumn(new_index)
            self.table.setHorizontalHeaderItem(new_index, QTableWidgetItem(name))
            for row in range(self.table.rowCount()):
                self.table.setItem(row, new_index, QTableWidgetItem(""))
        else:
            if name not in headers:
                return
            idx = headers.index(name)
            self.table.removeColumn(idx)

    def _save_to_csv(self):
        headers = self._table_headers()
        if "SMILES" not in headers or "code_name" not in headers:
            QMessageBox.warning(self, "Save Error", "SMILES and code_name columns are required.")
            return

        smiles_idx = headers.index("SMILES")
        code_name_idx = headers.index("code_name")
        charge_idx = headers.index("charge") if "charge" in headers else None
        mult_idx = headers.index("mult") if "mult" in headers else None

        code_names = []
        for row in range(self.table.rowCount()):
            smiles_item = self.table.item(row, smiles_idx)
            code_name_item = self.table.item(row, code_name_idx)
            if smiles_item is None or not smiles_item.text().strip():
                QMessageBox.warning(self, "Save Error", "SMILES values cannot be empty.")
                return
            if code_name_item is None or not code_name_item.text().strip():
                QMessageBox.warning(self, "Save Error", "code_name values cannot be empty.")
                return

            code_name = code_name_item.text().strip()
            code_names.append(code_name)

            if charge_idx is not None:
                item = self.table.item(row, charge_idx)
                val = item.text().strip() if item is not None else ""
                if not val:
                    QMessageBox.warning(self, "Save Error", "charge values cannot be empty.")
                    return
                if not val.lstrip("-").isdigit():
                    QMessageBox.warning(self, "Save Error", "charge values must be integers.")
                    return

            if mult_idx is not None:
                item = self.table.item(row, mult_idx)
                val = item.text().strip() if item is not None else ""
                if not val:
                    QMessageBox.warning(self, "Save Error", "mult values cannot be empty.")
                    return
                if not val.lstrip("-").isdigit():
                    QMessageBox.warning(self, "Save Error", "mult values must be integers.")
                    return

        duplicates = [name for name in set(code_names) if code_names.count(name) > 1]
        if duplicates:
            QMessageBox.warning(self, "Save Error", f"Duplicated code_name values: {', '.join(sorted(duplicates))}")
            return

        path, _ = QFileDialog.getSaveFileName(self, "Save CSV", "", "CSV Files (*.csv)")
        if not path:
            return

        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow([header for header in headers if header != "Image"])
            for row in range(self.table.rowCount()):
                values = []
                for col, header in enumerate(headers):
                    if header == "Image":
                        continue
                    item = self.table.item(row, col)
                    values.append(item.text() if item is not None else "")
                writer.writerow(values)

        self.saved_csv_path = path
        self.accept()


class DescriptorsTab(QWidget):
    accepted_drop_suffixes = (".sdf", ".xyz", ".pdb", ".csv", ".cdxml", ".cdx", ".mol")
    min_meaningful_smarts_atoms = 3

    def __init__(self, status_callback=None, open_csearch_callback=None, run_callback=None, parent=None):
        super().__init__(parent)
        self.status_callback = status_callback
        self.open_csearch_callback = open_csearch_callback
        self.run_callback = run_callback

        self.file_paths: list[str] = []
        self.csv_df = None
        self.smiles_column: str | None = None
        self.smarts_pattern: str = ""
        self.selected_atoms: list[int] = []
        self.atom_coords = []
        self.ambiguous_smarts = False
        self.atomic_mapping_enabled = False
        self.atomic_mapping_block_reason: str | None = None
        self.last_output_paths: list[str] = []
        self._mapped_csv_temp_path: str | None = None
        self.mcs_worker: MCSProcessWorker | None = None
        self._draw_width = 900
        self._draw_height = 440
        self.setAcceptDrops(True)

        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)

        columns_layout = QHBoxLayout()
        layout.addLayout(columns_layout, 1)

        # Left column: SMARTS + viewer
        left_col = QVBoxLayout()
        columns_layout.addLayout(left_col, 3)

        smarts_group = QGroupBox("SMARTS / Atomic Selection")
        smarts_layout = QVBoxLayout(smarts_group)

        pattern_row = QHBoxLayout()
        self.smarts_input = QLineEdit()
        self.smarts_input.setPlaceholderText("SMARTS pattern (optional)")
        detect_btn = QPushButton("Detect SMARTS")
        detect_btn.clicked.connect(self.detect_smarts)
        apply_btn = QPushButton("Apply SMARTS")
        apply_btn.clicked.connect(self.apply_manual_smarts)
        pattern_row.addWidget(self.smarts_input, 1)
        pattern_row.addWidget(detect_btn)
        pattern_row.addWidget(apply_btn)
        smarts_layout.addLayout(pattern_row)

        self.mol_viewer = ClickableMoleculeLabel()
        self.mol_viewer.setMinimumSize(560, 320)
        self.mol_viewer.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.mol_viewer.setStyleSheet("border: 1px solid #aaa; background: #f9f9f9;")
        self.mol_viewer.setText("Detect or apply a SMARTS pattern to enable clickable atom selection.")
        self.mol_viewer.clicked.connect(self._on_viewer_clicked)
        smarts_layout.addWidget(self.mol_viewer, 1)

        self.mol_info_label = QLabel("Atomic descriptors are optional. If none are selected, molecular descriptors only will be generated.")
        self.mol_info_label.setWordWrap(True)
        smarts_layout.addWidget(self.mol_info_label)

        self.smarts_guidance_label = QLabel(
            "Use SMARTS only when the input set shares a clear scaffold. Mixed or diverse datasets should run molecular descriptors only."
        )
        self.smarts_guidance_label.setWordWrap(True)
        self.smarts_guidance_label.setStyleSheet("color: #666; font-style: italic;")
        smarts_layout.addWidget(self.smarts_guidance_label)

        left_col.addWidget(smarts_group, 1)

        # Right column: files + settings + run + outputs drawer
        right_col = QVBoxLayout()
        columns_layout.addLayout(right_col, 2)

        input_group = QGroupBox("Input")
        input_layout = QVBoxLayout(input_group)

        self.file_list = FileListWidget(self, allowed_suffixes=self.accepted_drop_suffixes)
        self.file_list.setMinimumHeight(180)
        self.file_list.filesDropped.connect(self.add_files)
        input_layout.addWidget(self.file_list)

        file_buttons = QHBoxLayout()
        add_files_btn = QPushButton("Add Files...")
        add_files_btn.clicked.connect(self.browse_files)
        remove_files_btn = QPushButton("Remove Selected")
        remove_files_btn.clicked.connect(self.remove_selected_files)
        file_buttons.addWidget(add_files_btn)
        file_buttons.addWidget(remove_files_btn)
        input_layout.addLayout(file_buttons)

        flow_buttons = QHBoxLayout()
        start_csearch_btn = QPushButton("Start from SMILES / ChemDraw")
        start_csearch_btn.clicked.connect(self._open_csearch)
        open_table_btn = QPushButton("Open Molecule Table")
        open_table_btn.clicked.connect(self.open_molecule_table)
        flow_buttons.addWidget(start_csearch_btn)
        flow_buttons.addWidget(open_table_btn)
        input_layout.addLayout(flow_buttons)

        self.mode_label = QLabel("Add SDF, XYZ, PDB, or CSV files. Use CSEARCH if starting from SMILES/ChemDraw.")
        self.mode_label.setWordWrap(True)
        input_layout.addWidget(self.mode_label)
        right_col.addWidget(input_group)

        settings_group = QGroupBox("QDESCP Settings")
        settings_layout = QFormLayout(settings_group)
        self.charge_spin = QSpinBox()
        self.charge_spin.setRange(-99, 99)
        self.charge_spin.setValue(0)
        self.mult_spin = QSpinBox()
        self.mult_spin.setRange(1, 99)
        self.mult_spin.setValue(1)
        self.atoms_input = QLineEdit()
        self.atoms_input.setPlaceholderText("e.g. P, C=O")
        self.boltzmann_checkbox = QCheckBox("Boltzmann Averaging")
        self.boltzmann_checkbox.setChecked(True)
        self.robert_checkbox = QCheckBox("Robert-ready")
        settings_layout.addRow("Charge", self.charge_spin)
        settings_layout.addRow("Multiplicity", self.mult_spin)
        settings_layout.addRow("QDESCP Atoms", self.atoms_input)
        settings_layout.addRow("", self.boltzmann_checkbox)
        settings_layout.addRow("", self.robert_checkbox)
        right_col.addWidget(settings_group)

        run_row = QHBoxLayout()
        run_row.addStretch(1)
        self.run_button = QPushButton("Run QDESCP")
        self.run_button.setObjectName("run_descriptors_btn")
        self.run_button.setStyleSheet(stylesheets.RunButton)
        self.run_button.setFixedHeight(36)
        self.run_button.setMaximumWidth(180)
        self.run_button.clicked.connect(self._handle_run)
        run_row.addWidget(self.run_button)
        right_col.addLayout(run_row)

        self.outputs_toggle_btn = QPushButton("Show Outputs")
        self.outputs_toggle_btn.clicked.connect(self.toggle_output_drawer)
        right_col.addWidget(self.outputs_toggle_btn)

        self.output_group = QGroupBox("Outputs")
        self.output_group.setVisible(False)
        output_layout = QVBoxLayout(self.output_group)
        self.output_list = QListWidget()
        self.output_list.setMinimumHeight(100)
        output_layout.addWidget(self.output_list)
        output_btn_row = QHBoxLayout()
        self.open_output_folder_btn = QPushButton("Open Output Folder")
        self.open_output_folder_btn.setEnabled(False)
        self.open_output_folder_btn.clicked.connect(self.open_output_folder)
        output_btn_row.addStretch(1)
        output_btn_row.addWidget(self.open_output_folder_btn)
        output_layout.addLayout(output_btn_row)
        right_col.addWidget(self.output_group)

        right_col.addStretch(1)

    def _set_status(self, message: str):
        if self.status_callback is not None:
            self.status_callback(message)

    def shutdown_workers(self):
        if self.mcs_worker is not None and self.mcs_worker.isRunning():
            self.mcs_worker.requestInterruption()
            self.mcs_worker.wait(2000)

    def _open_csearch(self):
        if self.open_csearch_callback is not None:
            self.open_csearch_callback()

    def toggle_output_drawer(self):
        visible = not self.output_group.isVisible()
        self.output_group.setVisible(visible)
        self.outputs_toggle_btn.setText("Hide Outputs" if visible else "Show Outputs")

    def add_files(self, file_paths: Iterable[str]):
        added_count = 0
        for file_path in file_paths:
            if not Path(file_path).exists():
                continue
            if file_path in self.file_paths:
                continue
            self.file_paths.append(file_path)
            self.file_list.addItem(QListWidgetItem(file_path))
            added_count += 1
        self._update_mode_label()
        self._refresh_csv_context()
        if added_count:
            self._set_status(f"Added {added_count} file(s).")

    def browse_files(self):
        file_paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Select input files",
            "",
            "QDESCP Input (*.sdf *.xyz *.pdb *.csv *.cdxml *.cdx *.mol)",
        )
        if file_paths:
            self.add_files(file_paths)

    def dragEnterEvent(self, event):
        if not event.mimeData().hasUrls():
            event.ignore()
            return
        accepted = self._accepted_file_paths(event.mimeData().urls())
        if accepted:
            event.acceptProposedAction()
        else:
            event.ignore()

    def dropEvent(self, event):
        accepted = self._accepted_file_paths(event.mimeData().urls())
        if not accepted:
            event.ignore()
            return
        self.add_files(accepted)
        event.acceptProposedAction()

    def _accepted_file_paths(self, urls) -> list[str]:
        accepted = []
        for url in urls:
            path = url.toLocalFile()
            if not path:
                continue
            path_obj = Path(path)
            if path_obj.is_dir():
                for child in path_obj.rglob("*"):
                    if child.is_file() and child.suffix.lower() in self.accepted_drop_suffixes:
                        accepted.append(str(child))
                continue
            if path_obj.suffix.lower() in self.accepted_drop_suffixes:
                accepted.append(str(path_obj))
        return list(dict.fromkeys(accepted))

    def remove_selected_files(self):
        for item in self.file_list.selectedItems():
            row = self.file_list.row(item)
            self.file_list.takeItem(row)
            self.file_paths.pop(row)
        self._update_mode_label()
        self._refresh_csv_context()

    def _update_mode_label(self):
        if not self.file_paths:
            self.mode_label.setText("Add SDF, XYZ, PDB, or CSV files. Use CSEARCH if starting from SMILES/ChemDraw.")
            return
        if any(path.lower().endswith(".csv") for path in self.file_paths):
            self.mode_label.setText("CSV mode detected. You can detect SMARTS and click atoms for atomic descriptor mapping.")
        else:
            self.mode_label.setText("Direct structure mode detected. SMARTS mapping requires a CSV with a smiles column.")

    def _refresh_csv_context(self):
        self.csv_df = None
        self.smiles_column = None
        for path in self.file_paths:
            if path.lower().endswith(".csv"):
                try:
                    df = smart_read_csv(path)
                except Exception:
                    continue
                smiles_col = detect_smiles_column(df)
                if smiles_col:
                    self.csv_df = df
                    self.smiles_column = smiles_col
                    return

    def get_files(self) -> list[str]:
        return list(self.file_paths)

    def set_files(self, file_paths: Iterable[str]):
        self.file_paths = []
        self.file_list.clear()
        self.add_files(file_paths)

    def set_prefill(self, charge: int | None = None, mult: int | None = None):
        if charge is not None:
            self.charge_spin.setValue(charge)
        if mult is not None:
            self.mult_spin.setValue(mult)

    def open_molecule_table(self):
        source, _ = QFileDialog.getOpenFileName(
            self,
            "Open molecule source",
            "",
            "Molecule Input (*.csv *.sdf *.cdxml *.cdx *.mol)",
        )
        if not source:
            return

        try:
            if source.lower().endswith(".csv"):
                csv_df = smart_read_csv(source)
                smiles_column = detect_smiles_column(csv_df)
                if smiles_column is None:
                    QMessageBox.warning(self, "CSV Error", "No smiles column found in this CSV.")
                    return
                dialog = MoleculeTableDialog(self, csv_df=csv_df, smiles_column=smiles_column)
            else:
                mols = load_molecules_from_path(source)
                if not mols:
                    QMessageBox.warning(self, "Input Error", "No valid molecules found in the selected file.")
                    return
                dialog = MoleculeTableDialog(self, mols=mols)

            if dialog.exec() and dialog.saved_csv_path:
                self.add_files([dialog.saved_csv_path])
        except Exception as exc:
            QMessageBox.critical(self, "Import Error", str(exc))

    def detect_smarts(self):
        if self.csv_df is None or self.smiles_column is None:
            QMessageBox.warning(self, "SMARTS", "Load a CSV with a smiles column to detect SMARTS.")
            return

        smiles_values = self.csv_df[self.smiles_column].dropna().astype(str).tolist()
        if not smiles_values:
            QMessageBox.warning(self, "SMARTS", "No smiles entries found in the CSV.")
            return

        self.mcs_worker = MCSProcessWorker(smiles_values, timeout_ms=30000)
        self.mcs_worker.smartsReady.connect(self._on_smarts_detected)
        self.mcs_worker.error.connect(self._on_smarts_error)
        self.mcs_worker.timeout.connect(self._on_smarts_timeout)
        self.mcs_worker.start()
        self._set_status("Detecting common SMARTS pattern...")

    def _on_smarts_detected(self, smarts: str, atom_count: int, bond_count: int):
        self.smarts_input.setText(smarts)
        if atom_count < self.min_meaningful_smarts_atoms:
            self._activate_smarts(
                smarts,
                allow_atomic_selection=False,
                block_reason=(
                    f"Only a {atom_count}-atom shared core was found. Atomic descriptors are likely not useful for this set."
                ),
            )
            self._set_status("SMARTS detected, but the shared core is too small for reliable atomic descriptors.")
            return

        self._activate_smarts(smarts)
        self._set_status(f"SMARTS detected with a {atom_count}-atom shared core.")

    def _on_smarts_error(self, message: str):
        QMessageBox.critical(self, "SMARTS Error", message)
        self._set_status("SMARTS detection failed.")

    def _on_smarts_timeout(self):
        QMessageBox.warning(self, "SMARTS Timeout", "MCS detection timed out. Try a smaller input set or manual SMARTS.")
        self._set_status("SMARTS detection timed out.")

    def apply_manual_smarts(self):
        pattern = self.smarts_input.text().strip()
        if not pattern:
            QMessageBox.warning(self, "SMARTS", "Enter a SMARTS pattern first.")
            return
        self._activate_smarts(pattern)

    def _activate_smarts(self, pattern: str, allow_atomic_selection: bool = True, block_reason: str | None = None):
        mol = Chem.MolFromSmarts(pattern)
        if mol is None:
            QMessageBox.warning(self, "SMARTS", "Invalid SMARTS pattern.")
            return

        self.smarts_pattern = pattern
        self.selected_atoms = []
        self.ambiguous_smarts = self._has_ambiguous_matches(pattern)
        self.atomic_mapping_enabled = allow_atomic_selection and not self.ambiguous_smarts
        self.atomic_mapping_block_reason = block_reason

        if self.ambiguous_smarts:
            self.mol_info_label.setText(
                "Multiple SMARTS matches detected in at least one molecule. Atomic mapping will be blocked for this run."
            )
        elif not self.atomic_mapping_enabled and self.atomic_mapping_block_reason:
            self.mol_info_label.setText(self.atomic_mapping_block_reason + " Molecular descriptors will still run.")
        else:
            self.mol_info_label.setText("SMARTS active. Click atoms in the scaffold image to select atomic descriptors.")
        self._render_smarts_molecule(mol)

    def _has_ambiguous_matches(self, pattern: str) -> bool:
        if self.csv_df is None or self.smiles_column is None:
            return False

        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol is None:
            return False

        for smiles in self.csv_df[self.smiles_column].dropna().astype(str):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            matches = mol.GetSubstructMatches(pattern_mol)
            if len(matches) > 1:
                return True
        return False

    def _render_smarts_molecule(self, smarts_mol: Chem.Mol):
        rdDepictor.SetPreferCoordGen(True)
        drawer = rdMolDraw2D.MolDraw2DCairo(self._draw_width, self._draw_height)
        highlight = {idx: (0.13, 0.45, 0.95) for idx in self.selected_atoms}
        drawer.DrawMolecule(smarts_mol, highlightAtoms=self.selected_atoms, highlightAtomColors=highlight)
        drawer.FinishDrawing()

        png_bytes = drawer.GetDrawingText()
        pixmap = QPixmap()
        pixmap.loadFromData(png_bytes)
        self.mol_viewer.setPixmap(pixmap)
        self.atom_coords = [drawer.GetDrawCoords(i) for i in range(smarts_mol.GetNumAtoms())]

    def _on_viewer_clicked(self, x: float, y: float):
        if not self.smarts_pattern:
            return
        if self.ambiguous_smarts:
            QMessageBox.warning(
                self,
                "Atomic Mapping Blocked",
                "SMARTS matching is ambiguous for this input set. Atomic mapping is disabled for this run.",
            )
            return
        if not self.atomic_mapping_enabled:
            QMessageBox.warning(
                self,
                "Atomic Mapping Blocked",
                self.atomic_mapping_block_reason or "Atomic descriptors are disabled for this SMARTS pattern.",
            )
            return

        atom_index = self._atom_at_position(x, y)
        if atom_index is None:
            return

        if atom_index in self.selected_atoms:
            self.selected_atoms.remove(atom_index)
        else:
            self.selected_atoms.append(atom_index)

        smarts_mol = Chem.MolFromSmarts(self.smarts_pattern)
        if smarts_mol is not None:
            self._render_smarts_molecule(smarts_mol)

        if self.selected_atoms:
            self.mol_info_label.setText(f"{len(self.selected_atoms)} atom(s) selected for atomic descriptors.")
        else:
            self.mol_info_label.setText("No atoms selected. Run will produce molecular descriptors only.")

    def _atom_at_position(self, x: float, y: float) -> int | None:
        if not self.atom_coords:
            return None

        pixmap = self.mol_viewer.pixmap()
        if pixmap is None or pixmap.isNull():
            return None

        x_offset = (self.mol_viewer.width() - pixmap.width()) / 2.0
        y_offset = (self.mol_viewer.height() - pixmap.height()) / 2.0

        local_x = x - x_offset
        local_y = y - y_offset
        if local_x < 0 or local_y < 0 or local_x > pixmap.width() or local_y > pixmap.height():
            return None

        scale_x = self._draw_width / max(1.0, float(pixmap.width()))
        scale_y = self._draw_height / max(1.0, float(pixmap.height()))
        draw_x = local_x * scale_x
        draw_y = local_y * scale_y

        threshold = 250
        for idx, coord in enumerate(self.atom_coords):
            if (coord.x - draw_x) ** 2 + (coord.y - draw_y) ** 2 < threshold:
                return idx
        return None

    def _validate(self) -> bool:
        if not self.file_paths:
            QMessageBox.warning(self, "Input Error", "Please add at least one file.")
            return False

        missing_files = [path for path in self.file_paths if not Path(path).exists()]
        if missing_files:
            QMessageBox.warning(
                self,
                "Input Error",
                "The following input file(s) do not exist:\n\n" + "\n".join(missing_files),
            )
            return False

        atoms_text = self.atoms_input.text().strip()
        if atoms_text and any(not part.strip() for part in atoms_text.split(",")):
            QMessageBox.warning(self, "Input Error", "QDESCP Atoms must be a comma-separated list of non-empty values.")
            return False

        return True

    def prepare_run_payload(self) -> dict:
        atoms_text = self.atoms_input.text().strip()
        qdescp_atoms = [part.strip() for part in atoms_text.split(",") if part.strip()] if atoms_text else []

        params = {
            "files": self.get_files(),
            "charge": self.charge_spin.value(),
            "mult": self.mult_spin.value(),
            "qdescp_atoms": qdescp_atoms,
            "boltz": self.boltzmann_checkbox.isChecked(),
            "robert": self.robert_checkbox.isChecked(),
        }

        if self.ambiguous_smarts:
            self._set_status("Ambiguous SMARTS: running molecular descriptors only.")
            return params

        if not self.atomic_mapping_enabled:
            if self.atomic_mapping_block_reason:
                self._set_status("Shared core too small: running molecular descriptors only.")
            return params

        if not self.selected_atoms:
            self._set_status("No SMARTS atom selection detected: running molecular descriptors only.")
            return params

        if self.csv_df is None or self.smiles_column is None or not self.smarts_pattern:
            self._set_status("Atomic mapping requires CSV + SMARTS. Running molecular descriptors only.")
            return params

        mapped = generate_mapped_smiles(
            self.smarts_pattern,
            self.selected_atoms,
            self.csv_df[self.smiles_column].dropna().astype(str).tolist(),
        )

        if mapped.ambiguous_smiles:
            self._set_status("Ambiguous SMARTS matches found during mapping: atomic descriptors skipped.")
            return params

        mapped_df = self.csv_df.copy()
        mapped_df[self.smiles_column] = mapped.mapped_smiles
        with tempfile.NamedTemporaryFile(prefix="qdescp_mapped_", suffix=".csv", delete=False) as tmp:
            mapped_df.to_csv(tmp.name, index=False)
            self._mapped_csv_temp_path = tmp.name
            params["files"] = [tmp.name]

        self._set_status("Mapped CSV generated from SMARTS atom selection.")
        return params

    def _handle_run(self):
        if self.run_callback is not None:
            self.run_callback()

    def set_output_paths(self, output_paths: list[str]):
        self.last_output_paths = output_paths
        self.output_list.clear()
        for path in output_paths:
            self.output_list.addItem(QListWidgetItem(path))
        self.open_output_folder_btn.setEnabled(bool(output_paths))
        if output_paths:
            self.outputs_toggle_btn.setText(f"Show Outputs ({len(output_paths)})")
        else:
            self.outputs_toggle_btn.setText("Show Outputs")

    def open_output_folder(self):
        if not self.last_output_paths:
            return
        folder = str(Path(self.last_output_paths[0]).resolve().parent)
        QDesktopServices.openUrl(QUrl.fromLocalFile(folder))


class NMRTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        label = QLabel("NMR tab is unchanged in this phase.")
        label.setWordWrap(True)
        layout.addWidget(label)
        layout.addStretch(1)


class QDESCP(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.main_parent = parent
        self.worker: QDESCPWorker | None = None
        self.csearch_window: CSEARCH_window | None = None
        self.run_start_time: float = 0.0
        self.run_watch_dirs: set[Path] = set()

        self.setStyleSheet(stylesheets.QWidget)
        self.setWindowTitle("QDESCP")
        self.setMinimumSize(980, 760)
        self._build_ui()

    def _build_ui(self):
        self.status_bar = QStatusBar(self)
        self.setStatusBar(self.status_bar)

        central = QWidget(self)
        layout = QVBoxLayout(central)

        self.tab_widget = QTabWidget(central)
        self.descriptors_tab = DescriptorsTab(
            status_callback=self.set_status,
            open_csearch_callback=self.open_csearch,
            run_callback=self.run_qdescp,
            parent=self,
        )
        self.nmr_tab = NMRTab(parent=self)
        self.tab_widget.addTab(self.descriptors_tab, "Descriptors")
        self.tab_widget.addTab(self.nmr_tab, "NMR")
        layout.addWidget(self.tab_widget)

        self.setCentralWidget(central)

    def set_status(self, message: str):
        self.status_bar.showMessage(message)

    def set_input_payload(self, payload: dict):
        files = payload.get("files", [])
        self.descriptors_tab.set_files(files)
        self.descriptors_tab.set_prefill(charge=payload.get("charge"), mult=payload.get("mult"))
        self.set_status("QDESCP descriptors payload updated.")

    def open_csearch(self):
        if self.csearch_window is None or not self.csearch_window.isVisible():
            self.csearch_window = CSEARCH_window(self._on_csearch_sdf_ready)
            self.csearch_window.show()
            return
        self.csearch_window.raise_()
        self.csearch_window.activateWindow()

    def _on_csearch_sdf_ready(self, sdf_paths: list[str]):
        if not sdf_paths:
            return
        self.tab_widget.setCurrentWidget(self.descriptors_tab)
        self.descriptors_tab.add_files(sdf_paths)
        self.set_status(f"Imported {len(sdf_paths)} SDF file(s) from CSEARCH.")

    def run_qdescp(self):
        if self.worker is not None and self.worker.isRunning():
            self.worker.request_stop()
            self.worker.quit()
            self.worker.wait()
            self.worker = None
            self.descriptors_tab.run_button.setText("Run QDESCP")
            self.set_status("QDESCP run cancelled.")
            return

        if not self.descriptors_tab._validate():
            return

        params = self.descriptors_tab.prepare_run_payload()
        self.run_start_time = time()
        self.run_watch_dirs = self._watch_directories(params.get("files", []))

        self.worker = QDESCPWorker(params)
        self.worker.success_signal.connect(self.success)
        self.worker.failure_signal.connect(self.failure)
        self.worker.finished_signal.connect(self._on_worker_finished)
        self.descriptors_tab.run_button.setText("Cancel")
        self.worker.start()
        self.set_status("QDESCP run started...")

    def _watch_directories(self, file_paths: Iterable[str]) -> set[Path]:
        dirs = set()
        for path in file_paths:
            try:
                dirs.add(Path(path).resolve().parent)
            except Exception:
                continue
        return dirs

    def _discover_outputs(self) -> list[str]:
        outputs: list[str] = []
        for folder in self.run_watch_dirs:
            if not folder.exists():
                continue
            for extension in ("*.csv", "*.json"):
                for path in folder.rglob(extension):
                    try:
                        if path.stat().st_mtime >= self.run_start_time:
                            outputs.append(str(path))
                    except OSError:
                        continue
        return sorted(set(outputs))

    def _on_worker_finished(self):
        self.descriptors_tab.run_button.setText("Run QDESCP")
        output_paths = self._discover_outputs()
        self.descriptors_tab.set_output_paths(output_paths)
        self.set_status("QDESCP run complete.")
        self.worker = None

    def success(self, message: str):
        msg = QMessageBox(self)
        msg.setWindowTitle("Success")
        msg.setText(message)
        pixmap = QPixmap(Icons.green)
        if not pixmap.isNull():
            icon = QIcon(pixmap)
            msg.setWindowIcon(icon)
            msg.setIconPixmap(icon.pixmap(64, 64))
        else:
            msg.setIcon(QMessageBox.Icon.Information)
        msg.exec()

    def failure(self, message: str):
        msg = QMessageBox(self)
        msg.setWindowTitle("Error")
        msg.setText(message)
        pixmap = QPixmap(Icons.red)
        if not pixmap.isNull():
            icon = QIcon(pixmap)
            msg.setWindowIcon(icon)
            msg.setIconPixmap(icon.pixmap(64, 64))
        else:
            msg.setIcon(QMessageBox.Icon.Critical)
        msg.exec()

    def closeEvent(self, event):
        self.descriptors_tab.shutdown_workers()

        if self.worker and self.worker.isRunning():
            self.worker.request_stop()
            self.worker.quit()
            self.worker.wait()

        if self.main_parent is not None and hasattr(self.main_parent, "button_for_qdescp"):
            self.main_parent.button_for_qdescp.setEnabled(True)

        super().closeEvent(event)


class CSEARCH_window(QWidget):
    def __init__(self, on_sdf_ready=None):
        super().__init__()
        self.on_sdf_ready = on_sdf_ready
        self.setWindowTitle("CSEARCH for QDESCP")
        self.setStyleSheet(stylesheets.QWidget)
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        self.csearch = CSEARCH()
        self.csearch.worker.finished_signal.connect(self._on_csearch_finished)
        layout.addWidget(self.csearch)

        close_button = QPushButton("Close")
        close_button.setStyleSheet(stylesheets.QPushButton)
        close_button.clicked.connect(self.close)
        layout.addWidget(close_button)

    def _on_csearch_finished(self, _message: str):
        if self.on_sdf_ready is None:
            return
        sdf_paths = self.csearch.get_generated_sdf_paths()
        if sdf_paths:
            self.on_sdf_ready(sdf_paths)
