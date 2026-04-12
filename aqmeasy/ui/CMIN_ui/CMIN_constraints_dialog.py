"""Constraint picker dialog for CMIN xTB runs."""

import os
import tempfile
from typing import List

from PySide6.QtCore import QPoint, Qt, Signal
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)


class MoleculeCanvas(QLabel):
    """Clickable label used to pick atoms from 2D molecule depiction."""

    atomClicked = Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.atom_coords = []
        self.source_width = 1.0
        self.source_height = 1.0
        self.display_width = 1.0
        self.display_height = 1.0
        self.offset_x = 0.0
        self.offset_y = 0.0
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setMinimumSize(520, 380)
        self.setStyleSheet("border: 1px solid #d5d8dc; background: #ffffff;")

    def set_atom_coords(self, coords):
        self.atom_coords = coords or []

    def set_render_map(self, source_width, source_height, display_width, display_height):
        self.source_width = max(1.0, float(source_width))
        self.source_height = max(1.0, float(source_height))
        self.display_width = max(1.0, float(display_width))
        self.display_height = max(1.0, float(display_height))
        self.offset_x = max(0.0, (float(self.width()) - self.display_width) / 2.0)
        self.offset_y = max(0.0, (float(self.height()) - self.display_height) / 2.0)

    def mousePressEvent(self, event):
        if not self.atom_coords:
            return
        x = float(event.position().x())
        y = float(event.position().y())
        hit = self._find_atom(x, y)
        if hit is not None:
            self.atomClicked.emit(hit)

    def _find_atom(self, x, y):
        # Map from QLabel coordinates to RDKit drawing coordinates.
        if x < self.offset_x or y < self.offset_y:
            return None
        if x > self.offset_x + self.display_width or y > self.offset_y + self.display_height:
            return None

        x = ((x - self.offset_x) / self.display_width) * self.source_width
        y = ((y - self.offset_y) / self.display_height) * self.source_height

        click_radius_sq = 15.0 * 15.0
        for idx, point in enumerate(self.atom_coords):
            dx = float(point.x) - x
            dy = float(point.y) - y
            if (dx * dx + dy * dy) <= click_radius_sq:
                return idx + 1
        return None


class CMINConstraintDialog(QDialog):
    """Popup dialog to compose xTB constraints using numbered atoms."""

    def __init__(self, files, initial_constraints=None, selected_file=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("xTB Constraint Picker")
        self.setMinimumSize(1080, 700)

        self.files = list(files or [])
        self.selected_file = os.path.abspath(selected_file) if selected_file else None
        self.initial_constraints = initial_constraints or {}

        self.records = []
        self.active_mol = None
        self.selected_atoms = []

        self.constraints_atoms = self._parse_int_list(self.initial_constraints.get("constraints_atoms", ""))
        self.constraints_dist = self._parse_nested_list(self.initial_constraints.get("constraints_dist", ""), expected_len=3)
        self.constraints_angle = self._parse_nested_list(self.initial_constraints.get("constraints_angle", ""), expected_len=4)
        self.constraints_dihedral = self._parse_nested_list(self.initial_constraints.get("constraints_dihedral", ""), expected_len=5)

        self._build_ui()
        self._load_records()

    def _build_ui(self):
        root = QVBoxLayout(self)

        top_info = QLabel(
            "Atom numbering is 1-based. For XYZ, numbering follows file atom line order. "
            "One representative structure per file is shown (constraints are global for the run)."
        )
        top_info.setWordWrap(True)
        top_info.setStyleSheet("color: #566573;")
        root.addWidget(top_info)

        control_row = QHBoxLayout()
        control_row.addWidget(QLabel("File:"))
        self.structure_combo = QComboBox()
        self.structure_combo.currentIndexChanged.connect(self._on_structure_changed)
        control_row.addWidget(self.structure_combo, 1)
        root.addLayout(control_row)

        content = QHBoxLayout()
        root.addLayout(content, 1)

        left_panel = QVBoxLayout()
        self.canvas = MoleculeCanvas()
        self.canvas.atomClicked.connect(self._on_atom_clicked)
        left_panel.addWidget(self.canvas, 1)

        self.selection_label = QLabel("Selected atoms: none")
        self.selection_label.setWordWrap(True)
        left_panel.addWidget(self.selection_label)

        self.clear_selection_button = QPushButton("Clear Selection")
        self.clear_selection_button.clicked.connect(self._clear_selection)
        left_panel.addWidget(self.clear_selection_button)

        content.addLayout(left_panel, 3)

        right_column = QVBoxLayout()

        build_group = QGroupBox("Build Constraint")
        build_layout = QFormLayout()
        build_group.setLayout(build_layout)

        self.mode_combo = QComboBox()
        self.mode_combo.addItems([
            "Constrained Atoms",
            "Distance",
            "Angle",
            "Dihedral",
        ])
        self.mode_combo.currentTextChanged.connect(self._on_mode_changed)
        build_layout.addRow("Type:", self.mode_combo)

        self.value_spin = QDoubleSpinBox()
        self.value_spin.setDecimals(4)
        self.value_spin.setRange(0.0, 360.0)
        self.value_spin.setSingleStep(0.1)
        build_layout.addRow("Target:", self.value_spin)

        self.add_button = QPushButton("Add Constraint From Selection")
        self.add_button.clicked.connect(self._add_constraint_from_selection)
        build_layout.addRow(self.add_button)

        right_column.addWidget(build_group)

        self.atoms_list = QListWidget()
        self.dist_list = QListWidget()
        self.angle_list = QListWidget()
        self.dihedral_list = QListWidget()

        right_column.addWidget(self._list_group("Constrained Atoms", self.atoms_list, self._remove_atom_entry))
        right_column.addWidget(self._list_group("Distance Constraints", self.dist_list, self._remove_dist_entry))
        right_column.addWidget(self._list_group("Angle Constraints", self.angle_list, self._remove_angle_entry))
        right_column.addWidget(self._list_group("Dihedral Constraints", self.dihedral_list, self._remove_dihedral_entry))

        content.addLayout(right_column, 2)

        footer = QHBoxLayout()
        footer.addStretch(1)
        clear_all_button = QPushButton("Clear All")
        clear_all_button.clicked.connect(self._clear_all)
        footer.addWidget(clear_all_button)

        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        footer.addWidget(cancel)

        apply_button = QPushButton("Apply")
        apply_button.clicked.connect(self.accept)
        footer.addWidget(apply_button)
        root.addLayout(footer)

        self._on_mode_changed(self.mode_combo.currentText())
        self._refresh_constraint_lists()

    def _list_group(self, title, list_widget, remove_callback):
        group = QGroupBox(title)
        layout = QVBoxLayout()
        group.setLayout(layout)
        list_widget.setMinimumHeight(90)
        layout.addWidget(list_widget)
        remove_button = QPushButton("Remove Selected")
        remove_button.clicked.connect(remove_callback)
        layout.addWidget(remove_button)
        return group

    def _load_records(self):
        try:
            from rdkit import Chem
            from rdkit.Chem import rdDetermineBonds
        except Exception:
            QMessageBox.warning(self, "Missing dependency", "RDKit is required for the constraint picker.")
            return

        records = []
        unreadable_files = []
        for file_path in self.files:
            abs_path = os.path.abspath(file_path)
            ext = os.path.splitext(file_path)[1].lower()
            if ext == ".sdf":
                try:
                    supplier = Chem.SDMolSupplier(file_path, removeHs=False, sanitize=False)
                    picked = None
                    for mol in supplier:
                        if mol is None:
                            continue
                        picked = mol
                        break
                    if picked is not None:
                        records.append({
                            "label": os.path.basename(file_path),
                            "file_path": abs_path,
                            "mol": picked,
                        })
                    else:
                        unreadable_files.append(os.path.basename(file_path))
                except Exception:
                    unreadable_files.append(os.path.basename(file_path))
                    continue
            elif ext == ".xyz":
                xyz_blocks = self._read_xyz_blocks(file_path)
                if not xyz_blocks:
                    unreadable_files.append(os.path.basename(file_path))
                    continue
                xyz_block = xyz_blocks[0]
                try:
                    mol = Chem.MolFromXYZBlock(xyz_block)
                    if mol is None:
                        unreadable_files.append(os.path.basename(file_path))
                        continue
                    try:
                        rdDetermineBonds.DetermineBonds(mol)
                    except Exception:
                        pass
                    records.append({
                        "label": os.path.basename(file_path),
                        "file_path": abs_path,
                        "mol": mol,
                    })
                except Exception:
                    unreadable_files.append(os.path.basename(file_path))
                    continue

        self.records = records
        self.structure_combo.blockSignals(True)
        self.structure_combo.clear()
        for record in self.records:
            self.structure_combo.addItem(record["label"])
        self.structure_combo.blockSignals(False)

        if self.records:
            preferred_index = 0
            if self.selected_file:
                for idx, record in enumerate(self.records):
                    if record.get("file_path") == self.selected_file:
                        preferred_index = idx
                        break
            self.structure_combo.setCurrentIndex(preferred_index)
            self._on_structure_changed(preferred_index)
        else:
            self.canvas.setText("No valid SDF/XYZ structures found in selected files.")

        if unreadable_files:
            self.selection_label.setText(
                "Could not render some files: " + ", ".join(unreadable_files[:6])
                + ("..." if len(unreadable_files) > 6 else "")
            )

    def _read_xyz_blocks(self, file_path):
        blocks: List[str] = []
        try:
            with open(file_path, "r", encoding="utf-8", errors="ignore") as handle:
                lines = handle.readlines()
        except Exception:
            return blocks

        index = 0
        total = len(lines)
        while index < total:
            line = lines[index].strip()
            if not line:
                index += 1
                continue
            try:
                atom_count = int(line)
            except ValueError:
                index += 1
                continue
            block_end = index + atom_count + 2
            if atom_count > 0 and block_end <= total:
                blocks.append("".join(lines[index:block_end]))
                index = block_end
            else:
                index += 1
        return blocks

    def _on_structure_changed(self, index):
        if index < 0 or index >= len(self.records):
            return
        self.active_mol = self.records[index]["mol"]
        self._clear_selection()
        self._render_active_molecule()

    def _on_atom_clicked(self, atom_id):
        if atom_id in self.selected_atoms:
            self.selected_atoms.remove(atom_id)
        else:
            self.selected_atoms.append(atom_id)
        self.selected_atoms = sorted(set(self.selected_atoms))
        self._update_selection_label()
        self._render_active_molecule()

    def _update_selection_label(self):
        if not self.selected_atoms:
            self.selection_label.setText("Selected atoms: none")
        else:
            self.selection_label.setText("Selected atoms: " + ", ".join(str(v) for v in self.selected_atoms))

    def _render_active_molecule(self):
        if self.active_mol is None:
            return

        try:
            from rdkit import Chem
            from rdkit.Chem import rdDepictor
            from rdkit.Chem.Draw import rdMolDraw2D
        except Exception:
            return

        primary_error = None
        for fallback_round in (0, 1):
            try:
                mol = Chem.Mol(self.active_mol)
                if fallback_round == 1:
                    # A lighter molecule often draws when strict valence/aromatic information is problematic.
                    mol = Chem.RemoveHs(mol, sanitize=False)

                if mol.GetNumAtoms() == 0:
                    raise ValueError("Molecule has no atoms.")

                if mol.GetNumConformers() == 0:
                    rdDepictor.Compute2DCoords(mol)

                width = max(520, self.canvas.width())
                height = max(380, self.canvas.height())
                drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
                draw_options = drawer.drawOptions()
                draw_options.addAtomIndices = False
                draw_options.prepareMolsBeforeDrawing = False
                for atom_idx in range(mol.GetNumAtoms()):
                    draw_options.atomLabels[atom_idx] = str(atom_idx + 1)

                highlight = [atom - 1 for atom in self.selected_atoms if atom > 0]
                colors = {atom_idx: (0.95, 0.76, 0.06) for atom_idx in highlight}
                if highlight:
                    try:
                        drawer.DrawMolecule(mol, highlightAtoms=highlight, highlightAtomColors=colors)
                    except Exception:
                        drawer.DrawMolecule(mol)
                else:
                    drawer.DrawMolecule(mol)

                atom_coords = [drawer.GetDrawCoords(i) for i in range(mol.GetNumAtoms())]
                drawer.FinishDrawing()

                png_path = os.path.join(tempfile.gettempdir(), "aqmeasy_cmin_constraints_preview.png")
                drawer.WriteDrawingText(png_path)
                pixmap = QPixmap(png_path)
                if pixmap.isNull():
                    raise ValueError("Rendered image was empty.")

                self.canvas.setPixmap(
                    pixmap.scaled(
                        self.canvas.size(),
                        Qt.AspectRatioMode.KeepAspectRatio,
                        Qt.TransformationMode.SmoothTransformation,
                    )
                )
                shown_pixmap = self.canvas.pixmap()
                if shown_pixmap is not None:
                    self.canvas.set_render_map(
                        source_width=width,
                        source_height=height,
                        display_width=shown_pixmap.width(),
                        display_height=shown_pixmap.height(),
                    )
                self.canvas.set_atom_coords(atom_coords)
                return
            except Exception as exc:
                if primary_error is None:
                    primary_error = exc

        self.canvas.setText(f"Unable to render molecule: {primary_error}")
        self.canvas.set_render_map(1.0, 1.0, 1.0, 1.0)
        self.canvas.set_atom_coords([])

    def _on_mode_changed(self, mode):
        if mode == "Constrained Atoms":
            self.value_spin.setEnabled(False)
            return

        self.value_spin.setEnabled(True)
        if mode == "Distance":
            self.value_spin.setRange(0.0, 20.0)
            self.value_spin.setValue(1.8)
            self.value_spin.setSuffix(" A")
        elif mode == "Angle":
            self.value_spin.setRange(0.0, 360.0)
            self.value_spin.setValue(120.0)
            self.value_spin.setSuffix(" deg")
        else:
            self.value_spin.setRange(-180.0, 180.0)
            self.value_spin.setValue(180.0)
            self.value_spin.setSuffix(" deg")

    def _add_constraint_from_selection(self):
        mode = self.mode_combo.currentText()
        selected = list(self.selected_atoms)

        if mode == "Constrained Atoms":
            if not selected:
                QMessageBox.information(self, "Select atoms", "Pick at least one atom first.")
                return
            merged = sorted(set(self.constraints_atoms + selected))
            self.constraints_atoms = merged
            self._clear_selection()
            self._refresh_constraint_lists()
            return

        required = {"Distance": 2, "Angle": 3, "Dihedral": 4}[mode]
        if len(selected) != required:
            QMessageBox.information(self, "Select atoms", f"Select exactly {required} atoms for {mode.lower()} constraint.")
            return

        value = float(self.value_spin.value())
        if mode == "Distance":
            if value <= 0:
                QMessageBox.warning(self, "Invalid value", "Distance must be greater than 0.")
                return
            self.constraints_dist.append([selected[0], selected[1], value])
        elif mode == "Angle":
            if value <= 0 or value >= 360:
                QMessageBox.warning(self, "Invalid value", "Angle must be between 0 and 360.")
                return
            self.constraints_angle.append([selected[0], selected[1], selected[2], value])
        else:
            if value < -180 or value > 180:
                QMessageBox.warning(self, "Invalid value", "Dihedral must be between -180 and 180.")
                return
            self.constraints_dihedral.append([selected[0], selected[1], selected[2], selected[3], value])

        self._clear_selection()
        self._refresh_constraint_lists()

    def _remove_atom_entry(self):
        row = self.atoms_list.currentRow()
        if row >= 0:
            self.constraints_atoms.pop(row)
            self._refresh_constraint_lists()

    def _remove_dist_entry(self):
        row = self.dist_list.currentRow()
        if row >= 0:
            self.constraints_dist.pop(row)
            self._refresh_constraint_lists()

    def _remove_angle_entry(self):
        row = self.angle_list.currentRow()
        if row >= 0:
            self.constraints_angle.pop(row)
            self._refresh_constraint_lists()

    def _remove_dihedral_entry(self):
        row = self.dihedral_list.currentRow()
        if row >= 0:
            self.constraints_dihedral.pop(row)
            self._refresh_constraint_lists()

    def _clear_selection(self):
        self.selected_atoms = []
        self._update_selection_label()
        self._render_active_molecule()

    def _clear_all(self):
        self.constraints_atoms = []
        self.constraints_dist = []
        self.constraints_angle = []
        self.constraints_dihedral = []
        self._refresh_constraint_lists()

    def _refresh_constraint_lists(self):
        self.atoms_list.clear()
        for atom in self.constraints_atoms:
            self.atoms_list.addItem(str(int(atom)))

        self.dist_list.clear()
        for item in self.constraints_dist:
            self.dist_list.addItem(str(item))

        self.angle_list.clear()
        for item in self.constraints_angle:
            self.angle_list.addItem(str(item))

        self.dihedral_list.clear()
        for item in self.constraints_dihedral:
            self.dihedral_list.addItem(str(item))

    def export_constraints(self):
        """Return constraints ready to be assigned to CMIN xTB text fields."""
        atoms_text = ",".join(str(int(x)) for x in sorted(set(self.constraints_atoms)))
        return {
            "constraints_atoms": atoms_text,
            "constraints_dist": str(self.constraints_dist) if self.constraints_dist else "",
            "constraints_angle": str(self.constraints_angle) if self.constraints_angle else "",
            "constraints_dihedral": str(self.constraints_dihedral) if self.constraints_dihedral else "",
        }

    def _parse_int_list(self, text):
        import ast

        raw = str(text or "").strip()
        if not raw:
            return []

        if raw.startswith("["):
            try:
                values = ast.literal_eval(raw)
            except Exception:
                return []
        else:
            values = [v.strip() for v in raw.split(",") if v.strip()]

        parsed = []
        for value in values:
            try:
                parsed.append(int(value))
            except Exception:
                continue
        return sorted(set(v for v in parsed if v > 0))

    def _parse_nested_list(self, text, expected_len):
        import ast

        raw = str(text or "").strip()
        if not raw:
            return []
        try:
            values = ast.literal_eval(raw)
        except Exception:
            return []

        parsed = []
        if not isinstance(values, list):
            return parsed

        for entry in values:
            if not isinstance(entry, (list, tuple)) or len(entry) != expected_len:
                continue
            numeric = []
            ok = True
            for idx, item in enumerate(entry):
                try:
                    if idx < expected_len - 1:
                        numeric.append(int(item))
                    else:
                        numeric.append(float(item))
                except Exception:
                    ok = False
                    break
            if ok:
                parsed.append(numeric)
        return parsed
