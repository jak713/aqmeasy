"""Parameter controls for CMIN execution"""
import ast

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox,
    QComboBox, QLineEdit, QSpinBox, QDoubleSpinBox,
    QScrollArea, QFormLayout, QPushButton, QMessageBox
)

from aqmeasy.ui.CMIN_ui.CMIN_constraints_dialog import CMINConstraintDialog


class ParameterPanel(QWidget):
    """Widget for CMIN parameter configuration"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.init_ui()
        
    def init_ui(self):
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # Scroll area for all parameters
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll_content = QWidget()
        layout = QHBoxLayout(scroll_content)
        scroll.setWidget(scroll_content)
        main_layout.addWidget(scroll)
        
        # === GENERAL PARAMETERS ===
        general_group = QGroupBox("General Settings")
        general_group.setMaximumWidth(400)
        general_layout = QFormLayout()
        general_group.setLayout(general_layout)
        layout.addWidget(general_group)
        
        # Program selector
        self.program_combo = QComboBox()
        self.program_combo.addItems(['xtb', 'ani'])
        self.program_combo.setCurrentText('ani')
        self.program_combo.currentTextChanged.connect(self._on_program_changed)
        general_layout.addRow("Program:", self.program_combo)
        
        # Energy window
        self.ewin_spin = QDoubleSpinBox()
        self.ewin_spin.setRange(0.1, 100.0)
        self.ewin_spin.setValue(5.0)
        self.ewin_spin.setSuffix(" kcal/mol")
        general_layout.addRow("Energy Window:", self.ewin_spin)
        
        # Thresholds
        self.initial_e_threshold = QDoubleSpinBox()
        self.initial_e_threshold.setRange(0.0001, 10.0)
        self.initial_e_threshold.setValue(0.0001)
        self.initial_e_threshold.setDecimals(4)
        general_layout.addRow("Initial E Threshold:", self.initial_e_threshold)
        
        self.energy_threshold = QDoubleSpinBox()
        self.energy_threshold.setRange(0.01, 10.0)
        self.energy_threshold.setValue(0.25)
        general_layout.addRow("Energy Threshold:", self.energy_threshold)
        
        self.rms_threshold = QDoubleSpinBox()
        self.rms_threshold.setRange(0.01, 10.0)
        self.rms_threshold.setValue(0.25)
        general_layout.addRow("RMS Threshold:", self.rms_threshold)
        
        # Processors
        self.nprocs_spin = QSpinBox()
        self.nprocs_spin.setRange(1, 64)
        self.nprocs_spin.setValue(4)
        general_layout.addRow("Processors:", self.nprocs_spin)
        
        # Charge/Mult (optional overrides)
        self.charge_input = QLineEdit()
        self.charge_input.setPlaceholderText("Auto-detect")
        general_layout.addRow("Charge (optional):", self.charge_input)
        
        self.mult_input = QLineEdit()
        self.mult_input.setPlaceholderText("Auto-detect")
        general_layout.addRow("Multiplicity (optional):", self.mult_input)
        
        # Prefix/Suffix
        self.prefix_input = QLineEdit()
        general_layout.addRow("Prefix:", self.prefix_input)
        
        self.suffix_input = QLineEdit()
        general_layout.addRow("Suffix:", self.suffix_input)
        
        


        right_box = QVBoxLayout()


        # === xTB PARAMETERS ===
        self.xtb_group = QGroupBox("xTB Settings")
        self.xtb_group.setMinimumWidth(300)
        xtb_layout = QFormLayout()
        self.xtb_group.setLayout(xtb_layout)
        right_box.addWidget(self.xtb_group)
        
        self.xtb_keywords = QLineEdit()
        self.xtb_keywords.setPlaceholderText("e.g., --alpb ch2cl2 --gfn 1")
        xtb_layout.addRow("Additional Keywords:", self.xtb_keywords)
        
        self.constraints_atoms = QLineEdit()
        self.constraints_atoms.setPlaceholderText("e.g., 1,2,5")
        atoms_row = QHBoxLayout()
        atoms_row.addWidget(self.constraints_atoms)
        self.constraints_picker_button = QPushButton("Open Constraint Picker")
        self.constraints_picker_button.clicked.connect(self._open_constraints_picker)
        atoms_row.addWidget(self.constraints_picker_button)
        atoms_widget = QWidget()
        atoms_widget.setLayout(atoms_row)
        xtb_layout.addRow("Constrained Atoms:", atoms_widget)
        
        self.constraints_dist = QLineEdit()
        self.constraints_dist.setPlaceholderText("e.g., [[1,2,1.8],[4,5,2.0]]")
        xtb_layout.addRow("Distance Constraints:", self.constraints_dist)
        
        self.constraints_angle = QLineEdit()
        self.constraints_angle.setPlaceholderText("e.g., [[1,2,3,180],[4,5,6,120]]")
        xtb_layout.addRow("Angle Constraints:", self.constraints_angle)
        
        self.constraints_dihedral = QLineEdit()
        self.constraints_dihedral.setPlaceholderText("e.g., [[1,2,3,4,180]]")
        xtb_layout.addRow("Dihedral Constraints:", self.constraints_dihedral)
        
        self.stacksize = QSpinBox()
        self.stacksize.setRange(1, 10)
        self.stacksize.setValue(1)
        self.stacksize.setSuffix(" GB")
        xtb_layout.addRow("Stack Size:", self.stacksize)
        
        # === ANI PARAMETERS ===
        self.ani_group = QGroupBox("ANI Settings")
        self.ani_group.setMinimumWidth(300)
        ani_layout = QFormLayout()
        self.ani_group.setLayout(ani_layout)
        right_box.addWidget(self.ani_group)
        
        self.opt_steps = QSpinBox()
        self.opt_steps.setRange(100, 10000)
        self.opt_steps.setValue(1000)
        ani_layout.addRow("Max Opt Steps:", self.opt_steps)
        
        self.opt_fmax = QDoubleSpinBox()
        self.opt_fmax.setRange(0.001, 1.0)
        self.opt_fmax.setValue(0.05)
        self.opt_fmax.setDecimals(3)
        ani_layout.addRow("Opt Fmax:", self.opt_fmax)
        
        self.ani_method = QComboBox()
        self.ani_method.addItems(['ANI2x', 'ANI1x', 'ANI1ccx'])
        ani_layout.addRow("ANI Method:", self.ani_method)
        
        # Set initial visibility
        self._on_program_changed('ani')
        
        layout.addLayout(right_box)
        
    def _on_program_changed(self, program):
        """Show/hide parameter groups based on program selection"""
        if program == 'xtb':
            self.xtb_group.show()
            self.ani_group.hide()
        else:
            self.xtb_group.hide()
            self.ani_group.show()
    
    def get_parameters(self):
        """Collect all parameters as dictionary for CMIN"""
        params = {
            'program': self.program_combo.currentText(),
            'ewin_cmin': self.ewin_spin.value(),
            'initial_energy_threshold': self.initial_e_threshold.value(),
            'energy_threshold': self.energy_threshold.value(),
            'rms_threshold': self.rms_threshold.value(),
            'nprocs': self.nprocs_spin.value(),
        }
        
        # Optional charge/mult
        if self.charge_input.text().strip():
            try:
                params['charge'] = int(self.charge_input.text().strip())
            except ValueError:
                pass
        
        if self.mult_input.text().strip():
            try:
                params['mult'] = int(self.mult_input.text().strip())
            except ValueError:
                pass
        
        # Prefix/suffix
        if self.prefix_input.text():
            params['prefix'] = self.prefix_input.text()
        if self.suffix_input.text():
            params['suffix'] = self.suffix_input.text()
        
        # xTB parameters
        if self.program_combo.currentText() == 'xtb':
            if self.xtb_keywords.text():
                params['xtb_keywords'] = self.xtb_keywords.text()
            if self.constraints_atoms.text():
                parsed_atoms = self._parse_constraints_atoms(self.constraints_atoms.text())
                if parsed_atoms is None:
                    raise ValueError("Invalid constrained atoms. Use comma-separated integers, e.g. 1,2,5")
                if parsed_atoms:
                    params['constraints_atoms'] = parsed_atoms
            if self.constraints_dist.text():
                parsed_dist = self._parse_constraint_matrix(self.constraints_dist.text(), expected_len=3)
                if parsed_dist is None:
                    raise ValueError("Invalid distance constraints. Use [[a,b,value], ...]")
                if parsed_dist:
                    params['constraints_dist'] = parsed_dist
            if self.constraints_angle.text():
                parsed_angle = self._parse_constraint_matrix(self.constraints_angle.text(), expected_len=4)
                if parsed_angle is None:
                    raise ValueError("Invalid angle constraints. Use [[a,b,c,value], ...]")
                if parsed_angle:
                    params['constraints_angle'] = parsed_angle
            if self.constraints_dihedral.text():
                parsed_dihedral = self._parse_constraint_matrix(self.constraints_dihedral.text(), expected_len=5)
                if parsed_dihedral is None:
                    raise ValueError("Invalid dihedral constraints. Use [[a,b,c,d,value], ...]")
                if parsed_dihedral:
                    params['constraints_dihedral'] = parsed_dihedral
            params['stacksize'] = self.stacksize.text()
        
        # ANI parameters
        if self.program_combo.currentText() == 'ani':
            params['opt_steps'] = self.opt_steps.value()
            params['opt_fmax'] = self.opt_fmax.value()
            params['ani_method'] = self.ani_method.currentText()
        
        return params

    def _open_constraints_picker(self):
        """Open popup to graphically select xTB constraints by atom indices."""
        parent_window = self.parent()
        files = []
        selected_file = None
        if parent_window is not None and hasattr(parent_window, 'file_panel'):
            model = getattr(parent_window.file_panel, 'model', None)
            if model is not None:
                files = list(getattr(model, 'files', []) or [])
            selected_file = getattr(model, 'currently_selected_file', None)

        if not files:
            QMessageBox.information(self, "No input files", "Add at least one SDF/XYZ file before opening the constraint picker.")
            return

        dialog = CMINConstraintDialog(
            files=files,
            initial_constraints={
                'constraints_atoms': self.constraints_atoms.text(),
                'constraints_dist': self.constraints_dist.text(),
                'constraints_angle': self.constraints_angle.text(),
                'constraints_dihedral': self.constraints_dihedral.text(),
            },
            selected_file=selected_file,
            parent=self,
        )
        if dialog.exec() != dialog.DialogCode.Accepted:
            return

        payload = dialog.export_constraints()
        self.constraints_atoms.setText(payload.get('constraints_atoms', ''))
        self.constraints_dist.setText(payload.get('constraints_dist', ''))
        self.constraints_angle.setText(payload.get('constraints_angle', ''))
        self.constraints_dihedral.setText(payload.get('constraints_dihedral', ''))

    def _parse_constraints_atoms(self, text):
        raw = str(text or '').strip()
        if not raw:
            return []

        if raw.startswith('['):
            try:
                parsed = ast.literal_eval(raw)
            except Exception:
                return None
            if not isinstance(parsed, list):
                return None
            values = parsed
        else:
            values = [token.strip() for token in raw.split(',') if token.strip()]

        converted = []
        for value in values:
            try:
                atom_id = int(value)
            except Exception:
                return None
            if atom_id <= 0:
                return None
            converted.append(atom_id)
        return sorted(set(converted))

    def _parse_constraint_matrix(self, text, expected_len):
        raw = str(text or '').strip()
        if not raw:
            return []

        try:
            parsed = ast.literal_eval(raw)
        except Exception:
            return None

        if not isinstance(parsed, list):
            return None

        converted = []
        for row in parsed:
            if not isinstance(row, (list, tuple)) or len(row) != expected_len:
                return None

            new_row = []
            for idx, value in enumerate(row):
                try:
                    if idx < expected_len - 1:
                        atom_id = int(value)
                        if atom_id <= 0:
                            return None
                        new_row.append(atom_id)
                    else:
                        new_row.append(float(value))
                except Exception:
                    return None
            converted.append(new_row)

        return converted