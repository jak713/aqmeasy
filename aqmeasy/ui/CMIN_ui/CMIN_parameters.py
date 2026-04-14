"""Parameter controls for CMIN execution"""
import ast

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox,
    QComboBox, QLineEdit, QSpinBox, QDoubleSpinBox,
    QPushButton, QMessageBox,
    QGridLayout, QLabel, QSizePolicy
)
from PySide6.QtGui import QIcon

from aqmeasy.ui.CMIN_ui.CMIN_constraints_dialog import CMINConstraintDialog
from aqmeasy.ui.icons import Icons


class ParameterPanel(QWidget):
    """Widget for CMIN parameter configuration"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.init_ui()
        
    def init_ui(self):
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)

        # Always-visible core controls
        core_group = QGroupBox("Core Settings")
        core_layout = QGridLayout()
        core_layout.setHorizontalSpacing(12)
        core_layout.setVerticalSpacing(8)
        core_group.setLayout(core_layout)
        main_layout.addWidget(core_group)
        self.program_combo = QComboBox()
        self.program_combo.addItems(['xtb', 'ani'])
        self.program_combo.setCurrentText('ani')
        self.program_combo.currentTextChanged.connect(self._on_program_changed)

        self.ewin_spin = QDoubleSpinBox()
        self.ewin_spin.setRange(0.1, 100.0)
        self.ewin_spin.setValue(5.0)
        self.ewin_spin.setSuffix(" kcal/mol")

        self.initial_e_threshold = QDoubleSpinBox()
        self.initial_e_threshold.setRange(0.0001, 10.0)
        self.initial_e_threshold.setValue(0.0001)
        self.initial_e_threshold.setDecimals(4)
        self.initial_e_threshold.setSingleStep(0.01)

        self.energy_threshold = QDoubleSpinBox()
        self.energy_threshold.setRange(0.01, 10.0)
        self.energy_threshold.setValue(0.25)
        self.energy_threshold.setSingleStep(0.01)

        self.rms_threshold = QDoubleSpinBox()
        self.rms_threshold.setRange(0.01, 10.0)
        self.rms_threshold.setValue(0.25)
        self.rms_threshold.setSingleStep(0.01)
        
        core_labels = [
            ("Program:", self.program_combo),
            ("Energy Window:", self.ewin_spin),
            ("Initial E Threshold:", self.initial_e_threshold),
            ("Energy Threshold:", self.energy_threshold),
            ("RMS Threshold:", self.rms_threshold),
        ]
        for index, (label_text, widget) in enumerate(core_labels):
            row = index // 3
            column = (index % 3) * 2
            label = QLabel(label_text)
            label.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Preferred)
            widget.setMinimumWidth(140)
            core_layout.addWidget(label, row, column)
            core_layout.addWidget(widget, row, column + 1)

        self.advanced_toggle_button = QPushButton("Advanced Settings")
        self.advanced_toggle_button.setCheckable(True)
        self.advanced_toggle_button.setIcon(QIcon(Icons.eye))
        self.advanced_toggle_button.toggled.connect(self._toggle_advanced_panel)
        core_layout.addWidget(self.advanced_toggle_button, 1, 4, 1, 2)

        # Collapsible advanced settings in one compact group
        self.advanced_group = QGroupBox("Advanced Settings")
        self.advanced_layout = QGridLayout()
        self.advanced_layout.setHorizontalSpacing(12)
        self.advanced_layout.setVerticalSpacing(8)
        self.advanced_group.setLayout(self.advanced_layout)
        self.advanced_group.setVisible(False)
        main_layout.addWidget(self.advanced_group)

        self.nprocs_spin = QSpinBox()
        self.nprocs_spin.setRange(1, 64)
        self.nprocs_spin.setValue(4)
        self.charge_input = QLineEdit()
        self.charge_input.setPlaceholderText("Auto-detect")
        self.mult_input = QLineEdit()
        self.mult_input.setPlaceholderText("Auto-detect")
        self.prefix_input = QLineEdit()
        self.suffix_input = QLineEdit()

        self.general_container = QWidget()
        general_layout = QGridLayout(self.general_container)
        general_layout.setHorizontalSpacing(10)
        general_layout.setVerticalSpacing(8)
        general_title = QLabel("General Advanced Settings")
        general_title.setStyleSheet("font-weight: bold;")
        general_layout.addWidget(general_title, 0, 0, 1, 4)

        general_fields = [
            ("Processors:", self.nprocs_spin),
            ("Charge (optional):", self.charge_input),
            ("Multiplicity (optional):", self.mult_input),
            ("Prefix:", self.prefix_input),
            ("Suffix:", self.suffix_input),
        ]

        for index, (label_text, widget) in enumerate(general_fields):
            row = 1 + (index // 2)
            column = (index % 2) * 2
            general_layout.addWidget(QLabel(label_text), row, column)
            general_layout.addWidget(widget, row, column + 1)

        self.program_container = QWidget()
        program_layout = QVBoxLayout(self.program_container)
        program_layout.setContentsMargins(0, 0, 0, 0)
        program_layout.setSpacing(8)

        self.xtb_container = QWidget()
        xtb_layout = QGridLayout(self.xtb_container)
        xtb_layout.setHorizontalSpacing(10)
        xtb_layout.setVerticalSpacing(8)
        xtb_layout.setColumnStretch(1, 1)
        xtb_layout.setColumnStretch(3, 1)
        xtb_title = QLabel("xTB Settings")
        xtb_title.setStyleSheet("font-weight: bold;")
        xtb_layout.addWidget(xtb_title, 0, 0, 1, 4)

        self.xtb_keywords = QLineEdit()
        self.xtb_keywords.setPlaceholderText("e.g., --alpb ch2cl2 --gfn 1")
        
        keywords_stack_row = QWidget()
        keywords_stack_layout = QHBoxLayout(keywords_stack_row)
        keywords_stack_layout.setContentsMargins(0, 0, 0, 0)
        keywords_stack_layout.setSpacing(8)

        self.constraints_atoms = QLineEdit()
        self.constraints_atoms.setPlaceholderText("e.g., 1,2,5")
        atoms_row = QHBoxLayout()
        atoms_row.addWidget(self.constraints_atoms)
        self.constraints_picker_button = QPushButton("Open Constraint Picker")
        self.constraints_picker_button.clicked.connect(self._open_constraints_picker)
        atoms_row.addWidget(self.constraints_picker_button)
        atoms_widget = QWidget()
        atoms_widget.setLayout(atoms_row)
        xtb_layout.addWidget(QLabel("Constrained Atoms:"), 2, 0)
        xtb_layout.addWidget(atoms_widget, 2, 1, 1, 2)

        self.stacksize = QSpinBox()
        self.stacksize.setRange(1, 10)
        self.stacksize.setValue(1)
        self.stacksize.setSuffix(" GB")

        keywords_stack_layout.addWidget(self.xtb_keywords, 1)
        keywords_stack_layout.addWidget(QLabel("Stack Size:"))
        keywords_stack_layout.addWidget(self.stacksize)

        xtb_layout.addWidget(QLabel("Additional Keywords:"), 1, 0)
        xtb_layout.addWidget(keywords_stack_row, 1, 1, 1, 3)

        self.constraints_dist = QLineEdit()
        self.constraints_dist.setPlaceholderText("e.g., [[1,2,1.8],[4,5,2.0]]")
        xtb_layout.addWidget(QLabel("Distance Constraints:"), 3, 0)
        xtb_layout.addWidget(self.constraints_dist, 3, 1)
        
        self.constraints_angle = QLineEdit()
        self.constraints_angle.setPlaceholderText("e.g., [[1,2,3,180],[4,5,6,120]]")
        xtb_layout.addWidget(QLabel("Angle Constraints:"), 3, 2)
        xtb_layout.addWidget(self.constraints_angle, 3, 3)
        
        self.constraints_dihedral = QLineEdit()
        self.constraints_dihedral.setPlaceholderText("e.g., [[1,2,3,4,180]]")
        xtb_layout.addWidget(QLabel("Dihedral Constraints:"), 4, 0)
        xtb_layout.addWidget(self.constraints_dihedral, 4, 1, 1, 3)

        self.ani_container = QWidget()
        ani_layout = QGridLayout(self.ani_container)
        ani_layout.setHorizontalSpacing(10)
        ani_layout.setVerticalSpacing(8)
        ani_title = QLabel("ANI Settings")
        ani_title.setStyleSheet("font-weight: bold;")
        ani_layout.addWidget(ani_title, 0, 0, 1, 2)

        self.opt_steps = QSpinBox()
        self.opt_steps.setRange(100, 10000)
        self.opt_steps.setValue(1000)
        ani_layout.addWidget(QLabel("Max Opt Steps:"), 1, 0)
        ani_layout.addWidget(self.opt_steps, 1, 1)
        
        self.opt_fmax = QDoubleSpinBox()
        self.opt_fmax.setRange(0.001, 1.0)
        self.opt_fmax.setValue(0.05)
        self.opt_fmax.setDecimals(3)
        ani_layout.addWidget(QLabel("Opt Fmax:"), 2, 0)
        ani_layout.addWidget(self.opt_fmax, 2, 1)
        
        self.ani_method = QComboBox()
        self.ani_method.addItems(['ANI2x', 'ANI1x', 'ANI1ccx'])
        ani_layout.addWidget(QLabel("ANI Method:"), 3, 0)
        ani_layout.addWidget(self.ani_method, 3, 1)

        program_layout.addWidget(self.xtb_container)
        program_layout.addWidget(self.ani_container)

        self.advanced_layout.addWidget(self.general_container, 0, 0)
        self.advanced_layout.addWidget(self.program_container, 0, 1)
        self.advanced_layout.setColumnStretch(0, 1)
        self.advanced_layout.setColumnStretch(1, 1)

        main_layout.addStretch()

        # Set initial visibility
        self._on_program_changed('ani')

    def _toggle_advanced_panel(self, checked):
        """Expand/collapse advanced CMIN options."""
        self.advanced_group.setVisible(checked)
        self.advanced_toggle_button.setIcon(QIcon(Icons.eye_crossed if checked else Icons.eye))
        self._on_program_changed(self.program_combo.currentText())
        
    def _on_program_changed(self, program):
        """Show/hide parameter groups based on program selection"""
        xtb_visible = program == 'xtb'
        ani_visible = program == 'ani'
        self.xtb_container.setVisible(self.advanced_group.isVisible() and xtb_visible)
        self.ani_container.setVisible(self.advanced_group.isVisible() and ani_visible)
    
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
        file_panel = getattr(parent_window, 'file_panel', None) if parent_window is not None else None
        if file_panel is not None:
            model = getattr(file_panel, 'model', None)
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