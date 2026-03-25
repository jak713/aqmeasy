# aqmeasy/ui/QCORR_ui/QCORR_analysisswidget.py

from typing import Dict, Any, Optional, List
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTabWidget,
    QTableWidget, QTableWidgetItem, QLabel, QPushButton,
    QComboBox, QSlider, QGroupBox, QHeaderView, QSplitter,
    QDoubleSpinBox
)
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtCore import Qt, Signal
import json
import math


class QCORRAnalysisWidget(QWidget):
    """
    Widget for displaying and analyzing QCORR JSON output data.
    
    Features:
    - Structured metadata display across multiple tabs
    - Interactive vibrational mode visualization using py3Dmol
    """
    
    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.json_data: Optional[Dict[str, Any]] = None
        self.init_ui()
    
    def init_ui(self):
        """Initialize the user interface."""
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # Main tab widget for different data categories
        self.tab_widget = QTabWidget()
        layout.addWidget(self.tab_widget)
        
        # Create individual tabs
        self.general_tab = GeneralInfoTab()
        self.system_tab = SystemInfoTab()
        self.energy_tab = EnergyInfoTab()
        self.vibration_tab = VibrationTab()
        self.orbital_tab = OrbitalInfoTab()
        
        self.tab_widget.addTab(self.general_tab, "General")
        self.tab_widget.addTab(self.system_tab, "System")
        self.tab_widget.addTab(self.energy_tab, "Energies")
        self.tab_widget.addTab(self.vibration_tab, "Vibrations")
        self.tab_widget.addTab(self.orbital_tab, "Orbitals")
    
    def load_json_data(self, data: Dict[str, Any]):
        """
        Load and display JSON data across all tabs.
        
        Args:
            data: Parsed QCORR JSON dictionary
        """
        self.json_data = data
        
        # Populate each tab with relevant data
        self.general_tab.populate(data)
        self.system_tab.populate(data)
        self.energy_tab.populate(data)
        self.vibration_tab.populate(data)
        self.orbital_tab.populate(data)


class GeneralInfoTab(QWidget):
    """Tab displaying general metadata (QM program, functional, basis set, etc.)"""
    
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        self.table = QTableWidget()
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Property", "Value"])
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)
    
    def populate(self, data: Dict[str, Any]):
        """Populate table with general metadata."""
        metadata = data.get("metadata", {})
        
        items = [
            ("QM Program", metadata.get("QM program", "N/A")),
            ("Run Date", metadata.get("run date", "N/A")),
            ("Functional", metadata.get("functional", "N/A")),
            ("Basis Set", metadata.get("basis set", "N/A")),
            ("Solvation", metadata.get("solvation", "none")),
            ("Dispersion Model", metadata.get("dispersion model", "none")),
            ("Grid Type", metadata.get("grid type", "N/A")),
            ("Ground/TS", metadata.get("ground or transition state", "N/A")),
            ("Keywords", metadata.get("keywords line", "N/A")),
        ]
        
        self.table.setRowCount(len(items))
        for i, (key, value) in enumerate(items):
            self.table.setItem(i, 0, QTableWidgetItem(key))
            self.table.setItem(i, 1, QTableWidgetItem(str(value)))


class SystemInfoTab(QWidget):
    """Tab displaying system information (charge, multiplicity, atoms, etc.)"""
    
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        self.table = QTableWidget()
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Property", "Value"])
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)
    
    def populate(self, data: Dict[str, Any]):
        """Populate table with system information."""
        rotational = data.get("rotational", {})
        
        items = [
            ("Charge", data.get("charge", "N/A")),
            ("Multiplicity", data.get("mult", "N/A")),
            ("Number of Atoms", data.get("natom", "N/A")),
            ("Number of Basis Functions", data.get("nbasis", "N/A")),
            ("Number of MOs", data.get("nmo", "N/A")),
            ("Temperature (K)", data.get("temperature", "N/A")),
            ("Pressure (atm)", data.get("pressure", "N/A")),
            ("Symmetry Point Group", rotational.get("symmetry point group", "N/A")),
            ("Symmetry Number", rotational.get("symmetry number", "N/A")),
            ("Optimization Converged", "Yes" if data.get("optdone", False) else "No"),
        ]
        
        self.table.setRowCount(len(items))
        for i, (key, value) in enumerate(items):
            self.table.setItem(i, 0, QTableWidgetItem(key))
            self.table.setItem(i, 1, QTableWidgetItem(str(value)))


class EnergyInfoTab(QWidget):
    """Tab displaying energetic properties."""
    
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        self.table = QTableWidget()
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Property", "Value (Hartree)"])
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)
    
    def populate(self, data: Dict[str, Any]):
        """Populate table with energy information."""
        scf_energy = data.get("scfenergies", [None])[-1]  # Latest SCF energy
        
        items = [
            ("SCF Energy", scf_energy),
            ("Zero-Point Vibrational Energy (ZPVE)", data.get("zpve", "N/A")),
            ("Enthalpy", data.get("enthalpy", "N/A")),
            ("Free Energy", data.get("freeenergy", "N/A")),
            ("Entropy (Hartree/K)", data.get("entropy", "N/A")),
        ]
        
        self.table.setRowCount(len(items))
        for i, (key, value) in enumerate(items):
            self.table.setItem(i, 0, QTableWidgetItem(key))
            value_str = f"{value:.6f}" if isinstance(value, (int, float)) else str(value)
            self.table.setItem(i, 1, QTableWidgetItem(value_str))
# aqmeasy/ui/QCORR_ui/QCORR_analysiswidget.py

from typing import Dict, Any, Optional, List
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTabWidget,
    QTableWidget, QTableWidgetItem, QLabel, QPushButton,
    QComboBox, QSlider, QGroupBox, QHeaderView, QSplitter,
    QDoubleSpinBox, QMessageBox
)
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtCore import Qt, Signal
import py3Dmol
import math


class VibrationTab(QWidget):
    """
    Tab for vibrational mode display and visualization.
    
    Features:
    - Table of frequencies with IR intensities
    - 3D molecular viewer with mode animation (py3Dmol)
    """
    
    def __init__(self):
        super().__init__()
        self.json_data: Optional[Dict[str, Any]] = None
        self.init_ui()
    
    def init_ui(self):
        """Initialize vibration tab UI."""
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # Splitter for table and viewer
        splitter = QSplitter(Qt.Orientation.Horizontal)
        layout.addWidget(splitter)
        
        # Left: Frequency table
        table_widget = QWidget()
        table_layout = QVBoxLayout()
        table_widget.setLayout(table_layout)
        
        table_layout.addWidget(QLabel("<b>Vibrational Frequencies</b>"))
        
        self.freq_table = QTableWidget()
        self.freq_table.setColumnCount(3)
        self.freq_table.setHorizontalHeaderLabels(["Mode", "Frequency (cm⁻¹)", "IR Intensity"])
        self.freq_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.freq_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.freq_table.itemSelectionChanged.connect(self.on_mode_selected)
        table_layout.addWidget(self.freq_table)
        
        splitter.addWidget(table_widget)
        
        # Right: 3D Viewer
        viewer_widget = QWidget()
        viewer_layout = QVBoxLayout()
        viewer_widget.setLayout(viewer_layout)
        
        
        # Controls
        controls_group = QGroupBox("Visualization Controls")
        controls_layout = QVBoxLayout()
        controls_group.setLayout(controls_layout)
        
        # Mode selection
        mode_control_layout = QHBoxLayout()
        mode_control_layout.addWidget(QLabel("Mode:"))
        self.mode_combo = QComboBox()
        self.mode_combo.currentIndexChanged.connect(self.on_mode_combo_changed)
        mode_control_layout.addWidget(self.mode_combo)
        controls_layout.addLayout(mode_control_layout)
        
        # Amplitude control
        amplitude_layout = QHBoxLayout()
        amplitude_layout.addWidget(QLabel("Amplitude:"))
        self.amplitude_spin = QDoubleSpinBox()
        self.amplitude_spin.setRange(0.1, 5.0)
        self.amplitude_spin.setValue(1.0)
        self.amplitude_spin.setSingleStep(0.1)
        self.amplitude_spin.valueChanged.connect(self.update_visualization)
        amplitude_layout.addWidget(self.amplitude_spin)
        controls_layout.addLayout(amplitude_layout)
        
        # Style selection
        style_layout = QHBoxLayout()
        style_layout.addWidget(QLabel("Style:"))
        self.style_selector = QComboBox()
        self.style_selector.addItems(['Stick', 'Ball and Stick', 'VdW Spheres'])
        self.style_selector.setCurrentText('Ball and Stick')
        self.style_selector.currentTextChanged.connect(self.update_visualization)
        style_layout.addWidget(self.style_selector)
        controls_layout.addLayout(style_layout)
        
        # Animation buttons
        button_layout = QHBoxLayout()
        self.animate_button = QPushButton("Animate Mode")
        self.animate_button.clicked.connect(self.animate_mode)
        button_layout.addWidget(self.animate_button)
        
        self.stop_button = QPushButton("Stop Animation")
        self.stop_button.clicked.connect(self.stop_animation)
        button_layout.addWidget(self.stop_button)
        controls_layout.addLayout(button_layout)
        
        viewer_layout.addWidget(controls_group)
        
        # Web viewer for py3Dmol
        self.web_view = QWebEngineView()
        viewer_layout.addWidget(self.web_view)
        
        splitter.addWidget(viewer_widget)
        splitter.setSizes([300, 500])
    
    def populate(self, data: Dict[str, Any]):
        """Populate vibration data."""
        self.json_data = data
        
        vibfreqs = data.get("vibfreqs", [])
        vibirs = data.get("vibirs", [])
        vibsyms = data.get("vibsyms", [])
        
        # Populate frequency table
        self.freq_table.setRowCount(len(vibfreqs))
        for i, freq in enumerate(vibfreqs):
            ir_intensity = vibirs[i] if i < len(vibirs) else "N/A"
            symmetry = vibsyms[i] if i < len(vibsyms) else ""
            
            mode_label = f"{i+1} ({symmetry})" if symmetry else str(i+1)
            
            self.freq_table.setItem(i, 0, QTableWidgetItem(mode_label))
            self.freq_table.setItem(i, 1, QTableWidgetItem(f"{freq:.2f}"))
            
            ir_str = f"{ir_intensity:.4f}" if isinstance(ir_intensity, (int, float)) else str(ir_intensity)
            self.freq_table.setItem(i, 2, QTableWidgetItem(ir_str))
        
        # Populate mode dropdown
        self.mode_combo.clear()
        for i, freq in enumerate(vibfreqs):
            symmetry = vibsyms[i] if i < len(vibsyms) else ""
            label = f"Mode {i+1} ({symmetry}): {freq:.2f} cm⁻¹" if symmetry else f"Mode {i+1}: {freq:.2f} cm⁻¹"
            self.mode_combo.addItem(label)
        
        # Initialize viewer with structure
        if vibfreqs:
            self.update_visualization()
    
    def on_mode_selected(self):
        """Handle mode selection from table."""
        selected_rows = self.freq_table.selectionModel().selectedRows()
        if selected_rows:
            mode_idx = selected_rows[0].row()
            self.mode_combo.setCurrentIndex(mode_idx)
    
    def on_mode_combo_changed(self, index: int):
        """Handle mode selection from combo box."""
        if index >= 0:
            self.freq_table.selectRow(index)
            self.update_visualization()
    
    def update_visualization(self):
        """Update the 3D molecular visualization with static structure and arrows."""
        if not self.json_data:
            return
        
        mode_idx = self.mode_combo.currentIndex()
        if mode_idx < 0:
            return
        
        try:
            self.render_static_mode(mode_idx)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error rendering molecule: {str(e)}")
    
    def render_static_mode(self, mode_idx: int):
        """Render static molecule with displacement arrows."""
        amplitude = self.amplitude_spin.value()
        style = self.style_selector.currentText()
        
        atomcoords = self.json_data.get("atomcoords", [[]])
        atomnos = self.json_data.get("atomnos", [])
        vibdisps = self.json_data.get("vibdisps", [])
        
        if not atomcoords or not atomnos or mode_idx >= len(vibdisps):
            return
        
        # Use the last geometry (optimized structure)
        coords = atomcoords[-1] if len(atomcoords) > 0 else []
        displacements = vibdisps[mode_idx]

        n_atoms = min(len(coords), len(atomnos), len(displacements))
        if n_atoms == 0:
            return

        coords = coords[:n_atoms]
        atomnos = atomnos[:n_atoms]
        displacements = displacements[:n_atoms]
        
        # Create XYZ block
        xyz_block = self._create_xyz_block(coords, atomnos)
        
        # Create py3Dmol viewer
        viewer = py3Dmol.view(width='100%', height='100%')
        viewer.addModel(xyz_block, 'xyz')
        
        # Apply style
        if style == 'Stick':
            viewer.setStyle({'stick': {}})
        elif style == 'Ball and Stick':
            viewer.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.3}})
        elif style == 'VdW Spheres':
            viewer.setStyle({'stick': {}, 'sphere': {}})
        
        # Add displacement arrows
        for coord, disp in zip(coords, displacements):
            start = {'x': coord[0], 'y': coord[1], 'z': coord[2]}
            scale = amplitude * 0.5
            end = {
                'x': coord[0] + disp[0] * scale,
                'y': coord[1] + disp[1] * scale,
                'z': coord[2] + disp[2] * scale
            }
            viewer.addArrow({
                'start': start,
                'end': end,
                'radius': 0.15,
                'radiusRatio': 1.5,
                'mid': 0.8,
                'color': 'red'
            })
        
        viewer.zoomTo()
        html = viewer._make_html()
        self.web_view.setHtml(html)
    
    def animate_mode(self):
        """Animate the selected vibrational mode."""
        mode_idx = self.mode_combo.currentIndex()
        if mode_idx < 0:
            return
        
        try:
            self.render_animated_mode(mode_idx)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error animating mode: {str(e)}")
    
    def render_animated_mode(self, mode_idx: int):
        """Render animated vibrational mode using multi-frame XYZ data."""
        amplitude = self.amplitude_spin.value()
        style = self.style_selector.currentText()
        
        atomcoords = self.json_data.get("atomcoords", [[]])
        atomnos = self.json_data.get("atomnos", [])
        vibdisps = self.json_data.get("vibdisps", [])
        
        if not atomcoords or not atomnos or mode_idx >= len(vibdisps):
            return
        
        # Use the last geometry (optimized structure)
        coords = atomcoords[-1] if len(atomcoords) > 0 else []
        displacements = vibdisps[mode_idx]

        n_atoms = min(len(coords), len(atomnos), len(displacements))
        if n_atoms == 0:
            return

        coords = coords[:n_atoms]
        atomnos = atomnos[:n_atoms]
        displacements = displacements[:n_atoms]
        
        # Generate animation frames
        n_frames = 30
        xyz_frames = []
        
        for frame in range(n_frames):
            # Calculate phase: oscillate smoothly using cosine
            phase = math.cos(2 * math.pi * frame / n_frames)
            scale = amplitude * 0.3 * phase
            
            # Create displaced coordinates
            displaced_coords = [
                [c + d * scale for c, d in zip(coord, disp)]
                for coord, disp in zip(coords, displacements)
            ]

            xyz_frames.append(self._create_xyz_block(displaced_coords, atomnos))

        models_xyz = "\n".join(xyz_frames) + "\n"
        
        # Create py3Dmol viewer with animation
        viewer = py3Dmol.view(width='100%', height='100%')
        viewer.addModelsAsFrames(models_xyz, 'xyz')
        
        # Apply style to all frames
        if style == 'Stick':
            viewer.setStyle({}, {'stick': {}})
        elif style == 'Ball and Stick':
            viewer.setStyle({}, {'stick': {'radius': 0.15}, 'sphere': {'scale': 0.3}})
        elif style == 'VdW Spheres':
            viewer.setStyle({}, {'stick': {}, 'sphere': {}})
        
        viewer.zoomTo()
        viewer.animate({'loop': 'forward'})
        
        html = viewer._make_html()
        self.web_view.setHtml(html)
    
    def stop_animation(self):
        """Stop animation and return to static view."""
        self.update_visualization()
    
    def _create_xyz_block(self, coords: List[List[float]], atomnos: List[int]) -> str:
        """
        Create XYZ format block from coordinates and atomic numbers.
        
        Args:
            coords: List of [x, y, z] coordinates
            atomnos: List of atomic numbers
        
        Returns:
            XYZ format string
        """
        # Atomic number to symbol mapping
        atomic_numbers_to_symbols = {
            1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl',
            35: 'Br', 53: 'I', 11: 'Na', 12: 'Mg', 19: 'K', 20: 'Ca', 
            5: 'B', 14: 'Si', 26: 'Fe', 29: 'Cu', 30: 'Zn'
        }
        
        xyz_lines = []

        atom_lines = []
        for coord, atno in zip(coords, atomnos):
            if len(coord) < 3:
                continue
            symbol = atomic_numbers_to_symbols.get(int(atno), 'X')
            x, y, z = coord
            atom_lines.append(f"{symbol} {x:.6f} {y:.6f} {z:.6f}")

        xyz_lines.append(f"{len(atom_lines)}")
        xyz_lines.append("")
        xyz_lines.extend(atom_lines)
        
        return '\n'.join(xyz_lines)

class OrbitalInfoTab(QWidget):
    """Tab displaying molecular orbital information."""
    
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # Summary table
        summary_label = QLabel("<b>Orbital Summary</b>")
        layout.addWidget(summary_label)
        
        self.summary_table = QTableWidget()
        self.summary_table.setColumnCount(2)
        self.summary_table.setHorizontalHeaderLabels(["Property", "Value"])
        self.summary_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.summary_table)
        
        # Energy levels table
        energies_label = QLabel("<b>Orbital Energies (eV)</b>")
        layout.addWidget(energies_label)
        
        self.energies_table = QTableWidget()
        self.energies_table.setColumnCount(3)
        self.energies_table.setHorizontalHeaderLabels(["Orbital", "Symmetry", "Energy (eV)"])
        layout.addWidget(self.energies_table)
    
    def populate(self, data: Dict[str, Any]):
        """Populate orbital information."""
        homos = data.get("homos", [])
        moenergies = data.get("moenergies", [[]])
        mosyms = data.get("mosyms", [[]])
        
        # Summary
        homo_idx = homos[0] if homos else None
        lumo_idx = homo_idx + 1 if homo_idx is not None else None
        
        energies = moenergies[0] if moenergies else []
        homo_energy = energies[homo_idx] if homo_idx is not None and homo_idx < len(energies) else None
        lumo_energy = energies[lumo_idx] if lumo_idx is not None and lumo_idx < len(energies) else None
        gap = (lumo_energy - homo_energy) if (homo_energy is not None and lumo_energy is not None) else None
        
        summary_items = [
            ("HOMO Index", homo_idx if homo_idx is not None else "N/A"),
            ("LUMO Index", lumo_idx if lumo_idx is not None else "N/A"),
            ("HOMO Energy (eV)", f"{homo_energy:.4f}" if homo_energy is not None else "N/A"),
            ("LUMO Energy (eV)", f"{lumo_energy:.4f}" if lumo_energy is not None else "N/A"),
            ("HOMO-LUMO Gap (eV)", f"{gap:.4f}" if gap is not None else "N/A"),
        ]
        
        self.summary_table.setRowCount(len(summary_items))
        for i, (key, value) in enumerate(summary_items):
            self.summary_table.setItem(i, 0, QTableWidgetItem(key))
            self.summary_table.setItem(i, 1, QTableWidgetItem(str(value)))
        
        # Orbital energies (show HOMO-5 to LUMO+5)
        if energies and homo_idx is not None:
            start_idx = max(0, homo_idx - 5)
            end_idx = min(len(energies), homo_idx + 6)
            
            syms = mosyms[0] if mosyms else []
            
            display_orbitals = []
            for i in range(start_idx, end_idx):
                orbital_type = "HOMO" if i == homo_idx else ("LUMO" if i == lumo_idx else "")
                label = f"{i+1} ({orbital_type})" if orbital_type else str(i+1)
                sym = syms[i] if i < len(syms) else ""
                energy = energies[i]
                display_orbitals.append((label, sym, energy))
            
            self.energies_table.setRowCount(len(display_orbitals))
            for i, (label, sym, energy) in enumerate(display_orbitals):
                self.energies_table.setItem(i, 0, QTableWidgetItem(label))
                self.energies_table.setItem(i, 1, QTableWidgetItem(sym))
                self.energies_table.setItem(i, 2, QTableWidgetItem(f"{energy:.4f}"))