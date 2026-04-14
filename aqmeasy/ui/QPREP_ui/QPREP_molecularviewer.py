import os
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout,
    QLabel,  QComboBox,
    QCheckBox,
    QMessageBox, QGroupBox, QSlider, QFormLayout, QSizePolicy, QTextBrowser
)
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtCore import Qt, Signal, Slot, QObject, QThread
from rdkit import Chem
from rdkit.Chem import rdMolAlign
import py3Dmol

from aqmeasy.ui.stylesheets import stylesheets
from aqmeasy.utils import smiles2multiplicity


def build_aligned_overlay_entries(ensemble_entries, reference_entry):
    """Build aligned overlay blocks and RMSD values for one reference conformer."""
    if not ensemble_entries:
        return []

    reference_mol = reference_entry['mol']
    reference_atom_count = reference_mol.GetNumAtoms()
    reference_symbols = [atom.GetSymbol() for atom in reference_mol.GetAtoms()]
    atom_map = [(idx, idx) for idx in range(reference_atom_count)]
    overlay_entries = []

    for entry in ensemble_entries:
        entry_mol = Chem.Mol(entry['mol'])
        rmsd = float('inf')
        try:
            same_atom_order = entry_mol.GetNumAtoms() == reference_atom_count and [atom.GetSymbol() for atom in entry_mol.GetAtoms()] == reference_symbols
            if same_atom_order:
                rmsd = float(rdMolAlign.AlignMol(entry_mol, reference_mol, atomMap=atom_map))
            else:
                rmsd = float(rdMolAlign.GetBestRMS(reference_mol, entry_mol))
                rdMolAlign.AlignMol(entry_mol, reference_mol)
        except Exception:
            pass

        try:
            model_block = Chem.MolToMolBlock(entry_mol)
        except Exception:
            model_block = Chem.MolToMolBlock(entry['mol'])

        overlay_entries.append({'model_block': model_block, 'rmsd': rmsd})

    return overlay_entries


def compute_pairwise_average_rmsd(ensemble_entries):
    """Compute the mean RMSD across all conformer pairs in an ensemble."""
    if len(ensemble_entries) < 2:
        return None

    reference_mol = ensemble_entries[0]['mol']
    atom_count = reference_mol.GetNumAtoms()
    reference_symbols = [atom.GetSymbol() for atom in reference_mol.GetAtoms()]
    atom_map = [(idx, idx) for idx in range(atom_count)]

    rmsd_values = []
    for left_index in range(len(ensemble_entries) - 1):
        left_mol = ensemble_entries[left_index]['mol']
        left_symbols = [atom.GetSymbol() for atom in left_mol.GetAtoms()]

        for right_index in range(left_index + 1, len(ensemble_entries)):
            right_mol = Chem.Mol(ensemble_entries[right_index]['mol'])
            try:
                same_atom_order = (
                    left_mol.GetNumAtoms() == atom_count
                    and right_mol.GetNumAtoms() == atom_count
                    and left_symbols == reference_symbols
                    and [atom.GetSymbol() for atom in right_mol.GetAtoms()] == reference_symbols
                )
                if same_atom_order:
                    rmsd = float(rdMolAlign.AlignMol(right_mol, left_mol, atomMap=atom_map))
                else:
                    rmsd = float(rdMolAlign.GetBestRMS(left_mol, right_mol))
                rmsd_values.append(rmsd)
            except Exception:
                continue

    if not rmsd_values:
        return None

    return sum(rmsd_values) / len(rmsd_values)


class OverlayPreparationWorker(QObject):
    """Background worker that prepares overlay blocks and RMSD values."""

    finished = Signal(int, object, object, object)
    failed = Signal(int, object, str)

    def __init__(self, request_id, cache_key, ensemble_entries, reference_entry):
        super().__init__()
        self.request_id = request_id
        self.cache_key = cache_key
        self.ensemble_entries = ensemble_entries
        self.reference_entry = reference_entry

    @Slot()
    def run(self):
        thread = self.thread()

        try:
            if thread and thread.isInterruptionRequested():
                return

            overlay_entries = build_aligned_overlay_entries(self.ensemble_entries, self.reference_entry)
            pairwise_average_rmsd = compute_pairwise_average_rmsd(self.ensemble_entries)

            if thread and thread.isInterruptionRequested():
                return

            self.finished.emit(self.request_id, self.cache_key, overlay_entries, pairwise_average_rmsd)
        except Exception as error:
            self.failed.emit(self.request_id, self.cache_key, str(error))


class MoleculeViewer(QWidget):
    """
    Main widget for Molecule Viewer application.
    Displays molecules from SDF files, with the 3D view integrated directly.
    """

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Molecule Viewer")
        self.molecules = []
        self.loaded_files = []  # Track which files are loaded
        self._overlay_cache = {}
        self._overlay_summary_cache = {}
        self._syncing_selection = False
        self._overlay_thread = None
        self._overlay_worker = None
        self._overlay_request_id = 0
        self._overlay_pending_key = None
        self.setup_ui()
        # Initially hide the 3D viewer group
        self.web_view.setVisible(False)

    def _stop_overlay_preparation(self, wait=False):
        thread = self._overlay_thread
        if thread is not None and thread.isRunning():
            thread.requestInterruption()
            thread.quit()
            if wait:
                thread.wait(1000)

        if wait or thread is None or not thread.isRunning():
            self._overlay_thread = None
            self._overlay_worker = None
            self._overlay_pending_key = None

    def setup_ui(self):
        self.options_container = QWidget()
        options_layout = QVBoxLayout(self.options_container)
        options_layout.setContentsMargins(0, 0, 0, 0)
        options_layout.setSpacing(10)

        # File Information Section
        self.file_group = QGroupBox("File Information")
        file_layout = QVBoxLayout()
        self.file_info_label = QTextBrowser()
        self.file_info_label.setPlaceholderText("No files selected")
        self.file_info_label.setMinimumHeight(100)
        self.file_info_label.setMaximumHeight(150)
        file_layout.addWidget(self.file_info_label)
        self.file_group.setLayout(file_layout)
        options_layout.addWidget(self.file_group)

        # Display Options Section
        display_group = QGroupBox("Display Options")
        display_layout = QVBoxLayout()

        style_row = QHBoxLayout()
        style_row.addWidget(QLabel("Display Style:"))
        self.style_selector = QComboBox()
        self.style_selector.addItems(['Stick', 'Ball and Stick', 'VdW Spheres', 'Surface'])
        self.style_selector.currentIndexChanged.connect(self.render_selected_molecule)
        style_row.addWidget(self.style_selector)
        display_layout.addLayout(style_row)

        molecule_row = QHBoxLayout()
        molecule_row.addWidget(QLabel("Molecule:"))
        self.molecule_selector = QComboBox()
        self.molecule_selector.currentIndexChanged.connect(self.on_molecule_change)
        self.molecule_selector.setMaxVisibleItems(10)
        molecule_row.addWidget(self.molecule_selector)
        display_layout.addLayout(molecule_row)

        # List slider to navigate through the molecules to display
        list_slider_row = QHBoxLayout()
        list_slider_row.addWidget(QLabel("Molecule List:"))
        self.list_slider = QSlider(Qt.Horizontal)
        self.list_slider.setMinimum(0)
        self.list_slider.setMaximum(0)  # until file uploaded
        self.list_slider.setValue(0)
        self.list_slider.setTickInterval(1)
        self.list_slider.setSingleStep(1)
        self.list_slider.valueChanged.connect(self.on_molecule_change)
        self.list_slider.sliderReleased.connect(self.render_selected_molecule)

        list_slider_row.addWidget(self.list_slider)
        display_layout.addLayout(list_slider_row)

        overlay_row = QHBoxLayout()
        self.ensemble_overlay_checkbox = QCheckBox("Overlay ensemble from selected file")
        self.ensemble_overlay_checkbox.toggled.connect(self._on_overlay_toggled)
        overlay_row.addWidget(self.ensemble_overlay_checkbox)
        display_layout.addLayout(overlay_row)

        self.overlay_summary_label = QLabel("Average pairwise RMSD: N/A")
        self.overlay_summary_label.setWordWrap(True)
        display_layout.addWidget(self.overlay_summary_label)

        display_group.setLayout(display_layout)
        options_layout.addWidget(display_group)

        # Molecular Information Section
        self.info_group = QGroupBox("Molecular Information")
        self.info_layout = QFormLayout()

        self.source_file_label = QLabel("N/A")
        self.formula_label = QLabel("N/A")  # all states N/A before file imported.
        self.weight_label = QLabel("N/A")
        self.atoms_label = QLabel("N/A")
        self.bonds_label = QLabel("N/A")
        self.charge_label = QLabel("N/A")
        self.mult_label = QLabel("N/A")


        self.info_layout.addRow("Source File:", self.source_file_label)
        self.info_layout.addRow("Formula:", self.formula_label)
        self.info_layout.addRow("Weight:", self.weight_label)
        self.info_layout.addRow("Atoms:", self.atoms_label)
        self.info_layout.addRow("Bonds:", self.bonds_label)
        self.info_layout.addRow("Charge:", self.charge_label)
        self.info_layout.addRow("Multiplicity:", self.mult_label)

        self.info_group.setLayout(self.info_layout)
        options_layout.addWidget(self.info_group)

        options_layout.addStretch()

        self.viewer_group = QGroupBox("3D Viewer")
        viewer_layout = QVBoxLayout()
        self.web_view = QWebEngineView()
        self.web_view.setVisible(False)
        self.web_view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        viewer_layout.addWidget(self.web_view)
        self.viewer_group.setLayout(viewer_layout)

        # Main horizontal layout
        main_layout = QHBoxLayout(self)
        main_layout.addWidget(self.options_container)
        main_layout.addWidget(self.viewer_group, 1)

    def load_molecules_from_files(self, filenames):
        """Load molecules from multiple files"""
        if not filenames:
            self.clear_molecules()
            return

        # Filter for SDF files only
        sdf_files = [f for f in filenames if f.lower().endswith('.sdf')]
        
        # Filter for XYZ files
        xyz_files = [f for f in filenames if f.lower().endswith('.xyz')]

        if not sdf_files and not xyz_files:
            self.clear_molecules()
            return

        try:
            self.molecules = []
            self.loaded_files = []
            self._overlay_cache = {}
            self._overlay_summary_cache = {}
            total_molecules = 0

            for filename in sdf_files:
                try:
                    supplier = Chem.SDMolSupplier(filename, removeHs=False)
                    file_molecules = []
                    
                    for i, mol in enumerate(supplier):
                        if mol is not None:
                            # Store molecule with source file information
                            mol_data = {
                                'mol': mol,
                                'original_idx': i,
                                'source_file': filename,
                                    'file_basename': os.path.basename(filename),
                                    'global_idx': len(self.molecules),
                            }
                            file_molecules.append(mol_data)
                            self.molecules.append(mol_data)
                    
                    if file_molecules:
                        self.loaded_files.append({
                            'filename': filename,
                            'basename': os.path.basename(filename),
                            'molecule_count': len(file_molecules)
                        })
                        total_molecules += len(file_molecules)
                        
                except Exception as file_error:
                    print(f"Error loading SDF file {filename}: {str(file_error)}")
                    continue

            for filename in xyz_files:
                try:
                    with open(filename, 'r') as f:
                        lines = f.readlines()
                        num_atoms = int(lines[0].strip())
                        mol_block = ''.join(lines[0:2+num_atoms])
                        mol = Chem.MolFromXYZBlock(mol_block)

                        if mol is not None:
                            Chem.SanitizeMol(mol)
                            mol_data = {
                                'mol': mol,
                                'original_idx': 0,
                                'source_file': filename,
                                'file_basename': os.path.basename(filename),
                                'global_idx': len(self.molecules),
                            }
                            self.molecules.append(mol_data)
                            self.loaded_files.append({
                                'filename': filename,
                                'basename': os.path.basename(filename),
                                'molecule_count': 1
                            })
                            total_molecules += 1

                except Exception as file_error:
                    print(f"Error loading XYZ file {filename}: {str(file_error)}")
                    continue

            if not self.molecules:
                QMessageBox.warning(self, "Error", "No valid molecules found in the selected files.")
                self.clear_molecules()
                return

            self.web_view.setVisible(True)
            
            # Update file info
            file_info = f"Loaded {len(sdf_files)} SDF file(s):\n"
            for file_data in self.loaded_files:
                file_info += f"• {file_data['basename']}: {file_data['molecule_count']} molecules\n"
            file_info += f"\nTotal molecules: {total_molecules}"
            self.file_info_label.setText(file_info)
            
            self.list_slider.setMaximum(len(self.molecules) - 1)

            # Update molecule selector
            self.molecule_selector.currentIndexChanged.disconnect()
            self.molecule_selector.clear()

            for i, mol_data in enumerate(self.molecules):
                mol = mol_data['mol']
                original_idx = mol_data['original_idx']
                
                mol_name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"Mol {original_idx + 1}"
                if not mol_name.strip():
                    mol_name = f"Mol {original_idx + 1}"

                try:
                    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                    label = f"{mol_name} ({formula})"
                except:
                    label = f"{mol_name}"

                self.molecule_selector.addItem(label)

            self.molecule_selector.currentIndexChanged.connect(self.on_molecule_change)

            if self.molecules:
                self.molecule_selector.setCurrentIndex(0)
                self.on_molecule_change(0)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load molecules from files.\n\n{str(e)}")
            self.clear_molecules()
            self.web_view.setVisible(False)

    def load_molecules_from_file(self, filename):
        """
        Loads molecules from a single file, clearing previous molecules.
        This method is designed to handle a single file selection from the list.
        """
        self.clear_molecules()
        if filename:
            self.load_molecules_from_files([filename])
        else:
            self.clear_molecules()


    def clear_molecules(self):
        self._stop_overlay_preparation()
        self.molecules = []
        self.loaded_files = []
        self._overlay_cache = {}
        self._overlay_summary_cache = {}
        self.file_info_label.setText("No files selected")
        self.molecule_selector.clear()
        self.ensemble_overlay_checkbox.setChecked(False)
        self.web_view.setHtml("")
        self.web_view.setVisible(False)
        self.overlay_summary_label.setText("Average pairwise RMSD: N/A")
        self.formula_label.setText("N/A")
        self.weight_label.setText("N/A")
        self.atoms_label.setText("N/A")
        self.bonds_label.setText("N/A")
        self.charge_label.setText("N/A")
        self.mult_label.setText("N/A")
        self.source_file_label.setText("N/A")

    def on_molecule_change(self, index):
        if self._syncing_selection:
            return

        if index < 0 or not self.molecules or index >= len(self.molecules):
            self.formula_label.setText("N/A")
            self.weight_label.setText("N/A")
            self.atoms_label.setText("N/A")
            self.bonds_label.setText("N/A")
            self.charge_label.setText("N/A")
            self.mult_label.setText("N/A")
            self.source_file_label.setText("N/A")
            return

        self._syncing_selection = True
        try:
            self.list_slider.blockSignals(True)
            self.molecule_selector.blockSignals(True)
            self.list_slider.setValue(index)
            self.molecule_selector.setCurrentIndex(index)
        finally:
            self.list_slider.blockSignals(False)
            self.molecule_selector.blockSignals(False)
            self._syncing_selection = False

        mol_data = self.molecules[index]
        mol = mol_data['mol']
        original_idx = mol_data['original_idx']
        source_file = mol_data['file_basename']

        try:
            # Display basic molecule info
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            mol_weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)
            num_atoms = mol.GetNumAtoms()
            num_bonds = mol.GetNumBonds()
            charge = Chem.GetFormalCharge(mol)
            mult = smiles2multiplicity(Chem.MolToSmiles(mol))

            mol_name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"Molecule {original_idx + 1}"
            if not mol_name.strip():
                mol_name = f"Molecule {original_idx + 1}"

            self.formula_label.setText(formula)
            self.weight_label.setText(f"{mol_weight:.2f} g/mol")
            self.atoms_label.setText(str(num_atoms))
            self.bonds_label.setText(str(num_bonds))
            self.charge_label.setText(str(charge))
            self.mult_label.setText(str(mult))
            self.source_file_label.setText(source_file)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error getting molecule details: {str(e)}")
            self.formula_label.setText("Error")
            self.weight_label.setText("Error")
            self.atoms_label.setText("Error")
            self.bonds_label.setText("Error")
            self.charge_label.setText("Error")
            self.mult_label.setText("Error")
            self.source_file_label.setText("Error")

        if self.sender() is not self.list_slider:
            self.render_selected_molecule()

    def render_selected_molecule(self):
        if not self.molecules:
            return
        
        # check if selected molecules ends with xyz
        xyz = any(mol_data['source_file'].lower().endswith('.xyz') for mol_data in self.molecules)

        index = self.list_slider.value()
        if index < 0 or index >= len(self.molecules):
            return

        mol_data = self.molecules[index]
        mol = mol_data['mol']
        style = self.style_selector.currentText()

        if self.ensemble_overlay_checkbox.isChecked() and not xyz:
            ensemble = [entry for entry in self.molecules if entry['source_file'] == mol_data['source_file']]
            cache_key = self._overlay_cache_key(mol_data)
            ordered_ensemble = self._overlay_cache.get(cache_key)

            if ordered_ensemble is None:
                self._show_overlay_loading_message(mol_data)
                self._start_overlay_preparation(cache_key, ensemble, mol_data)
                return

            source_file = mol_data['source_file']
            pairwise_average_rmsd = self._overlay_summary_cache.get(source_file)
            if pairwise_average_rmsd is not None:
                self.overlay_summary_label.setText(f"Average pairwise RMSD: {pairwise_average_rmsd:.3f} Å")
            else:
                self.overlay_summary_label.setText("Average pairwise RMSD: N/A")

            self._render_overlay_entries(ordered_ensemble, style)
            return

        self.overlay_summary_label.setText("Average pairwise RMSD: N/A")
        self._render_single_molecule(mol, style, xyz)

    def _render_single_molecule(self, mol, style, xyz):
        viewer = py3Dmol.view(width='100%', height='100%')

        if xyz:
            try:
                xyz_block = Chem.MolToXYZBlock(mol)
                viewer.addModel(xyz_block, 'xyz')
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error converting molecule to XYZ format: {str(e)}")
                return
        else:
            mol_block = Chem.MolToMolBlock(mol) 
            viewer.addModel(mol_block, 'mol')

        if style == 'Stick':
            viewer.setStyle({'stick': {}})
        elif style == 'Ball and Stick':
            viewer.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.3}})
        elif style == 'VdW Spheres':
            viewer.setStyle({'stick': {}, 'sphere': {}})
        elif style == 'Surface':
            viewer.setStyle({'stick': {}})
            viewer.addSurface(py3Dmol.VDW, {'opacity': 0.8})

        viewer.zoomTo()
        html = viewer._make_html()
        self.web_view.setHtml(html)

    def _render_overlay_entries(self, overlay_entries, style):
        viewer = py3Dmol.view(width='100%', height='100%')
        ordered_entries = list(overlay_entries)

        for model_idx, ensemble_entry in enumerate(ordered_entries):
            viewer.addModel(ensemble_entry['model_block'], 'mol')
            style_spec = self._style_for_overlay_model(style)
            viewer.setStyle({'model': model_idx}, style_spec)

        viewer.zoomTo()
        self.web_view.setHtml(viewer._make_html())

    def _show_overlay_loading_message(self, mol_data):
        mol_name = mol_data['mol'].GetProp('_Name') if mol_data['mol'].HasProp('_Name') else 'selected molecule'
        if not str(mol_name).strip():
            mol_name = 'selected molecule'

        html = f"""
                <html>
                    <head>
                        <style>
                            html, body {{
                                margin: 0;
                                width: 100%;
                                height: 100%;
                                background: #f7f7f9;
                                font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
                            }}
                            .overlay-loading {{
                                width: 100%;
                                height: 100%;
                                display: flex;
                                align-items: center;
                                justify-content: center;
                                text-align: center;
                                color: #222;
                            }}
                            .panel {{
                                max-width: 28rem;
                                padding: 1.25rem 1.5rem;
                                border-radius: 14px;
                                background: rgba(255, 255, 255, 0.92);
                                box-shadow: 0 10px 30px rgba(0, 0, 0, 0.10);
                            }}
                            .title {{
                                font-size: 18px;
                                font-weight: 600;
                                margin-bottom: 0.45rem;
                            }}
                            .subtitle {{
                                font-size: 13px;
                                line-height: 1.45;
                                color: #555;
                            }}
                        </style>
                    </head>
                    <body>
                        <div class="overlay-loading">
                            <div class="panel">
                                <div class="title">Preparing overlay for {mol_name}</div>
                                <div class="subtitle">This can take a moment for large conformer sets.</div>
                            </div>
                        </div>
                    </body>
                </html>
                """
        self.web_view.setHtml(html)

    def _overlay_cache_key(self, mol_data):
        return (mol_data['source_file'], mol_data.get('global_idx', mol_data['original_idx']))

    def _start_overlay_preparation(self, cache_key, ensemble_entries, reference_entry):
        if self._overlay_pending_key == cache_key:
            return

        if self._overlay_thread is not None and self._overlay_thread.isRunning():
            self._overlay_thread.requestInterruption()

        self._overlay_request_id += 1
        request_id = self._overlay_request_id
        self._overlay_pending_key = cache_key

        thread = QThread()
        worker = OverlayPreparationWorker(request_id, cache_key, ensemble_entries, reference_entry)
        worker.moveToThread(thread)

        thread.started.connect(worker.run)
        worker.finished.connect(self._on_overlay_preparation_finished)
        worker.failed.connect(self._on_overlay_preparation_failed)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)

        self._overlay_thread = thread
        self._overlay_worker = worker
        thread.start()

    def _on_overlay_preparation_finished(self, request_id, cache_key, overlay_entries, pairwise_average_rmsd):
        if request_id != self._overlay_request_id:
            return

        self._overlay_cache[cache_key] = list(overlay_entries)
        self._overlay_summary_cache[cache_key[0]] = pairwise_average_rmsd
        if self._overlay_pending_key == cache_key:
            self._overlay_pending_key = None

        self._overlay_thread = None
        self._overlay_worker = None

        if not self.molecules:
            return

        index = self.list_slider.value()
        if index < 0 or index >= len(self.molecules):
            return

        current_key = self._overlay_cache_key(self.molecules[index])
        if current_key == cache_key and self.ensemble_overlay_checkbox.isChecked():
            self.render_selected_molecule()

    def _on_overlay_preparation_failed(self, request_id, cache_key, error_message):
        if request_id != self._overlay_request_id:
            return

        if self._overlay_pending_key == cache_key:
            self._overlay_pending_key = None

        self._overlay_thread = None
        self._overlay_worker = None
        QMessageBox.critical(self, "Error", f"Failed to prepare overlay ensemble.\n\n{error_message}")

    def _on_overlay_toggled(self, checked):
        """Render the overlay state when overlay mode changes."""
        self.render_selected_molecule()

    def _style_for_overlay_model(self, style):
        """Return a py3Dmol style for one model in an ensemble overlay."""
        opacity = 0.85

        if style == 'Stick':
            return {'stick': {'opacity': opacity}}
        if style == 'Ball and Stick':
            return {
                'stick': {'radius': 0.15, 'opacity': opacity},
                'sphere': {'scale': 0.3, 'opacity': opacity},
            }
        if style == 'VdW Spheres':
            return {'sphere': {'opacity': opacity}}
        if style == 'Surface':
            return {'stick': {'opacity': opacity}}
        return {'stick': {'opacity': opacity}}

    def _prepare_aligned_overlay(self, ensemble_entries, reference_entry):
        """Align overlay molecules to the reference, and optionally sort by RMSD."""
        if not ensemble_entries:
            return []

        cache_key = self._overlay_cache_key(reference_entry)
        if cache_key not in self._overlay_cache:
            self._overlay_cache[cache_key] = build_aligned_overlay_entries(ensemble_entries, reference_entry)

        return list(self._overlay_cache.get(cache_key, []))

    def closeEvent(self, event):
        self._stop_overlay_preparation(wait=True)
        super().closeEvent(event)