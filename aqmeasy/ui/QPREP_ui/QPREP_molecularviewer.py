import os
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout,
    QLabel,  QComboBox,
    QCheckBox,
    QMessageBox, QGroupBox,  QApplication, QSlider, QFormLayout, QSizePolicy, QTextBrowser
)
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtCore import  Qt
from rdkit import Chem
from rdkit.Chem import rdMolAlign
import py3Dmol

from aqmeasy.ui.stylesheets import stylesheets
from aqmeasy.utils import smiles2multiplicity


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
        self.setup_ui()
        # Initially hide the 3D viewer group
        self.web_view.setVisible(False)

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

        order_row = QHBoxLayout()
        self.order_by_rmsd_checkbox = QCheckBox("Order overlay by RMSD to reference")
        self.order_by_rmsd_checkbox.setEnabled(False)
        self.order_by_rmsd_checkbox.toggled.connect(self.render_selected_molecule)
        order_row.addWidget(self.order_by_rmsd_checkbox)
        display_layout.addLayout(order_row)

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
                                'file_basename': os.path.basename(filename)
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
                                'file_basename': os.path.basename(filename)
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
        self.molecules = []
        self.loaded_files = []
        self.file_info_label.setText("No files selected")
        self.molecule_selector.clear()
        self.ensemble_overlay_checkbox.setChecked(False)
        self.order_by_rmsd_checkbox.setChecked(False)
        self.order_by_rmsd_checkbox.setEnabled(False)
        self.web_view.setHtml("")
        self.web_view.setVisible(False)
        self.formula_label.setText("N/A")
        self.weight_label.setText("N/A")
        self.atoms_label.setText("N/A")
        self.bonds_label.setText("N/A")
        self.charge_label.setText("N/A")
        self.mult_label.setText("N/A")
        self.source_file_label.setText("N/A")

    def on_molecule_change(self, index):
        if index < 0 or not self.molecules or index >= len(self.molecules):
            self.formula_label.setText("N/A")
            self.weight_label.setText("N/A")
            self.atoms_label.setText("N/A")
            self.bonds_label.setText("N/A")
            self.charge_label.setText("N/A")
            self.mult_label.setText("N/A")
            self.source_file_label.setText("N/A")
            return

        self.list_slider.setValue(index)
        self.molecule_selector.setCurrentIndex(index)

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

            QApplication.processEvents()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error getting molecule details: {str(e)}")
            self.formula_label.setText("Error")
            self.weight_label.setText("Error")
            self.atoms_label.setText("Error")
            self.bonds_label.setText("Error")
            self.charge_label.setText("Error")
            self.mult_label.setText("Error")
            self.source_file_label.setText("Error")

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

        viewer = py3Dmol.view(width='100%', height='100%')

        if self.ensemble_overlay_checkbox.isChecked() and not xyz:
            ensemble = [entry for entry in self.molecules if entry['source_file'] == mol_data['source_file']]
            ordered_ensemble = self._prepare_aligned_overlay(ensemble, mol)
            rmsd_sorted = self.order_by_rmsd_checkbox.isChecked()
            for model_idx, ensemble_entry in enumerate(ordered_ensemble):
                model_block = Chem.MolToMolBlock(ensemble_entry['mol'])
                viewer.addModel(model_block, 'mol')
                style_spec = self._style_for_overlay_model(style, model_idx, rmsd_sorted)
                viewer.setStyle({'model': model_idx}, style_spec)

            viewer.zoomTo()
            self.web_view.setHtml(viewer._make_html())
            return

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

    def _on_overlay_toggled(self, checked):
        """Enable RMSD ordering only when overlay mode is active."""
        self.order_by_rmsd_checkbox.setEnabled(checked)
        if not checked and self.order_by_rmsd_checkbox.isChecked():
            self.order_by_rmsd_checkbox.blockSignals(True)
            self.order_by_rmsd_checkbox.setChecked(False)
            self.order_by_rmsd_checkbox.blockSignals(False)
        self.render_selected_molecule()

    def _style_for_overlay_model(self, style, model_idx, rmsd_sorted=False):
        """Return a py3Dmol style for one model in an ensemble overlay."""
        if model_idx == 0:
            stick_radius = 0.22
            sphere_scale = 0.3
            color = '#1f77b4'
            opacity = 1.0
        else:
            stick_radius = 0.12
            sphere_scale = 0.2
            if rmsd_sorted:
                palette = ['#2ecc71', '#27ae60', '#f1c40f', '#e67e22', '#e74c3c', '#8e44ad']
                color = palette[(model_idx - 1) % len(palette)]
            else:
                color = '#95a5a6'
            opacity = 0.45

        if style == 'Stick':
            return {'stick': {'radius': stick_radius, 'color': color, 'opacity': opacity}}
        if style == 'Ball and Stick':
            return {
                'stick': {'radius': stick_radius, 'color': color, 'opacity': opacity},
                'sphere': {'scale': sphere_scale, 'color': color, 'opacity': opacity},
            }
        if style == 'VdW Spheres':
            return {'sphere': {'color': color, 'opacity': opacity}}
        if style == 'Surface':
            return {'stick': {'radius': stick_radius, 'color': color, 'opacity': opacity}}
        return {'stick': {'radius': stick_radius, 'color': color, 'opacity': opacity}}

    def _prepare_aligned_overlay(self, ensemble_entries, reference_mol):
        """Align overlay molecules to the reference, and optionally sort by RMSD."""
        if not ensemble_entries:
            return []

        overlay_entries = []
        for entry in ensemble_entries:
            entry_mol = Chem.Mol(entry['mol'])
            rmsd = float('inf')
            try:
                rmsd = float(rdMolAlign.GetBestRMS(reference_mol, entry_mol))
                rdMolAlign.AlignMol(entry_mol, reference_mol)
            except Exception:
                pass
            overlay_entries.append({'mol': entry_mol, 'rmsd': rmsd})

        if self.order_by_rmsd_checkbox.isChecked():
            overlay_entries.sort(key=lambda item: item['rmsd'])

        return overlay_entries