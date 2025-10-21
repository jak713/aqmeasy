import os
import sys
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QComboBox, QDialog,
    QMessageBox, QFileDialog, QGroupBox, QTextEdit, QApplication, QSlider, QFormLayout
)
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtCore import QEventLoop, QUrl, Signal, Qt
from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter
import py3Dmol
import math


class MoleculeViewer(QWidget):
    """
    Main widget for Molecule Viewer application.
    Displays molecules from SDF files, with the 3D view integrated directly.
    Now supports multiple files.
    """

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Molecule Viewer")
        self.molecules = []
        self.loaded_files = []  # Track which files are loaded
        self.setup_ui()
        # Initially hide the 3D viewer group
        self.viewer_group.setVisible(False)

    def setup_ui(self):
        # This container holds all the options, but NOT the 3D viewer itself
        self.options_container = QWidget()
        options_layout = QVBoxLayout(self.options_container)
        options_layout.setContentsMargins(0, 0, 0, 0)
        options_layout.setSpacing(10)

        # File Information Section
        file_group = QGroupBox("File Information")
        file_layout = QVBoxLayout()
        self.file_info_label = QLabel("No files selected")
        self.file_info_label.setWordWrap(True)
        file_layout.addWidget(self.file_info_label)
        file_group.setLayout(file_layout)
        options_layout.addWidget(file_group)

        # Display Options Section
        display_group = QGroupBox("Display Options")
        display_layout = QVBoxLayout()

        combo_style = """
            QComboBox { 
                background-color: #edf2f4; 
                color: #2b2d42; 
                border: 1px solid #8d99ae; 
                border-radius: 3px; 
                padding: 2px; 
                combobox-popup: 0;
            }
        """

        style_row = QHBoxLayout()
        style_row.addWidget(QLabel("Display Style:"))
        self.style_selector = QComboBox()
        self.style_selector.setStyleSheet(combo_style)
        self.style_selector.addItems(['Stick', 'Ball and Stick', 'Surface'])
        self.style_selector.currentIndexChanged.connect(self.render_selected_molecule)
        style_row.addWidget(self.style_selector)
        display_layout.addLayout(style_row)

        molecule_row = QHBoxLayout()
        molecule_row.addWidget(QLabel("Molecule:"))
        self.molecule_selector = QComboBox()
        self.molecule_selector.setStyleSheet(combo_style)
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

        display_group.setLayout(display_layout)
        options_layout.addWidget(display_group)

        # Buttons Section
        button_layout = QHBoxLayout()

        self.export_sdf_button = QPushButton("Export Selected")
        self.export_sdf_button.setMinimumHeight(30)
        self.export_sdf_button.clicked.connect(self.on_export_sdf)
        self.export_sdf_button.setEnabled(False)

        # Button for exporting lowest energy molecules
        self.export_lowest_energy_molecules_btn = QPushButton("Export 10 Lowest Energy")
        self.export_lowest_energy_molecules_btn.setMinimumHeight(30)
        self.export_lowest_energy_molecules_btn.clicked.connect(self.export_lowest_energy_molecules)
        self.export_lowest_energy_molecules_btn.setEnabled(False)

        button_layout.addWidget(self.export_sdf_button)
        button_layout.addWidget(self.export_lowest_energy_molecules_btn)
        options_layout.addLayout(button_layout)

        # Molecular Information Section
        self.info_group = QGroupBox("Molecular Information")
        self.info_layout = QFormLayout()

        self.formula_label = QLabel("N/A")  # all states N/A before file imported.
        self.weight_label = QLabel("N/A")
        self.atoms_label = QLabel("N/A")
        self.bonds_label = QLabel("N/A")
        self.charge_label = QLabel("N/A")
        self.mult_label = QLabel("N/A")
        self.energy_label = QLabel("N/A")
        self.source_file_label = QLabel("N/A")

        self.info_layout.addRow("Source File:", self.source_file_label)
        self.info_layout.addRow("Formula:", self.formula_label)
        self.info_layout.addRow("Weight:", self.weight_label)
        self.info_layout.addRow("Atoms:", self.atoms_label)
        self.info_layout.addRow("Bonds:", self.bonds_label)
        self.info_layout.addRow("Charge:", self.charge_label)
        self.info_layout.addRow("Multiplicity:", self.mult_label)
        self.info_layout.addRow("Energy:", self.energy_label)

        self.info_group.setLayout(self.info_layout)
        options_layout.addWidget(self.info_group)

        options_layout.addStretch()

        # 3D Viewer Section - This is the part that will be hidden/shown
        self.viewer_group = QGroupBox("3D Viewer")
        viewer_layout = QVBoxLayout()
        self.web_view = QWebEngineView()
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
        
        if not sdf_files:
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
                    print(f"Error loading {filename}: {str(file_error)}")
                    continue

            if not self.molecules:
                QMessageBox.warning(self, "Error", "No valid molecules found in the SDF files.")
                self.clear_molecules()
                return

            self.viewer_group.setVisible(True)
            
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

            self.export_sdf_button.setEnabled(True)
            self.export_lowest_energy_molecules_btn.setEnabled(True)

            if self.molecules:
                self.molecule_selector.setCurrentIndex(0)
                self.on_molecule_change(0)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load molecules from SDF files.\n\n{str(e)}")
            self.clear_molecules()
            self.viewer_group.setVisible(False)

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
        self.export_sdf_button.setEnabled(False)
        self.export_lowest_energy_molecules_btn.setEnabled(False)
        self.web_view.setHtml("")
        self.viewer_group.setVisible(False)
        self.formula_label.setText("N/A")
        self.weight_label.setText("N/A")
        self.atoms_label.setText("N/A")
        self.bonds_label.setText("N/A")
        self.charge_label.setText("N/A")
        self.mult_label.setText("N/A")
        self.energy_label.setText("N/A")
        self.source_file_label.setText("N/A")

    def on_molecule_change(self, index):
        if index < 0 or not self.molecules or index >= len(self.molecules):
            self.formula_label.setText("N/A")
            self.weight_label.setText("N/A")
            self.atoms_label.setText("N/A")
            self.bonds_label.setText("N/A")
            self.charge_label.setText("N/A")
            self.mult_label.setText("N/A")
            self.energy_label.setText("N/A")
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
            mult = self.getMult(mol)

            mol_name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"Molecule {original_idx + 1}"
            if not mol_name.strip():
                mol_name = f"Molecule {original_idx + 1}"

            self.formula_label.setText(formula)
            self.weight_label.setText(f"{mol_weight:.2f} g/mol")
            self.atoms_label.setText(str(num_atoms))
            self.bonds_label.setText(str(num_bonds))
            self.charge_label.setText(str(charge))
            self.mult_label.setText(str(mult))
            self.energy_label.setText("Calculating...")
            self.source_file_label.setText(source_file)

            QApplication.processEvents()

            # Calculate and display the energy
            try:
                # First check for 3D coordinates, if not present, embed them
                if mol.GetNumConformers() == 0:
                    AllChem.EmbedMolecule(mol, AllChem.ETKDG())

                # Perform a quick energy minimization
                ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
                ff.Minimize()
                energy = ff.CalcEnergy()

                # Update the displayed info with the energy
                self.energy_label.setText(f"{energy:.4f} kcal/mol")

            except Exception as energy_e:
                self.energy_label.setText("Error")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error getting molecule details: {str(e)}")
            self.formula_label.setText("Error")
            self.weight_label.setText("Error")
            self.atoms_label.setText("Error")
            self.bonds_label.setText("Error")
            self.charge_label.setText("Error")
            self.mult_label.setText("Error")
            self.energy_label.setText("Error")
            self.source_file_label.setText("Error")

        self.render_selected_molecule()

    def getMult(self, mol):
        """2S + 1, where S is the total spin, 1 unpaired electron = 1/2"""
        unpaired_electrons = 0
        for atom in mol.GetAtoms():
            unpaired_electrons += atom.GetNumRadicalElectrons()
        return int(unpaired_electrons / 2 + 1)

    def render_selected_molecule(self):
        if not self.molecules:
            return

        index = self.list_slider.value()
        if index < 0 or index >= len(self.molecules):
            return

        mol_data = self.molecules[index]
        mol = mol_data['mol']
        style = self.style_selector.currentText()

        if not mol.GetNumConformers():
            try:
                AllChem.EmbedMolecule(mol)
                AllChem.MMFFOptimizeMolecule(mol)
            except:
                QMessageBox.critical(self, "Error", "Could not generate 3D coordinates for this molecule.")
                return

        mol_block = Chem.MolToMolBlock(mol)

        viewer = py3Dmol.view(width='100%', height='100%')
        viewer.addModel(mol_block, 'mol')

        if style == 'Stick':
            viewer.setStyle({'stick': {}})
        elif style == 'Ball and Stick':
            viewer.setStyle({'stick': {}, 'sphere': {}})
        elif style == 'Surface':
            viewer.setStyle({'stick': {}})
            viewer.addSurface(py3Dmol.VDW, {'opacity': 0.8})

        viewer.zoomTo()
        html = viewer._make_html()
        self.web_view.setHtml(html)

    def on_export_sdf(self):
        if not self.molecules:
            QMessageBox.warning(self, "Error", "No molecules available. Please select SDF files in File Operations.")
            return

        current_index = self.molecule_selector.currentIndex()
        if current_index == -1:
            QMessageBox.warning(self, "Error", "No molecule selected. Please select a molecule to export.")
            return

        mol_data = self.molecules[current_index]
        mol = mol_data['mol']

        filename, _ = QFileDialog.getSaveFileName(self, "Save SDF File", "molecule.sdf", "SDF Files (*.sdf)")

        if filename:
            try:
                writer = SDWriter(filename)
                writer.write(mol)
                writer.close()

                export_message = f"SUCCESS: SDF file exported!\n\n"
                export_message += f"File saved to: {filename}\n"
                export_message += f"Molecule: {self.molecule_selector.currentText()}\n"
                export_message += f"Source: {mol_data['source_file']}"

                QMessageBox.information(self, "Export Successful", export_message)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to export SDF file.\n\n{str(e)}")

    def export_lowest_energy_molecules(self):
        if not self.molecules:
            QMessageBox.warning(self, "Error", "No molecules available. Please select SDF files.")
            return

        QMessageBox.information(self, "Calculation",
                                f"Calculating energies for {len(self.molecules)} molecules... This may take a moment.")
        QApplication.processEvents()

        molecules_with_energy = []
        for mol_data in self.molecules:
            mol = mol_data['mol']
            try:
                # Ensure the molecule has 3D coordinates
                if mol.GetNumConformers() == 0:
                    AllChem.EmbedMolecule(mol, AllChem.ETKDG())

                # Perform energy minimization and get the energy
                ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
                ff.Minimize()
                energy = ff.CalcEnergy()

                # Store the molecule data and its calculated energy
                molecules_with_energy.append((mol_data, energy))
            except Exception:
                # If energy calculation fails, skip this molecule
                continue

        if not molecules_with_energy:
            QMessageBox.critical(self, "Error", "No molecules could be processed for energy calculation.")
            return

        # Sort molecules by energy
        molecules_with_energy.sort(key=lambda x: x[1])

        # Get the top 10 lowest energy molecules, or fewer if not available
        lowest_10_molecules = [item[0]['mol'] for item in molecules_with_energy[:10]]

        filename, _ = QFileDialog.getSaveFileName(self, "Save Lowest Energy Molecules", "lowest_energy_molecules.sdf",
                                                  "SDF Files (*.sdf)")

        if filename:
            try:
                writer = SDWriter(filename)
                for mol in lowest_10_molecules:
                    writer.write(mol)
                writer.close()

                export_message = f"SUCCESS: Exported {len(lowest_10_molecules)} lowest energy molecules!\n\n"
                export_message += f"File saved to: {filename}\n"
                export_message += f"Selected from {len(self.molecules)} total molecules across {len(self.loaded_files)} files"
                QMessageBox.information(self, "Export Successful", export_message)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to export SDF file.\n\n{str(e)}")
        else:
            QMessageBox.information(self, "Export Cancelled", "Export cancelled by user.")