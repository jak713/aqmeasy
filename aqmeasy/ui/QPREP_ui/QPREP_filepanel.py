import os
import sys
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QPushButton, QTextEdit, QFileDialog, QListWidget, QListWidgetItem, QSizePolicy, QApplication, QStyle
)
from PySide6.QtCore import QThread, Signal
from PySide6.QtGui import QIcon
from aqmeasy.controllers.QPREP_worker import QPrepWorker
from aqmeasy.ui.stylesheets import stylesheets
from aqmeasy.ui.icons import Icons
FILE_FILTERS = [
    "Structured Data Files (*.sdf)",
    "XYZ (*.xyz)",
    "Gaussian Input Files (*.com *.gjf)",
]


class FilePanel(QWidget):
    """Panel to manage input file selection, output directory, and generate input scripts."""

    # Signal to notify when new files are selected
    filesSelected = Signal(list)
    fileSelected = Signal(str)  # Keep for backward compatibility

    # New signal for when a single file is selected from the list
    singleFileSelected = Signal(str)

    def __init__(self, model=None, parameter_panel=None, molecular_viewer=None):
        super().__init__()
        self.setAcceptDrops(True)
        self.model = model
        self.parameter_panel = parameter_panel
        self.molecular_viewer = molecular_viewer
        self.selected_files = []
        self.setup_ui()

        if self.molecular_viewer:
            self.filesSelected.connect(self.molecular_viewer.load_molecules_from_files)
            # Keep backward compatibility for single file
            self.fileSelected.connect(self.molecular_viewer.load_molecules_from_file)
            # Connect the new signal to the molecular viewer
            self.singleFileSelected.connect(self.molecular_viewer.load_molecules_from_file)


    def setup_ui(self):
        layout = QVBoxLayout()

        # Input File Section
        input_group = QGroupBox("Input Files")
        input_layout = QVBoxLayout()
        
        # File selection buttons
        button_row = QHBoxLayout()
        self.browse_multiple_btn = QPushButton()
        self.browse_multiple_btn.setIcon(QIcon(Icons.file_open))
        self.browse_multiple_btn.clicked.connect(self.get_multiple_filenames)
        self.clear_files_btn = QPushButton()
        self.clear_files_btn.setIcon(QIcon(Icons.trash))
        self.clear_files_btn.clicked.connect(self.clear_files)
        
        button_row.addWidget(self.browse_multiple_btn)
        button_row.addWidget(self.clear_files_btn)
        input_layout.addLayout(button_row)
        
        # File list display
        self.file_list = QListWidget()
        self.file_list.setMaximumHeight(120)
        # Connect the itemClicked signal to the new method
        self.file_list.itemClicked.connect(self.on_file_list_item_clicked)
        input_layout.addWidget(self.file_list)
        
        input_group.setLayout(input_layout)

        # Output Directory Section
        output_group = QGroupBox("Output")
        output_layout = QVBoxLayout()
        output_row = QHBoxLayout()
        output_row.addWidget(QLabel("Output Directory"))
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setPlaceholderText("Select output directory")
        output_row.addWidget(self.output_dir_edit)
        self.browse_output_btn = QPushButton()
        self.browse_output_btn.setIcon(QIcon(Icons.folder_open))
        self.browse_output_btn.clicked.connect(self.get_output_directory)
        output_row.addWidget(self.browse_output_btn)
        output_layout.addLayout(output_row)
        output_group.setLayout(output_layout)

        # Buttons
        self.make_btn = QPushButton("Run QPREP")
        self.make_btn.setStyleSheet(stylesheets.RunButton)
        self.make_btn.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum)
        self.make_btn.setIcon(QApplication.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay))
        self.make_btn.setMinimumHeight(40)
        self.make_btn.clicked.connect(self.run_qprep)
        self.preview_btn = QPushButton("Preview")
        self.preview_btn.setMinimumHeight(30)
        self.preview_btn.clicked.connect(self.preview_input)
        button_layout = QHBoxLayout()
        button_layout.addWidget(self.preview_btn)
        button_layout.addWidget(self.make_btn)
        self.slurm_btn = QPushButton("SLURM")
        self.slurm_btn.setMinimumHeight(30)
        self.slurm_btn.clicked.connect(self.generate_slurm_script)
        button_layout.addWidget(self.slurm_btn)

        # Status Display
        status_group = QGroupBox("Status / Preview")
        status_layout = QVBoxLayout()

        self.status_text = QTextEdit()
        self.status_text.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.status_text.setPlaceholderText("Ready to create input files")
        self.status_text.setReadOnly(True)
        self.status_text.setMinimumHeight(300)
        status_layout.addWidget(self.status_text)
        status_group.setLayout(status_layout)

        # Final Layout
        layout.addWidget(input_group)
        layout.addWidget(output_group)
        layout.addLayout(button_layout)
        layout.addWidget(status_group)
        layout.addStretch()
        self.setLayout(layout)


    def setup_fallback_timer(self):
        from PySide6.QtCore import QTimer
        self.last_software = None
        self.timer = QTimer()
        self.timer.timeout.connect(self.check_software_change)
        self.timer.start(500)

    def check_software_change(self):
        if self.model:
            current_software = self.model.software().lower()
            if current_software != self.last_software:
                self.last_software = current_software
                self.update_button_color()


    def clear_files(self):
        """Clear all selected files"""
        self.selected_files.clear()
        self.file_list.clear()
        self.status_text.setPlainText("All files cleared. Ready to select new files.")
        self.filesSelected.emit([])

    def get_files_from_csearch(self, file_list):
        """Receive files from CSEARCH"""
        self.selected_files = file_list
        self.update_file_list()
        self.display_files_info()
        self.filesSelected.emit(self.selected_files)

    def get_single_filename(self):
        """Select a single file"""
        initial_filter = FILE_FILTERS[0]
        filters = ';;'.join(FILE_FILTERS)
        filename, selected_filter = QFileDialog.getOpenFileName(
            self,
            caption="Select Single File",
            filter=filters,
            selectedFilter=initial_filter
        )
        if filename:
            self.selected_files = [filename]  # Replace with single file
            self.update_file_list()
            self.display_files_info()
            # Emit both signals for compatibility
            self.fileSelected.emit(filename)
            self.filesSelected.emit(self.selected_files)

    def get_multiple_filenames(self):
        """Select multiple files"""
        initial_filter = FILE_FILTERS[0]
        filters = ';;'.join(FILE_FILTERS)
        filenames, selected_filter = QFileDialog.getOpenFileNames(
            self,
            caption="Select Multiple Files",
            filter=filters,
            selectedFilter=initial_filter
        )
        if filenames:
            self.selected_files = filenames
            self.update_file_list()
            self.display_files_info()
            self.filesSelected.emit(self.selected_files)
            # For backward compatibility and to display something
            if self.selected_files:
                self.fileSelected.emit(self.selected_files[0])
                # Automatically select the first item and load it
                first_item = self.file_list.item(0)
                self.file_list.setCurrentItem(first_item)
                self.on_file_list_item_clicked(first_item)

    def update_file_list(self):
        """Update the file list widget"""
        self.file_list.clear()
        for filename in self.selected_files:
            item = QListWidgetItem(os.path.basename(filename))
            item.setToolTip(filename)  # Show full path on hover
            self.file_list.addItem(item)
    
    def on_file_list_item_clicked(self, item):
        """Method to handle clicks on the file list widget."""
        index = self.file_list.row(item)
        if 0 <= index < len(self.selected_files):
            file_path = self.selected_files[index]
            self.singleFileSelected.emit(file_path)


    def display_files_info(self):
        """Display information about selected files"""
        if not self.selected_files:
            self.status_text.setPlainText("No files selected.")
            return

        try:
            status_text = f"Selected {len(self.selected_files)} file(s):\n"
            total_size = 0
            sdf_count = 0
            
            for i, filename in enumerate(self.selected_files):
                file_basename = os.path.basename(filename)
                file_size = os.path.getsize(filename)
                file_ext = os.path.splitext(filename)[1].lower()
                total_size += file_size
                
                status_text += f"{i+1}. {file_basename} ({file_size:,} bytes, {file_ext})\n"
                
                if file_ext == '.sdf':
                    sdf_count += 1
                    try:
                        from rdkit import Chem
                        supplier = Chem.SDMolSupplier(filename)
                        mol_count = len([mol for mol in supplier if mol is not None])
                        status_text += f"   → {mol_count} molecules\n"
                    except Exception:
                        status_text += "   → Could not read molecule count\n"

            status_text += f"\nTotal size: {total_size:,} bytes\n"
            if sdf_count > 0:
                status_text += f"SDF files: {sdf_count}\n"
            
            status_text += "-" * 50 + "\n"
            status_text += "Click 'Preview Input' to see what the generated input files will look like,\n"
            status_text += "or click 'Make Input Files' to create the actual files for all selected files."

            if sdf_count > 0:
                status_text += f"\n\nFor SDF files: View molecules in the Molecular Analysis panel on the right."

            self.status_text.setPlainText(status_text)
        except Exception as e:
            error_msg = f"Error reading files: {str(e)}"
            self.status_text.setPlainText(error_msg)

    def preview_input(self):
        params = self.collect_input_params()
        if not params['input_files']:
            self.status_text.setPlainText("Error: No input files selected. Please select files first.")
            self.preview_btn.setEnabled(True)
            return
            
        try:
            preview_content = ""
            for i, input_file in enumerate(params['input_files']):
                file_params = params.copy()
                file_params['input_file'] = input_file
                file_params['molecule_name'] = os.path.splitext(os.path.basename(input_file))[0]
                
                preview_content += f"=== FILE {i+1}: {os.path.basename(input_file)} ===\n\n"
                
                if params['software'].lower() == 'orca':
                    file_preview = self.generate_orca_preview(file_params)
                elif params['software'].lower() == 'gaussian':
                    file_preview = self.generate_gaussian_preview(file_params)
                else:
                    file_preview = "Error: Unsupported software type."
                
                preview_content += file_preview + "\n\n"
                
                if i < len(params['input_files']) - 1:  # Add separator between files
                    preview_content += "=" * 70 + "\n\n"
            
            preview_content += f"{'=' * 70}\n\nThis is a preview of the input files that will be generated.\n"
            preview_content += f"Total files to process: {len(params['input_files'])}\n"
            preview_content += "Click 'Make Input Files' to create the actual files."
            
            self.status_text.setPlainText(preview_content)
        except Exception as e:
            self.status_text.setPlainText(f"Error generating preview: {str(e)}")
        self.preview_btn.setEnabled(True)

################################################################################
    # Drag and drop events
    def dragEnterEvent(self, event):
        urls = event.mimeData().urls()
        if not urls:
            event.ignore()
            return
        for url in urls:
            if not url.isLocalFile():
                event.ignore()
                return
        event.acceptProposedAction()

    def dropEvent(self, event):
        urls = event.mimeData().urls()
        dirs = [url.toLocalFile() for url in urls if url.isLocalFile() and os.path.isdir(url.toLocalFile())]
        if not dirs:
            event.ignore()
            return
        folder_path = dirs[0] if len(dirs) == 1 else dirs
        event.acceptProposedAction()

################################################################################

    def generate_orca_preview(self, params):
        calc_line = f"! {params['functional']} {params['basis']}"
        if params['keywords']:
            calc_line += f" {params['keywords']}"
        if params['solvent_block']:
            calc_line += f" {params['solvent_block']}"
        preview = f"{calc_line}\n\n"
        preview += f"%pal nprocs {params['nprocs']} end\n"
        preview += f"%maxcore {params['mem']}000\n\n"
        preview += f"* xyz {params['charge']} {params['multiplicity']}\n"
        preview += "[Molecular coordinates will be extracted from input file]\n*\n"
        return preview

    def generate_gaussian_preview(self, params):
        preview = f"%nprocshared={params['nprocs']}\n"
        preview += f"%mem={params['mem']}GB\n"
        route = f"#p {params['functional']}/{params['basis']}"
        if params['keywords']:
            route += f" {params['keywords']}"
        if params['solvent_block']:
            route += f" {params['solvent_block']}"
        preview += f"{route}\n\n"
        preview += f"{params['molecule_name']}\n\n"
        preview += f"{params['charge']} {params['multiplicity']}\n"
        preview += "[Molecular coordinates will be extracted from input file]\n\n"
        return preview

    def run_qprep(self):
        params = self.collect_input_params()
        if not params['input_files']:
            self.status_text.setPlainText("Error: No input files selected. Please select files.")
            return
            
        file_count = len(params['input_files'])
        self.status_text.setPlainText(
            f"Processing... Generating input files for {file_count} file(s). \n\n"
            f"This may take a moment depending on the size and number of your input files.")
        
        self.thread = QThread()
        self.worker = QPrepWorker(params)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.worker.resultReady.connect(self.on_qprep_success)
        self.worker.errorReady.connect(self.on_qprep_error)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.start()
        self.make_btn.setEnabled(False)
        self.preview_btn.setEnabled(False)
        self.thread.finished.connect(lambda: self.make_btn.setEnabled(True))
        self.thread.finished.connect(lambda: self.preview_btn.setEnabled(True))

    def on_qprep_success(self, message):
        params = self.collect_input_params()
        output_dir = params.get('output_dir', 'Current directory')
        file_count = len(params['input_files'])
        
        success_message = f"SUCCESS: Input files generated successfully!\n\n"
        success_message += f"Files processed: {file_count}\n"
        success_message += f"Output location: {output_dir}\n"
        success_message += f"Software: {params['software']}\n"
        success_message += f"Method: {params['functional']}/{params['basis']}\n"
        success_message += f"Memory: {params['mem']}GB\n"
        success_message += f"Processors: {params['nprocs']}\n"
        if params['solvent_model'] != 'None':
            success_message += f"🧪 Solvent: {params['solvent']} ({params['solvent_model']})\n"
        success_message += f"\n{message}"
        self.status_text.setPlainText(success_message)

    def on_qprep_error(self, error_message):
        self.status_text.setPlainText(f"ERROR: Failed to generate input files.\n\n{error_message}")

    def collect_input_params(self):
        params = {}
        if self.model:
            params['software'] = self.model.software().lower()
            params['functional'] = self.model.functional()
            params['basis'] = self.model.basisSet()
            params['nprocs'] = self.model.nprocs()
            params['mem'] = self.model.mem()
            
        # Use selected_files list instead of single file
        params['input_files'] = self.selected_files
        # Keep input_file for backward compatibility (first file)
        params['input_file'] = self.selected_files[0] if self.selected_files else ''
        
        output_dir = self.output_dir_edit.text().strip()
        if output_dir != '':
            params['output_dir'] = output_dir
        else:
            # If no output directory specified, use directory of the first input file
            params['output_dir'] = os.path.dirname(params['input_files'][0]) if params['input_files'] else ''

        if self.parameter_panel:
            try:
                params['charge'] = int(self.parameter_panel.get_current_charge())
            except Exception:
                params['charge'] = 0
            try:
                params['multiplicity'] = int(self.parameter_panel.get_current_multiplicity())
            except Exception:
                params['multiplicity'] = 1
        else:
            params['charge'] = 0
            params['multiplicity'] = 1

        if self.parameter_panel:
            params['job_type'] = self.parameter_panel.get_current_job_type()
            params['solvent_model'] = self.parameter_panel.get_current_solvent_model()
            params['solvent'] = self.parameter_panel.get_current_solvent()
        else:
            params['job_type'] = 'Geometry Optimization'
            params['solvent_model'] = 'None'
            params['solvent'] = 'None'
            
        if params['solvent_model'] != 'None' and params['solvent'] != 'None':
            solvent_blocks = {
                'orca': {
                    'CPCM': f"CPCM({params['solvent']})",
                    'SMD': f"SMD({params['solvent']})",
                    'COSMORS': f"COSMORS({params['solvent']})",
                    'DRACO': f"CPCM({params['solvent']}) DRACO",
                },
                'gaussian': {
                    'PCM': f"SCRF=(PCM,Solvent={params['solvent']})",
                    'SMD': f"SCRF=(SMD,Solvent={params['solvent']})",
                    'IEFPCM': f"SCRF=(IEFPCM,Solvent={params['solvent']})",
                    'CPCM': f"SCRF=(CPCM,Solvent={params['solvent']})",
                }
            }
            params["solvent_block"] = solvent_blocks.get(params['software'], {}).get(params['solvent_model'], '')
        else:
            params["solvent_block"] = ''
            
        job_types = {
            'Single Point': '',
            'Geometry Optimization': 'opt',
            'Frequency': 'freq',
            'Opt+Freq': 'opt freq',
        }
        params["keywords"] = job_types.get(params['job_type'], "")
        
        # Set molecule_name to first file for general use
        first_file = params['input_files'][0] if params['input_files'] else 'input'
        params['molecule_name'] = os.path.splitext(os.path.basename(first_file))[0]
        
        return params

    def get_output_directory(self):
        directory = QFileDialog.getExistingDirectory(
            self,
            "Select Output Directory",
            "",
            QFileDialog.ShowDirsOnly
        )
        if directory:
            self.output_dir_edit.setText(directory)
            self.update_output_status(directory)

    def update_output_status(self, directory):
        try:
            dir_basename = os.path.basename(directory) or directory
            current_text = self.status_text.toPlainText()
            output_info = f"\nOutput directory: {dir_basename}\n"
            output_info += f"Full path: {directory}\n"
            output_info += "-" * 50 + "\n"
            output_info += "Ready to generate input files!"
            if current_text and not current_text.startswith("Ready to create"):
                updated_text = current_text + "\n" + output_info
            else:
                updated_text = output_info
            self.status_text.setPlainText(updated_text)
        except Exception as e:
            error_msg = f"Error with output directory: {directory}\nError: {str(e)}"
            self.status_text.append(error_msg)

    def generate_slurm_script(self):
        params = self.collect_input_params()
        input_files = params.get('input_files', [])
        nprocs = params.get('nprocs', 1)
        software = params.get('software', 'orca').lower()
        
        if not input_files:
            self.status_text.setPlainText("Error: No input files selected for SLURM script generation.")
            return

        output_dir = params.get('output_dir', '').strip()
        if output_dir:
            script_dir = output_dir
        else:
            script_dir = os.path.dirname(input_files[0])

        if software == 'gaussian':
            script_filename = "submitgaussian_batch.txt"
        else:
            script_filename = "submitorca_batch.txt"
        script_path = os.path.join(script_dir, script_filename)

        # Generate batch script for multiple files
        if software == 'gaussian':
            script_content = f'''#!/bin/bash --login
#SBATCH -p multicore # (or --partition) Single-node multicore
#SBATCH -n {nprocs} # (or --ntasks=) Number of cores (2--40)
#SBATCH -t 4-0
#SBATCH --array=1-{len(input_files)}

# Load g16 for the CPU type our job is running on
module load gaussian/g16c01_em64t_detectcpu

## Set up scratch dir (please do this!)
export GAUSS_SCRDIR=/scratch/$USER/gau_temp_${{SLURM_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}
mkdir -p $GAUSS_SCRDIR

## Say how much memory to use (4GB per core)
export GAUSS_MDEF=$((SLURM_NTASKS*4))GB

## Inform Gaussian how many cores to use
export GAUSS_PDEF=$SLURM_NTASKS

# Define input files array
input_files=(
'''
            for filename in input_files:
                base_name = os.path.splitext(os.path.basename(filename))[0]
                inp_basename = base_name + '.inp'
                script_content += f'    "{inp_basename}"\n'
            
            script_content += ''')

# Get current file based on array index
current_input="${input_files[$SLURM_ARRAY_TASK_ID-1]}"
current_output="${current_input%.*}.out"

$g16root/g16/g16 < "$current_input" > "$current_output"
'''
        else:
            script_content = f'''#!/bin/bash --login

#SBATCH -p multicore    # (or --partition=) Submit to the AMD Genoa nodes
#SBATCH -n {nprocs}            # (or --ntasks=) Number of cores (2--168). Must match the number in your ORCA input file!
#SBATCH -t 4-0          # Wallclock timelimit (4-0 is 4 days, max permitted is 7-0)
#SBATCH --array=1-{len(input_files)}

module purge
module load apps/binapps/orca/6.0.1-avx2

# Define input files array
input_files=(
'''
            for filename in input_files:
                base_name = os.path.splitext(os.path.basename(filename))[0]
                inp_basename = base_name + '.inp'
                script_content += f'    "{inp_basename}"\n'
            
            script_content += ''')

# Get current file based on array index
current_input="${input_files[$SLURM_ARRAY_TASK_ID-1]}"

$ORCA_HOME/orca "$current_input" > results.${{SLURM_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}.txt
'''

        try:
            with open(script_path, 'w') as f:
                f.write(script_content)
            self.status_text.setPlainText(
                f"SLURM batch script generated: {script_path}\n"
                f"This script will process {len(input_files)} files using job arrays.\n"
                f"Submit with: sbatch {script_filename}")
        except Exception as e:
            self.status_text.setPlainText(f"Error writing SLURM script: {str(e)}")

    # Keep original methods for backward compatibility
    def get_filename(self):
        """Backward compatibility method"""
        self.get_single_filename()

    def display_file_info(self, filename):
        """Backward compatibility method"""
        self.selected_files = [filename]
        self.update_file_list()
        self.display_files_info()