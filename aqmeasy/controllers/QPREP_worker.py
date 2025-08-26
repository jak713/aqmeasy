from PySide6.QtCore import Signal, QObject, Slot
from aqme.qprep import qprep
import os
import shutil


class QPrepWorker(QObject):
    resultReady = Signal(str)
    errorReady = Signal(str)
    finished = Signal()
    progressUpdate = Signal(str, int, int)  # message, current, total
    
    def __init__(self, params):
        super().__init__()
        self.params = params

    @Slot()
    def run(self):
        try:
            # Handle both single file (backward compatibility) and multiple files
            input_files = self.params.get('input_files', [])
            if not input_files and self.params.get('input_file'):
                input_files = [self.params['input_file']]
            
            if not input_files:
                self.errorReady.emit("No input files provided.")
                return

            # Create testfiles directory
            testfiles_dir = os.path.join(os.path.dirname(__file__), 'testfiles')
            if not os.path.exists(testfiles_dir):
                os.makedirs(testfiles_dir)

            # Copy all files to testfiles directory and track them
            processed_files = []
            copied_files = []
            
            for input_path in input_files:
                if not os.path.abspath(input_path).startswith(os.path.abspath(testfiles_dir)):
                    basename = os.path.basename(input_path)
                    dest_path = os.path.join(testfiles_dir, basename)
                    
                    # Handle duplicate filenames by adding a counter
                    counter = 1
                    original_dest = dest_path
                    while os.path.exists(dest_path):
                        name, ext = os.path.splitext(original_dest)
                        dest_path = f"{name}_{counter}{ext}"
                        counter += 1
                    
                    shutil.copy2(input_path, dest_path)
                    processed_files.append(dest_path)
                    copied_files.append(dest_path)
                else:
                    processed_files.append(input_path)

            print("***********************")
            print(f"Debug: Running qprep with {len(processed_files)} files:", processed_files)
            print("Debug: Parameters:", self.params)
            print("***********************")

            # Process files in batch
            total_files = len(processed_files)
            
            for i, file_path in enumerate(processed_files):
                self.progressUpdate.emit(f"Processing file {i+1} of {total_files}: {os.path.basename(file_path)}", 
                                       i+1, total_files)
                
                # Create individual parameters for each file
                file_params = self.params.copy()
                file_params['input_file'] = file_path
                
                try:
                    qprep(files=file_path,
                          program=self.params['software'],
                          qm_input=self.params["keywords"] + ' ' + self.params['functional'] + " " + 
                                  self.params['basis'] + " " + self.params['solvent_block'],
                          mem=f"{self.params['mem']}GB",
                          nprocs=self.params['nprocs'],
                    )
                except Exception as file_error:
                    print(f"Error processing {file_path}: {str(file_error)}")
                    # Continue with other files even if one fails
                    continue

            # Move all generated files to user-selected output directory if specified
            output_dir = self.params.get('output_dir', '').strip()
            qcalc_dir = os.path.join(os.path.dirname(__file__), 'QCALC')
            
            moved_files_count = 0
            if output_dir and os.path.isdir(output_dir) and os.path.abspath(output_dir) != os.path.abspath(qcalc_dir):
                if os.path.exists(qcalc_dir):
                    for fname in os.listdir(qcalc_dir):
                        src = os.path.join(qcalc_dir, fname)
                        dst = os.path.join(output_dir, fname)
                        if os.path.isfile(src):
                            try:
                                # Handle duplicate filenames in output directory
                                counter = 1
                                original_dst = dst
                                while os.path.exists(dst):
                                    name, ext = os.path.splitext(original_dst)
                                    dst = f"{name}_{counter}{ext}"
                                    counter += 1
                                
                                shutil.move(src, dst)
                                moved_files_count += 1
                            except Exception as move_error:
                                print(f"Error moving file {src}: {str(move_error)}")
                    
                    # Clean up empty QCALC directory
                    try:
                        if os.path.exists(qcalc_dir) and not os.listdir(qcalc_dir):
                            shutil.rmtree(qcalc_dir, ignore_errors=True)
                    except Exception:
                        pass

            # Clean up copied files from testfiles directory
            for copied_file in copied_files:
                try:
                    if os.path.exists(copied_file):
                        os.remove(copied_file)
                except Exception:
                    pass

            self.finished.emit()
            
            result_message = f"Successfully processed {total_files} input file(s)."
            if moved_files_count > 0:
                result_message += f"\n{moved_files_count} output files moved to specified directory."
            
            self.resultReady.emit(result_message)
            
        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            self.errorReady.emit(f"Error: {str(e)}\n\nTraceback:\n{tb}")
            self.finished.emit()