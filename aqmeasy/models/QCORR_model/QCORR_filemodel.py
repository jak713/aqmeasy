import os
import glob
from PySide6.QtCore import QObject, Signal

class FileModel(QObject):
    """Signals for the QCORR file model"""
    
    filesChanged = Signal(list)
    wdirChanged = Signal(str)
    currentlySelectedFileChanged = Signal(str)
    isSelectableChanged = Signal(bool)
    fileStatusChanged = Signal(dict)

    def __init__(self):
        super().__init__()
        self.files = []
        self.w_dir_main = ""
        self.currently_selected_file = ""
        self.isSelectable = True
        self.file_statuses = {} 
        self.original_to_new_path = {}

    def __get__files__(self):
        return self.files

    def __get__w_dir_main__(self):
        return self.w_dir_main

    def __get__currently_selected_file__(self):
        return self.currently_selected_file
    
    def set_file_status(self, filepath, status):
        """
        Set status for a file.
        
        Args:
            filepath: Full path to file
            status: 'success', 'failed', or 'pending'
        """
        self.file_statuses[filepath] = status
        self.fileStatusChanged.emit(self.file_statuses)

    def update_files_after_qcorr(self, output_dir):
        """
        After QCORR runs, discover where files were reorganized to.
        
        Args:
            output_dir: Base output directory (contains QCORR_X folders)
        """
        import glob
        
        # Find the latest QCORR_X directory
        qcorr_dirs = glob.glob(os.path.join(output_dir, 'QCORR*'))
        if not qcorr_dirs:
            return
        
        # Sort to get latest (QCORR_1, QCORR_2, etc.)
        latest_qcorr = sorted(qcorr_dirs)[-1]
        
        # Map old filenames to new paths
        new_files = []
        self.file_statuses.clear()
        
        # Scan success directory recursively (e.g., success/SP_calcs)
        success_dir = os.path.join(latest_qcorr, 'success')
        if os.path.isdir(success_dir):
            for root, dirs, files in os.walk(success_dir):
                if 'json_files' in dirs:
                    dirs.remove('json_files')

                for filename in files:
                    filepath = os.path.join(root, filename)
                    if self._is_output_file(filepath):
                        new_files.append(filepath)
                        self.file_statuses[filepath] = 'success'
        
        # Scan failed directory and its subdirectories
        failed_dir = os.path.join(latest_qcorr, 'failed')
        if os.path.isdir(failed_dir):
            # Walk through error subdirectories
            for root, dirs, files in os.walk(failed_dir):
                # ADDED: Skip json_files directories entirely
                if 'json_files' in dirs:
                    dirs.remove('json_files')  # Don't descend into json_files
                
                for filename in files:
                    filepath = os.path.join(root, filename)
                    # UPDATED: Only track output files (not JSON, not new inputs)
                    if self._is_output_file(filepath):
                        new_files.append(filepath)
                        # Derive error type from subdirectory name
                        rel_path = os.path.relpath(root, failed_dir)
                        error_type = rel_path if rel_path != '.' else 'failed'
                        self.file_statuses[filepath] = f'failed_{error_type}'
        
        # Update model
        self.files = new_files
        self.w_dir_main = latest_qcorr
        
        # Emit signals
        self.filesChanged.emit(self.files)
        self.wdirChanged.emit(self.w_dir_main)
        self.fileStatusChanged.emit(self.file_statuses)
    

    def _is_output_file(self, filepath):
        """
        Determine if a file is a QM output file (not JSON, not new input).
        
        Returns:
            True if file is a QM output (.log, .out, etc.)
            False if file is JSON, input file, or other
        """
        filename = os.path.basename(filepath)
        
        # Skip JSON files
        if filename.endswith('.json'):
            return False
        
        # Skip new input files created by QCORR
        if filename.endswith(('.com', '.inp', '.gjf', '.in')):
            return False
        
        # Skip hidden files
        if filename.startswith('.'):
            return False
        
        # Common QM output extensions
        output_extensions = ('.log', '.out', '.output', '.xyz', '.sdf', '.mol')
        
        # Accept if it has a known output extension
        if filename.endswith(output_extensions):
            return True
        
        # If no extension match, check if it's NOT in json_files directory
        # (belt and suspenders approach)
        if 'json_files' in filepath:
            return False
        
        # Default: accept it (catches non-standard output files)
        return True
        
    def get_file_status(self, filepath):
        """Get status for a file."""
        return self.file_statuses.get(filepath, 'unknown')

    def as_dict(self):
        """Returns the current state of choices as a dict, and is once again updated upon changing, so that correct user selection is created in the file."""
        return {
            'files': self.files,
            'w_dir_main': self.w_dir_main,
            'file_statuses': self.file_statuses
        }