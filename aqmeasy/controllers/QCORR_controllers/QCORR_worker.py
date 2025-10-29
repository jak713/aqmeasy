import glob
from importlib.metadata import files
from PySide6.QtCore import QObject, Slot, QRunnable, QThreadPool, Signal
from aqme.qcorr import qcorr
from aqmeasy.models.QCORR_model.QCORR_parammodel import default_values
import os
import io, contextlib, traceback

class QCORRWorker(QObject):
    result = Signal(str)
    error = Signal(str)
    finished = Signal()

    def __init__(self, parent):
        super().__init__()
        self.param_model = parent.parameter_model
        self.file_model = parent.file_model
        self.threadpool = QThreadPool()

    def run_qcorr(self):
        worker = Worker(self,self.param_model, self.file_model)
        self.threadpool.start(worker)


class Worker(QRunnable):
    def __init__(self, parent,param_model, file_model):
        super().__init__()
        self.param_model = param_model
        self.file_model = file_model
        self.parent = parent

    @Slot()
    def run(self):
        """
        Works as follows:
        1. Collects non-default parameters from ParamModel
        2. Collects files and working directory from FileModel
        3. Sets up the files in the working directory (for qcorr to simply accept *.log / *.out files)
        4. Runs the qcorr command with the collected parameters and files
        """
        params = self.collect_qcorr_params()
        print("QCORR parameters to be used:", params)
        w_dir_main  = self.set_up_files_in_wdir()
        print("QCORR working directory:", w_dir_main)

        # w_dir_main = w_dir_main.strip('/') # for some reason qcorr dislikes leading slash?
        # files = w_dir_main+'/*.*'
        # above does not work, hence we SET CURRENT WORKING DIRECTORY to w_dir_main and give files as '*.*'
        # print("QCORR files to be processed:", files)
        files = '*.*'
        os.chdir(w_dir_main)
        # check files in files directory to see whether they exist:
        if not self.check_files_exist(files):
            print("No files found in the specified working directory. Aborting QCORR run.")
            return
            
        buf_out = io.StringIO()
        buf_err = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf_out), contextlib.redirect_stderr(buf_err):
                qcorr(files='*.*', **params)
        except Exception:
            buf_err.write(traceback.format_exc())

        stdout_text = buf_out.getvalue()
        stderr_text = buf_err.getvalue()

        result = ""
        if stdout_text:
            result += stdout_text + "\n"
        if stderr_text:
            result += "ERROR STREAM:\n" + stderr_text

        self.qcorr_result = result

        # emit if a Signal was provided/attached
        self.parent.result.emit(result)
        # After running qcorr, we can set the working directory back
        self.file_model.w_dir_main = w_dir_main
        # emit w dir change signal
        self.file_model.wdirChanged.emit(w_dir_main)

    def check_files_exist(self, files_pattern):
        """
        Checks if there are files matching the given pattern.
        """
        matched_files = glob.glob(files_pattern)
        print("Matched files:", matched_files)
        return len(matched_files) > 0
    
    def collect_qcorr_params(self):
        """
        Collects qcorr parameters from ParamModel as dict, compares them to default_values, if different stores them as attributes of self.
        """
        params = getattr(self, "param_model", None)
        if params is None:
            raise AttributeError("self.param_model is missing")
        params = params.as_dict()

        changed = {}
        for key, value in params.items():
            if key in default_values and default_values[key] != value:
                # if hasattr(self, key):
                #     pass
                # setattr(self, key, value)
                changed[key] = value
        return changed

    def collect_qcorr_files(self):
        """
        Collects the list of files from FileModel.
        """
        files = getattr(self, "file_model", None)
        if files is None:
            raise AttributeError("self.file_model is missing")
        return files.files
    
    def collect_qcorr_wdir(self):
        """
        Collects the working directory from FileModel. If not set, use directory of first file.
        """
        files = getattr(self, "file_model", None)
        if files is None:
            raise AttributeError("self.file_model is missing")
        return files.w_dir_main or os.path.dirname(files.files[0])
    
    def set_up_files_in_wdir(self):
        """
        Sets (moves over) the files in the working directory to be processed by QCORR.

        need to see whether we want to move files or copy...
        """
        wdir = self.collect_qcorr_wdir() + '/QCORR'
        files = self.collect_qcorr_files()
        if os.path.exists(wdir):
            # if exists, create new dir with incremented number
            i = 1
            new_wdir = f"{wdir}_{i}"
            while os.path.exists(new_wdir):
                i += 1
                new_wdir = f"{wdir}_{i}"
            os.makedirs(new_wdir)
            wdir = new_wdir
        else:
            os.makedirs(wdir)
        # Move files to wdir (this may have to be changed, feels hacky at the moment)
        for file in files:
            filename = os.path.basename(file)
            dest = os.path.join(wdir, filename)
            if not os.path.exists(dest):
                os.rename(file, dest)
        return wdir