import os
import sys
import pytest
from pytestqt import qtbot
from PySide6.QtCore import Qt

from aqmeasy.models.CSEARCH_model.CSEARCH_command import general_command_default
from aqmeasy.ui.CSEARCH_ui.CSEARCH import CSEARCH
from aqmeasy.controllers.CSEARCH_controller import Worker

test_smiles = ["CCO", "C1=CC=CC=C1", "C(C(=O)O)N", "CC(=O)OC1=CC=CC=C1C(=O)O", "C1CCCCC1", "C1=CC=CN=C1", "CCN(CC)CC", "C(CCl)Br", "CC(C)O", "C1=CC=C(C=C1)O"]

test_code_names = ["ethanol", "benzene", "glycine", "aspirin", "cyclohexane", "pyridine", "triethylamine", "bromochloromethane", "isopropanol", "phenol"]

@pytest.fixture
def csearch_widget(qtbot):
    widget = CSEARCH()
    widget.hide()  # Prevent the widget from showing during tests
    qtbot.addWidget(widget)
    return widget

# To test some components directly:
@pytest.fixture
def main_widget(csearch_widget):
    return csearch_widget.main_widget

@pytest.fixture
def csv_table(csearch_widget, qtbot, main_widget):
    """This has to be initialised before using."""
    main_widget.control.show_csv()
    csv = csearch_widget.main_widget.control.csv
    csv.hide()  # Hide the CSV table window during tests
    return csv

@pytest.fixture
def control(csearch_widget):
    return csearch_widget.main_widget.control

@pytest.fixture
def csv_model(csearch_widget):
    return csearch_widget.main_widget.control.model

class TestUIElements:
    """
    For dynamic creation of widgets:
    1. Initialisation from the CSEARCH class. (model, worker, main_widget, stylesheet, window title)
    2. Show All button opens csv table window.
    """

    def test_csearch_initialisation(self, csearch_widget):
        assert csearch_widget.model is not None
        assert csearch_widget.worker is not None
        assert csearch_widget.main_widget is not None
        assert csearch_widget.styleSheet() is not None
        assert csearch_widget.windowTitle() == "CSEARCH"

    def test_csearch_ui_show_all_button(self, main_widget, csv_table):
        """Pressing the button should trigger self.show_all_button.clicked.connect(self.control.show_csv)"""
        assert main_widget.show_all_button is not None
        assert csv_table is not None

class TestSyncToCSVModel:
    """
    Testing signals/slots between UI and model. 
    The CSV model is a dictionary model holding the csv data in the form of key (column name) : list (column data).

    Keys: "SMILES", "code_name", "charge", "multiplicity", "constraints_atoms", "constraints_dist", "constraints_angle", "constraints_dihedral", "complex_type", "geom"

    1. SMILES input fields updates the model's SMILES entry.
    2. Model's SMILES entry updates the SMILES input field. (e.g., when navigating between entries in the CSV table)
    3. 
    """
    def test_smiles_input_updates_model(self, main_widget, control, qtbot):
        smiles_input = main_widget.smiles_input
        for smiles in test_smiles:
            with qtbot.waitSignal(control.model.signals.updated, timeout=1000):
                control.new_molecule()  # ensure the row exists
                smiles_input.setText(smiles)
                assert control.model.__getitem__("SMILES")[control.current_index - 1] == smiles
        
        for _ in test_smiles:
            control.delete_molecule()  # clean up the model after test
        
        assert len(control.model.__getitem__("SMILES")) == 1  # only the initial empty entry remains

    def test_model_smiles_updates_smiles_input(self, main_widget, control, qtbot):
        smiles_input = main_widget.smiles_input
        for smiles in test_smiles: # note update_ui blocks signals to prevent infinite looping
            control.new_molecule()  # ensure the row exists
            control.model.__setitem__("SMILES", control.model.__getitem__("SMILES")[:-1]+ [smiles])  # update the last entry
            control.current_index = len(control.model.__getitem__("SMILES"))  # move to the last entry
            main_widget.update_ui()
            assert smiles_input.toPlainText() == smiles
        
        for _ in test_smiles:
            control.delete_molecule()  # clean up the model after test
        
        assert len(control.model.__getitem__("SMILES")) == 1  

class TestCSVTable:
    """ 
    1. Check that the CSV table displays the correct number of rows and columns based on number of entries in the model.

    2. Check that editing a code_name cell in the CSV table updates the model.

    3. Check that updating the UI with SMILES updates the CSV table cell.

    4. Check that deleting an entry in the model updates the CSV table.


    """
    def test_csv_table_row_column_count(self, control, csv_table, csv_model):
        for _ in test_smiles:
            control.new_molecule()  # add new entry to the model
            expected_rows = len(csv_model["SMILES"])
            expected_columns = len(csv_model.keys()) # should always be the same

            assert csv_table.get_row_count() == expected_rows
            assert csv_table.get_column_count() == expected_columns

        for _ in test_smiles:
            control.delete_molecule()  # remove entry from the model
            expected_rows = len(csv_model["SMILES"])
            expected_columns = len(csv_model.keys()) # should always be the same

            assert csv_table.get_row_count() == expected_rows
            assert csv_table.get_column_count() == expected_columns

    def test_csv_table_cell_edit_updates_model(self, csv_table, csv_model, qtbot):
        for code_name in test_code_names:
            row = 0 # Edit the first row
            column = list(csv_model.keys()).index("code_name")
            item = csv_table.get_item(row, column)

            with qtbot.waitSignal(csv_model.signals.updated, timeout=1000):
                item.setText(code_name)
                assert csv_model["code_name"][row] == code_name

    def test_ui_smiles_updates_csv_table(self, csv_table, main_widget, qtbot):
        # Remove all existing entries (leave one empty)
        main_widget.control.model.__setitem__("SMILES", [""])
        assert csv_table.get_row_count() == 1
        row_count = 1 # csv_table should have 1 row to begin with
        for smiles in test_smiles:
            main_widget.smiles_input.setText(smiles)
            assert csv_table.get_row_count() == row_count
            row_count += 1
            with qtbot.waitSignal(main_widget.new_molecule_button.clicked, timeout = 1000):
                qtbot.mouseClick(main_widget.new_molecule_button, Qt.MouseButton.LeftButton)

class TestPropertiesTable:
    """
    Properties table is located in the main widget. It displays details about the currently viewed molecule (or currently viewed row in the csv file). 
    """


class TestCSEARCHWorker:   
    """
    1. collect_research_params() returns dictionary with correct keys and default values when no UI inputs are set.
    2-4. collect_research_params() returns dictionary with correct keys and updated values when UI inputs are set.
    """

    def test_collect_research_params_defaults(self, csearch_widget):
        """
        For default inputs, there should only be one UI-imposed change, which is the program. Without specifying the program, aqme cannot run. 
        """
        CSEARCHworker = csearch_widget.worker
        worker = Worker(CSEARCHworker)
        params = worker.collect_csearch_params()

        assert params == {"program": "rdkit"}

    def test_collect_research_params_with_input(self, csearch_widget):
        CSEARCHworker = csearch_widget.worker
        worker = Worker(CSEARCHworker)
        CSEARCHworker.model.__setitem__("input", "test.csv")
        params = worker.collect_csearch_params()

        assert params["input"] == "test.csv"
        assert params["program"] == "rdkit"

    def test_collect_research_params_with_UI_changes(self, csearch_widget):
        """
        The model should be updated 
        """
        ...
        
        
class TestCSEARCHsetupUI:
    """
    0. Clear model and UI before running tests.
    1. Test importing the test file from test_files (test.csv) - should populate the model and CSV table correctly.
    2. Check that Run AQME button triggers user to save the file and does not proceed without saving.
    3. Check after saving that Run AQME triggers the command_model to start processing.
    4. Check that new directory is created for CSEARCH outputs.
    """

    def test_clear_model_and_ui(self, control, main_widget):
        # delete all current entries:
        control.model.__setitem__("SMILES", [""])
        assert control.model.__getitem__("SMILES") == [""]
        assert len(control.model["SMILES"]) == 1
        assert main_widget.smiles_input.toPlainText() == ""
        assert control.current_index == 1
        assert control.get_total_index() == 1


    def test_import_test_file(self, main_widget, csv_table):
        test_file_path = os.path.join(os.path.dirname(__file__), "test2.csv")
        main_widget.import_file(file_name=test_file_path)

        assert os.path.exists(test_file_path)
        assert len(main_widget.csv_model["SMILES"]) == 4  # 4 entries in test.csv
        assert csv_table.get_row_count() == 4

    def test_run_aqme_without_saving(self, csearch_widget, qtbot):
        """Test that running AQME without saving shows an error message."""
        # Capture the error signal
        error_message = None
        def capture_error(msg):
            nonlocal error_message
            error_message = msg
        
        csearch_widget.worker.error.connect(capture_error)
        
        # clear the command model input to simulate unsaved state
        csearch_widget.worker.model.__setitem__("input", "")
        csearch_widget.worker.run()
        
        qtbot.wait(100)
        assert error_message == "Please save the CSV file before running AQME."
