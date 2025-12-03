import os
import sys
import pytest
from pytestqt import qtbot
from PySide6.QtCore import Qt

from aqmeasy.models.CSEARCH_model.CSEARCH_command import command_model
from aqmeasy.ui.CSEARCH_ui.CSEARCH import CSEARCH

test_smiles = ["CCO", "C1=CC=CC=C1", "C(C(=O)O)N", "CC(=O)OC1=CC=CC=C1C(=O)O", "C1CCCCC1", "C1=CC=CN=C1", "CCN(CC)CC", "C(CCl)Br", "CC(C)O", "C1=CC=C(C=C1)O"]

test_code_names = ["ethanol", "benzene", "glycine", "aspirin", "cyclohexane", "pyridine", "triethylamine", "bromochloromethane", "isopropanol", "phenol"]

@pytest.fixture
def csearch_widget(qtbot):
    widget = CSEARCH()
    qtbot.addWidget(widget)
    return widget

# To test some components directly:
@pytest.fixture
def main_widget(csearch_widget):
    return csearch_widget.main_widget

@pytest.fixture
def csv_table(csearch_widget, qtbot, main_widget):
    """This has to be initialised before using."""
    with qtbot.waitSignal(main_widget.show_all_button.clicked, timeout=1000):
            qtbot.mouseClick(main_widget.show_all_button, Qt.MouseButton.LeftButton)
    return csearch_widget.main_widget.control.csv

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
        assert csv_table.isVisible()

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
                assert control.model["SMILES"][control.current_index - 1] == smiles
        
        for _ in test_smiles:
            control.delete_molecule()  # clean up the model after test
        
        assert len(control.model["SMILES"]) == 1  # only the initial empty entry remains

    # def test_model_updates_smiles_input(self, main_widget, control, qtbot):
    #     smiles_input = main_widget.smiles_input
    #     for idx, smiles in enumerate(test_smiles):
    #             control.new_molecule()  # ensure the row exists
    #             control.model["SMILES"][idx+1] = smiles
    #             assert smiles_input.toPlainText() == smiles

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
        row_count = 1 # csv_table should have 1 row to begin with
        for smiles in test_smiles:
            # Set the smiles
            main_widget.smiles_input.setText(smiles)
            # assert that csv_table has updated
            assert csv_table.get_row_count() == row_count
            row_count += 1
            # press Add in UI for more entries
            with qtbot.waitSignal(main_widget.new_molecule_button.clicked, timeout = 1000):
                qtbot.mouseClick(main_widget.new_molecule_button, Qt.MouseButton.LeftButton)

class TestPropertiesTable:
    """
    Properties table is located in the main widget. It displays details about the currently viewed molecule (or currently viewed row in the csv file). 
    """
# class TestKeyboardShortcuts:
#     """ 
#     1. """
#     ...