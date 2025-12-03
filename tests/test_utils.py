import pytest
from aqmeasy import utils

valid_smiles = ["CCC", "CCO", "c1ccccc1", "C1CCCCC1", "CC(=O)O", "CNC", "CCN(CC)CC", "C1=CC=CN=C1", "C[H+]C", "***********(*****)******","C.C.C", "CC[13C]CC" ""]

invalid_smiles = ["ccc", "CcS", "c1cccc1", "1111", "C+"]

valid_query_smiles = {
    "2244":"CC(=O)OC1=CC=CC=C1C(=O)O",
    "Aspirin":"CC(=O)OC1=CC=CC=C1C(=O)O",
    }

invalid_queries = ["", "Bad Query"]

stereoisomeric_smiles = {
            "C/C=C/C":"[C:1](/[C:2](=[C:3](/[C:4]([H:10])([H:11])[H:12])[H:9])[H:8])([H:5])([H:6])[H:7]",

            "N[C@@H](C)C(=O)O": "[N:1]([C@@:2]([C:3]([H:10])([H:11])[H:12])([C:4](=[O:5])[O:6][H:13])[H:9])([H:7])[H:8]", 

            "N[CH](C)C(=O)O":"[N:1]([C:2]([C:3]([H:10])([H:11])[H:12])([C:4](=[O:5])[O:6][H:13])[H:9])([H:7])[H:8]"
            }

bonds_smiles = {
                "C-C-O":"[C:1]([C:2]([O:3][H:9])([H:7])[H:8])([H:4])([H:5])[H:6]", 

                "O=C=O": "[O:1]=[C:2]=[O:3]", 

                "C#N": "[C:1](#[N:2])[H:3]", 

                "[Na+].[Cl-]": "[Cl-:2].[Na+:1]", 

                "CCO.O": "[C:1]([C:2]([O:3][H:10])([H:8])[H:9])([H:5])([H:6])[H:7].[O:4]([H:11])[H:12]"
                }

enumerated_smiles = ["[C:1]([C:2]([C:3]([C:4]([H:12])([H:13])[H:14])([H:10])[H:11])([H:8])[H:9])([H:5])([H:6])[H:7]"]

explicit_smiles = ["Li", "[2H]O[2H]", "N[C@@H](C)C(=O)O", "N[CH](C)C(=O)O"]

neutral_smiles = {
        "CCC" : 0,
        "CCO" : 0,
        "c1ccccc1" : 0,

}
positive_smiles = {
        "[C+]CS" : 1,
        "[Co++]" : 2,
        "[Co+2]" : 2,
        "[NH4+]" : 1,
        "[Co++]CC" : 2,
}
negative_smiles = {
        "CC(=O)[O-]" : -1,
        "[Cl-]" : -1,
        "[O-]C=O" : -1,
}

non_bonded_smiles = {
        "CCC[C+].CCC[O-]": 0,
        "[K+]c1ccccc1":1,
}

singlet_smiles = {
        "CC" : 1,
        "N1C=C[NH+]=C1" : 1,
}

doublet_smiles = {
        "[CH3]" : 2,
}

triplet_smiles = {
        "[O][O]" : 3,
        "[CH2]" : 3
}

atoms_smiles = {
    "CCO": 1+3+1+2+1+1,
    "c1ccccc1": 6+6*1,
    "C1CCCCC1": 6+6*2,
    "CC(=O)O": 1+3+1+1+1+1,
    "CNC": 1+3+1+1+1+3,
}

electron_smiles = {
    "[N][O]": 7+8,
    "[O][O]": 8+8,
    "CCCCCCC": 6*7+3*2*1+2*5*1,
    "[Ga+].[As-]": 31-1 + 33+1,
    "[Ti+4].[O-]": 22-4 + 8+1,
    "[O-1][C][N+2]": 8+1 + 6 + 7-2,
    "": 0,
}

class TestUtils:
    """
    smiles2pixmap:
    1. Valid smiles string should return a pixmap
    2. Empty smiles string should return empty pixmap
    3. Wrong smiles string should raise ValueError

    pubchem2smiles:
    4. Valid CID/name -> returns smiles (this is only really testing the API which should be thoroughly tested in itself)
    5. Invalid CID/name -> raises ValueError
    6. Empty query -> raises ValueError

    smiles2enumerate:
    This poses issues for some valid smiles which are unable to be enumerated, especially for old ChemDraw file formats.
    7. Stereoisomeric smiles
    8. Bond-specific smiles

    smiles_enumerated:
    9. Valid enumerated smiles
    10. smiles with [] but not :

    smiles2charge:
    11. Neutral smiles
    12. Positive charge smiles
    13. Negative charge smiles
    14. Non-bonded smiles
    15. Empty smiles -> raises ValueError
    16. Invalid smiles -> raises ValueError 

    smiles2multiplicity:
    17. Singlet 
    18. Doublet
    19. Triplet
    20. Invalid smiles -> raises ValueError

    smiles2numatoms:
    21. Valid smiles -> returns correct number of atoms
    22. Empty smiles -> returns 0
    23. Invalid smiles -> raises ValueError

    smiles2numelectrons:
    24. Valid smiles -> returns correct number of electrons
    25. Empty smiles -> returns 0
    26. Invalid smiles -> raises ValueError

    smiles2findmetal ?
    
    """

    def test_valid_smiles2pixmap(self):
        for smiles in valid_smiles:
            print(f"Input SMILES: {smiles}")
            assert utils.smiles2pixmap(smiles)is not None

    def test_empty_smiles2pixmap(self):
            assert utils.smiles2pixmap("")is not None

    def test_invalid_smiles2pixmap(self):
            for smiles in invalid_smiles:
                with pytest.raises(ValueError): 
                    utils.smiles2pixmap(smiles)

    def test_valid_query_pubchem2smiles(self):
        for query, smiles in valid_query_smiles.items():
            assert utils.pubchem2smiles(query) == smiles

    def test_invalid_query_pubchem2smiles(self):
        for query in invalid_queries:
            with pytest.raises(ValueError):
                utils.pubchem2smiles(query)

    def test_empty_query_pubchem2smiles(self):
        with pytest.raises(ValueError):
            utils.pubchem2smiles("")

    def test_stereoisomeric_smiles2enumerate(self):
        for smiles, enumerated in stereoisomeric_smiles.items():
            assert utils.smiles_enumerated(smiles) is False
            result = utils.smiles2enumerate(smiles)
            assert utils.smiles_enumerated(result) is True
            assert result == enumerated

    def test_bonds_smiles2enumerate(self):
        for smiles, enumerated in bonds_smiles.items():
            assert utils.smiles_enumerated(smiles) is False
            result = utils.smiles2enumerate(smiles)
            assert utils.smiles_enumerated(result) is True
            assert result == enumerated

    def test_valid_smiles_enumerated(self):
        for smiles in enumerated_smiles:
            assert utils.smiles_enumerated(smiles) is True

    def test_invalid_smiles_enumerated(self):
        for smiles in explicit_smiles:
            assert utils.smiles_enumerated(smiles) is False

    def test_smiles2charge_neutral(self):
        for smiles, charge in neutral_smiles.items():
            assert utils.smiles2charge(smiles) == charge

    def test_smiles2charge_positive(self):
        for smiles, charge in positive_smiles.items():
            print(f"Input SMILES: {smiles}")
            assert utils.smiles2charge(smiles) == charge

    def test_smiles2charge_negative(self):
        for smiles, charge in negative_smiles.items():
            assert utils.smiles2charge(smiles) == charge

    def test_smiles2charge_non_bonded(self):
        for smiles, charge in non_bonded_smiles.items():
            assert utils.smiles2charge(smiles) == charge

    def test_empty_smiles2charge(self):
        with pytest.raises(ValueError):
            utils.smiles2charge("")
    
    def test_invalid_smiles2charge(self):
        for smiles in invalid_smiles:
            with pytest.raises(ValueError):
                utils.smiles2charge(smiles)

    def test_smiles2multiplicity_singlet(self):
        for smiles, multiplicity in singlet_smiles.items():
            assert utils.smiles2multiplicity(smiles) == multiplicity

    def test_smiles2multiplicity_doublet(self):
        for smiles, multiplicity in doublet_smiles.items():
            assert utils.smiles2multiplicity(smiles) == multiplicity
    
    def test_smiles2multiplicity_triplet(self):
        for smiles, multiplicity in triplet_smiles.items():
            assert utils.smiles2multiplicity(smiles) == multiplicity

    def test_invalid_smiles2multiplicity(self):
        for smiles in invalid_smiles:
            with pytest.raises(ValueError):
                utils.smiles2multiplicity(smiles)

    def test_smiles2numatoms_valid(self):
        for smiles, num_atoms in atoms_smiles.items():
            assert utils.smiles2numatoms(smiles) == num_atoms

    def test_empty_smiles2numatoms(self):
        assert utils.smiles2numatoms("") == 0

    def test_invalid_smiles2numatoms(self):
        for smiles in invalid_smiles:
            with pytest.raises(ValueError):
                utils.smiles2numatoms(smiles)

    def test_smiles2numelectrons_valid(self):
        for smiles, num_electrons in electron_smiles.items():
            print(f"Input SMILES: {smiles}")
            assert utils.smiles2numelectrons(smiles) == num_electrons

    def test_empty_smiles2numelectrons(self):
        assert utils.smiles2numelectrons("") == 0

    def test_invalid_smiles2numelectrons(self):
        for smiles in invalid_smiles:
            with pytest.raises(ValueError):
                utils.smiles2numelectrons(smiles)