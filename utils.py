from PySide6.QtGui import QPixmap
from PySide6.QtCore import Qt

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D

import pubchempy as pcp


def smiles2pixmap(smiles):
    """Convert a SMILES string to a QPixmap image of the molecule.
    Args:
        smiles (str): The SMILES string of the molecule.
    Returns:
        QPixmap: A QPixmap image of the molecule."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")
    drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    drawer.WriteDrawingText("/tmp/molecule.png")
    pixmap = QPixmap("/tmp/molecule.png")
    pixmap = pixmap.scaled(120, 120, Qt.KeepAspectRatio, Qt.SmoothTransformation)
    return pixmap

def pubchem2smiles(search_text):
    """Convert a search text (CAS, CID, or name) to a SMILES string using PubChem.
    Args:
        search_text (str): The search text (CAS, CID, or name).
    Returns:
        str: The SMILES string of the compound."""
    search_text = search_text.strip()
    if not search_text:
        
        return
    if "-" in search_text:
        compounds = pcp.get_compounds(search_text, 'name')
        if compounds:
            cid = compounds[0].cid
            compound = pcp.get_compounds(cid, 'cid')[0]
        else:
            raise ValueError("No CID match for the given CAS.")
    elif search_text.isdigit():
        compound = pcp.get_compounds(int(search_text), 'cid')[0]
    else:
        compound = pcp.get_compounds(search_text, 'name')[0]

    smiles = compound.isomeric_smiles
    if not smiles:
        raise ValueError("No SMILES found for the given input.")

    return smiles

def smiles2enumerate(smiles):
    """Convert a SMILES string to an enumerated SMILES string.
    Args:
        smiles (str): The SMILES string of the molecule.
    Returns:
        str: The enumerated SMILES string of the molecule."""
    if not smiles:
        return
    elif smiles == "":
        return
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        for i, atom in enumerate(mol.GetAtoms()):
            atom.SetAtomMapNum(i + 1)
        enumerated_smiles = Chem.MolToSmiles(mol)                
        return enumerated_smiles
    except Exception:
        return

def smiles2charge(smiles):
    """Convert a SMILES string to the formal charge of the molecule.
    Args:
        smiles (str): The SMILES string of the molecule.
    Returns:
        int: The formal charge of the molecule."""
    if not smiles:
        return
    elif smiles == "":
        return
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")
    charge = Chem.GetFormalCharge(mol)
    return charge

# if not hasattr(self, 'user_defined_multiplicity'):
#             self.user_defined_multiplicity = {}
        
#         if self.current_index in self.user_defined_multiplicity:
#             return

# if not hasattr(self, 'user_defined_charge'):
#             self.user_defined_charge = {}
        
#         if self.current_index in self.user_defined_charge:
#             return

def smiles2multiplicity(smiles):
    """Convert a SMILES string to the multiplicity of the molecule.
    Args:
        smiles (str): The SMILES string of the molecule.
    Returns:
        int: The multiplicity of the molecule."""
    if not smiles:
        return
    elif smiles == "":
        return
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return
        # raise ValueError("Invalid SMILES string.")
    multiplicity = Descriptors.NumRadicalElectrons(mol) + 1
    return multiplicity

def smiles2numatoms(smiles):
    """Convert a SMILES string to the number of atoms in the molecule.
    Args:
        smiles (str): The SMILES string of the molecule.
    Returns:
        int: The number of atoms in the molecule."""
    if not smiles:
        return 0
    elif smiles == "":
        return 0
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return
        # raise ValueError("Invalid SMILES string.")
    num_atoms = mol.GetNumAtoms()
    return num_atoms

def smiles2numelectrons(smiles):
    """Convert a SMILES string to the number of electrons in the molecule.
    Args:
        smiles (str): The SMILES string of the molecule.
    Returns:
        int: The number of electrons in the molecule."""
    if not smiles: 
        return 0
    elif smiles == "":
        return 0 
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")
    num_electrons = sum([atom.GetAtomicNum() for atom in mol.GetAtoms()])
    return num_electrons

def smiles2findmetal(smiles):
    """Find the transition metal atoms in a SMILES string.
    Args:
        smiles (str): The SMILES string of the molecule.
    Returns:
        list: A list of transition metal symbols found in the molecule."""
    mol = Chem.MolFromSmiles(smiles)
    metal_atoms = [] 
    transition_metals = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co','Ni', 'Cu', 'Zn', 
                        'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 
                            'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 
                            'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in transition_metals:
            metal_atoms.append(atom.GetSymbol())
        if len(metal_atoms) > 0:
            return metal_atoms
    else:
        return 
    
def import_file():
    pass