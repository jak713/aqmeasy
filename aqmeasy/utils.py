from pathlib import Path
from PySide6.QtGui import QPixmap
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

import pubchempy as pcp
import os
import sys
import tempfile

def smiles2pixmap(smiles:str) -> QPixmap:
    """Convert a SMILES string to a QPixmap image of the molecule. Only used in csv_table.
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
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp_file:
        tmp_path = tmp_file.name
    drawer.WriteDrawingText(tmp_path)
    pixmap = QPixmap(tmp_path)
    try:
        os.unlink(tmp_path)
    except:
        pass
    pixmap = pixmap.scaled(120, 120, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation)
    return pixmap

def pubchem2smiles(search_text: str) -> str | None:
    """Convert a search text (CID, or name) to a SMILES string using PubChem.
    Args:
        search_text (str): The search text (CID or name).
    Returns:
        str: The SMILES string of the compound."""
    search_text = search_text.strip()
    if not search_text:
        raise ValueError("Query cannot be empty.")
    
    try:
        if search_text.isdigit():
            compound = pcp.get_compounds(int(search_text), 'cid')[0]
        else:
            compound = pcp.get_compounds(search_text, 'name')[0]

        smiles = compound.smiles
    except:
        raise ValueError("No SMILES found for the given input.")
    return smiles

def smiles2enumerate(smiles:str) -> str:
    """Convert a SMILES string to an enumerated SMILES string.
    Args:
        smiles (str): The SMILES string of the molecule.
    Returns:
        str: The enumerated SMILES string of the molecule."""
    if not smiles or not smiles.strip():
        raise ValueError("SMILES string cannot be empty.")
    
    if smiles_enumerated(smiles):
        return smiles
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        for i, atom in enumerate(mol.GetAtoms()):
            atom.SetAtomMapNum(i + 1)
        enumerated_smiles = Chem.MolToSmiles(mol)                
        return enumerated_smiles
    except Exception:
        return "Invalid SMILES string."
    
def smiles_enumerated(smiles:str) -> bool:
    """Checks if a given smiles string is enumerated by checking for square brackets, colons and numbers.
    Args:
        smiles (str)
    Returns:
        bool
    """
    return ("]" and "[" and ":")  in smiles

def smiles2charge(smiles):
    """Convert a SMILES string to the formal charge of the molecule.
    Args:
        smiles (str): The SMILES string of the molecule.
    Returns:
        int: The formal charge of the molecule."""
    if not smiles or not smiles.strip():
        raise ValueError("SMILES string cannot be empty.")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")
    charge = Chem.GetFormalCharge(mol)
    return charge

def smiles2multiplicity(smiles: str) -> int:
    """Convert a SMILES string to the multiplicity of the molecule.
    Args:
        smiles (str): The SMILES string of the molecule.
    Returns:
        int: 2S + 1"""
    if not smiles or not smiles.strip():
        raise ValueError("SMILES string cannot be empty.")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")
    unpaired_electrons = 0
    for atom in mol.GetAtoms():
        unpaired_electrons += atom.GetNumRadicalElectrons()
    mult = unpaired_electrons + 1 
    return mult


def smiles2numatoms(smiles: str) -> int:
    """Convert a SMILES string to the number of atoms in the molecule. Adds implicit hydrogens.
    Args:
        smiles (str): The SMILES string of the molecule.
    Returns:
        int: The number of atoms in the molecule."""
    if not smiles or not smiles.strip():
        return 0
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
    except:
        raise ValueError("Invalid SMILES string.")
    
    num_atoms = mol.GetNumAtoms()
    return num_atoms

def smiles2numelectrons(smiles: str) -> int:
    """Convert a SMILES string to the number of electrons in the molecule.
    Args:
        smiles (str): The SMILES string of the molecule.
    Returns:
        int: The number of electrons in the molecule."""
    if not smiles or not smiles.strip():
        return 0
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
    except:
        raise ValueError("Invalid SMILES string.")
    
    num_electrons = sum([atom.GetAtomicNum() for atom in mol.GetAtoms()] ) - smiles2charge(smiles)
    return num_electrons

def smiles2findmetal(smiles: str) -> list:
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
    return metal_atoms

def smiles2ismetalcomplex(smiles: str) -> bool:
    """Check if a SMILES string contains transition metal atoms.
    Args:
        smiles (str): The SMILES string of the molecule.
    Returns:
        bool: True if the molecule contains transition metal atoms, False otherwise."""
    mol = Chem.MolFromSmiles(smiles)
    transition_metals = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co','Ni', 'Cu', 'Zn', 
                        'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 
                            'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 
                            'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in transition_metals:
            return True
    return False
    
def command2clipboard(command) -> bool:
    """Copy a command to the clipboard.
    Args:
        command (str): The command to copy.
    Returns:
        bool: True if the command was copied successfully, False otherwise."""
    clipboard = QApplication.clipboard()
    command = str(command)
    if command == "None":
        return False
    try:
        clipboard.setText(command)
        return True
    except Exception:
        return False


def resource_path(relative_path):
    """Get absolute path to resource, works for dev and for PyInstaller"""
    # try:
    #     # PyInstaller creates a temp folder and stores path in _MEIPASS
    #     base_path = sys._MEIPASS
    # except AttributeError:
    
    # If not running as bundled app, use the script's directory
    base_path = Path(__file__).resolve().parent
    
    return os.path.join(base_path, relative_path)