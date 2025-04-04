from PySide6.QtGui import QPixmap
from PySide6.QtCore import Qt

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

import pubchempy as pcp


def smiles2pixmap(smiles):
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