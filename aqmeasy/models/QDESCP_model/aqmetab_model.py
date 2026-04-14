from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.rdmolfiles import MolsFromCDXMLFile
from rdkit.Chem.rdmolops import GetMolFrags


METAL_ATOMIC_NUMBERS = {
    3, 11, 19, 37, 55, 87,
    4, 12, 20, 38, 56, 88,
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    72, 73, 74, 75, 76, 77, 78, 79, 80,
    13, 49, 50, 81, 82, 83,
}


@dataclass(slots=True)
class MappedSmilesResult:
    mapped_smiles: list[Optional[str]]
    ambiguous_smiles: list[str]


@dataclass(slots=True)
class CommonSmartsResult:
    smarts: str
    atom_count: int
    bond_count: int


def smart_read_csv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path)


def detect_smiles_column(df: pd.DataFrame) -> Optional[str]:
    return next((column for column in df.columns if str(column).strip().lower() == "smiles"), None)


def load_molecules_from_path(path: str | Path) -> list[Chem.Mol]:
    file_path = Path(path)
    suffix = file_path.suffix.lower()

    if suffix == ".cdxml":
        mols = MolsFromCDXMLFile(str(file_path), sanitize=False, removeHs=False)
        valid_mols: list[Chem.Mol] = []
        for mol in mols:
            if mol is None:
                continue
            fragments = GetMolFrags(mol, asMols=True, sanitizeFrags=False)
            valid_mols.extend(fragment for fragment in fragments if fragment is not None)
        return valid_mols

    if suffix == ".sdf":
        return [mol for mol in Chem.SDMolSupplier(str(file_path)) if mol is not None]

    if suffix == ".mol":
        mol = Chem.MolFromMolFile(str(file_path))
        return [mol] if mol is not None else []

    if suffix == ".cdx":
        raise ValueError(
            "CDX import is not supported for automatic parsing. Export the file as CDXML first."
        )

    mol = Chem.MolFromMolFile(str(file_path))
    return [mol] if mol is not None else []


def find_common_smarts_details(smiles_list: Iterable[str], timeout_seconds: int = 30) -> CommonSmartsResult:
    mols = []
    for smiles in smiles_list:
        if not smiles:
            continue
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is not None:
            mols.append(mol)

    if not mols:
        raise ValueError("No valid SMILES strings were provided.")

    timeout_seconds = max(1, int(timeout_seconds))
    result = rdFMCS.FindMCS(mols, timeout=timeout_seconds)

    if getattr(result, "canceled", False):
        raise TimeoutError("MCS search timed out before a common SMARTS pattern was found.")

    if not result.smartsString:
        raise ValueError("No common SMARTS pattern could be found.")

    return CommonSmartsResult(
        smarts=result.smartsString,
        atom_count=int(getattr(result, "numAtoms", 0)),
        bond_count=int(getattr(result, "numBonds", 0)),
    )


def find_common_smarts(smiles_list: Iterable[str], timeout_seconds: int = 30) -> str:
    return find_common_smarts_details(smiles_list, timeout_seconds=timeout_seconds).smarts


def generate_mapped_smiles(
    smarts_pattern: str,
    selected_pattern_indices: Iterable[int],
    smiles_list: Iterable[str],
) -> MappedSmilesResult:
    pattern_mol = Chem.MolFromSmarts(smarts_pattern)
    if pattern_mol is None:
        raise ValueError("Invalid SMARTS pattern")

    mapped_smiles: list[Optional[str]] = []
    ambiguous_smiles: list[str] = []

    for smiles in smiles_list:
        mol = Chem.AddHs(Chem.MolFromSmiles(str(smiles)))
        if mol is None:
            mapped_smiles.append(None)
            continue

        matches = mol.GetSubstructMatches(pattern_mol)
        if len(matches) > 1:
            mapped_smiles.append(None)
            ambiguous_smiles.append(str(smiles))
            continue

        if not matches:
            mapped_smiles.append(None)
            continue

        match = matches[0]
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)

        for map_index, pattern_index in enumerate(selected_pattern_indices):
            if 0 <= pattern_index < len(match):
                mol.GetAtomWithIdx(match[pattern_index]).SetAtomMapNum(map_index + 1)

        mapped_smiles.append(Chem.MolToSmiles(mol))

    return MappedSmilesResult(mapped_smiles=mapped_smiles, ambiguous_smiles=ambiguous_smiles)


def extract_qdescp_prefill_from_sdf(file_paths: Iterable[str]) -> dict[str, int]:
    for file_path in file_paths:
        supplier = Chem.SDMolSupplier(str(file_path))
        for mol in supplier:
            if mol is None:
                continue
            charge = int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))
            radical_electrons = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
            multiplicity = max(1, radical_electrons + 1)
            return {"charge": charge, "mult": multiplicity}

    return {"charge": 0, "mult": 1}
