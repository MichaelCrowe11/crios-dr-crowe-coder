from __future__ import annotations
from rdkit import Chem
from rdkit.Chem import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from typing import Optional, Dict, Any
import logging

logger = logging.getLogger(__name__)

class MoleculeValidator:
    """Validation + standardization (salts, tautomers, canonical SMILES)."""

    def __init__(self, remove_salts: bool = True, standardize: bool = True):
        self.remove_salts = remove_salts
        self.standardize = standardize
        self._remover = SaltRemover.SaltRemover() if remove_salts else None
        self._cleaner = rdMolStandardize.CleanupParameters()

    def from_smiles(self, smiles: str) -> Chem.Mol:
        m = Chem.MolFromSmiles(smiles, sanitize=True)
        if m is None:
            raise ValueError(f"Invalid SMILES: {smiles!r}")
        Chem.SanitizeMol(m)
        if self.remove_salts and self._remover:
            m = self._remover.StripMol(m, dontRemoveEverything=True)
        if self.standardize:
            m = rdMolStandardize.Cleanup(m, self._cleaner)
            m = rdMolStandardize.Normalize(m)
            m = rdMolStandardize.Reionize(m)
            try:
                te = rdMolStandardize.TautomerEnumerator()
                m = te.Canonicalize(m)
            except Exception:
                pass
        return m

    @staticmethod
    def canonical_smiles(m: Chem.Mol) -> str:
        return Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)

    def validate_and_canonicalize(self, smiles: str) -> str:
        m = self.from_smiles(smiles)
        return self.canonical_smiles(m)