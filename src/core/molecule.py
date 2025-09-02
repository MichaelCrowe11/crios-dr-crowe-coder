from __future__ import annotations
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from typing import Optional, Dict, List
import hashlib
import logging

from .descriptors import DescriptorCalculator
from .validators import MoleculeValidator

logger = logging.getLogger(__name__)

class CriOSMolecule:
    """Lazy RDKit molecule with cached descriptors & fingerprint."""
    def __init__(self, smiles: str, mol_id: Optional[str] = None):
        self.smiles_in = smiles
        self._validator = MoleculeValidator()
        self._mol: Optional[Chem.Mol] = None
        self._canonical_smiles: Optional[str] = None
        self._fingerprint = None
        self._descriptors: Dict[str, float] = {}
        self.mol_id = mol_id or self._gen_id(smiles)

    @staticmethod
    def _gen_id(s: str) -> str:
        return f"MOL_{hashlib.sha1(s.encode('utf-8')).hexdigest()[:12]}"

    @property
    def mol(self) -> Chem.Mol:
        if self._mol is None:
            self._mol = self._validator.from_smiles(self.smiles_in)
        return self._mol

    @property
    def canonical_smiles(self) -> str:
        if self._canonical_smiles is None:
            self._canonical_smiles = self._validator.canonical_smiles(self.mol)
        return self._canonical_smiles

    def calculate_descriptors(self, names: Optional[List[str]] = None) -> Dict[str, float]:
        calc = DescriptorCalculator(names)
        self._descriptors.update(calc.calculate(self.mol))
        return dict(self._descriptors)

    def fingerprint(self, radius: int = 2, nbits: int = 2048):
        if self._fingerprint is None:
            self._fingerprint = AllChem.GetMorganFingerprintAsBitVect(self.mol, radius, nBits=nbits)
        return self._fingerprint

    def tanimoto(self, other: "CriOSMolecule") -> float:
        return DataStructs.TanimotoSimilarity(self.fingerprint(), other.fingerprint())

    def passes_filter(self, filter_name: str = "Lipinski") -> bool:
        from .filters import MolecularFilter  # lazy import to avoid cycles
        return MolecularFilter(filter_name).evaluate(self)

    def to_row(self) -> Dict[str, object]:
        row = {"mol_id": self.mol_id, "smiles": self.canonical_smiles}
        row.update(self._descriptors)
        return row