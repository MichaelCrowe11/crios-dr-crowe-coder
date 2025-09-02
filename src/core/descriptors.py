from __future__ import annotations
from typing import List, Dict, Any, Optional, Iterable
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
import functools
import logging

logger = logging.getLogger(__name__)

# Keep names user-facing & stable
_DESCRIPTOR_FUNCS = {
    "MW": lambda m: Descriptors.MolWt(m),
    "LogP": lambda m: Crippen.MolLogP(m),
    "TPSA": lambda m: rdMolDescriptors.CalcTPSA(m),
    "HBA": lambda m: rdMolDescriptors.CalcNumHBA(m),
    "HBD": lambda m: rdMolDescriptors.CalcNumHBD(m),
    "RotatableBonds": lambda m: rdMolDescriptors.CalcNumRotatableBonds(m),
    "FractionCSP3": lambda m: rdMolDescriptors.CalcFractionCSP3(m),
    "HeavyAtomCount": lambda m: rdMolDescriptors.CalcNumHeavyAtoms(m),
}

class DescriptorCalculator:
    """Cached descriptor computation with batch support."""

    def __init__(self, names: Optional[List[str]] = None):
        self.names = names or list(_DESCRIPTOR_FUNCS.keys())

    def calculate(self, mol: Chem.Mol) -> Dict[str, float]:
        out: Dict[str, float] = {}
        for name in self.names:
            fn = _DESCRIPTOR_FUNCS.get(name)
            if fn is None:
                logger.warning("Unknown descriptor: %s", name); continue
            try:
                out[name] = float(fn(mol))
            except Exception as e:
                logger.exception("Descriptor %s failed: %s", name, e)
                out[name] = float("nan")
        return out

    def batch(self, mols: Iterable[Chem.Mol]) -> List[Dict[str, float]]:
        return [self.calculate(m) for m in mols]