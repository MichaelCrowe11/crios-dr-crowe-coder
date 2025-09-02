from __future__ import annotations
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from typing import Iterable, Tuple, List

def fp_from_smiles(smiles: str, radius: int = 2, nbits: int = 2048):
    m = Chem.MolFromSmiles(smiles)
    if m is None: return None
    return AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=nbits)

def tanimoto_search(query_smiles: str, candidates: Iterable[Tuple[str,str]], threshold: float = 0.7) -> List[Tuple[str, float]]:
    """candidates: iterable of (mol_id, smiles)"""
    qfp = fp_from_smiles(query_smiles)
    if qfp is None: raise ValueError("Invalid query SMILES")
    hits: List[Tuple[str,float]] = []
    for mol_id, smi in candidates:
        fp = fp_from_smiles(smi)
        if fp is None: continue
        sim = DataStructs.TanimotoSimilarity(qfp, fp)
        if sim >= threshold:
            hits.append((mol_id, float(sim)))
    hits.sort(key=lambda x: x[1], reverse=True)
    return hits

# Keep backward compatibility
class CompoundDiscovery:
    """Legacy compatibility wrapper"""
    def find_similar_compounds(self, query_smiles: str, database: List[str], threshold: float = 0.7) -> List[Tuple[str, float]]:
        candidates = [(f"MOL_{i}", s) for i, s in enumerate(database)]
        hits = tanimoto_search(query_smiles, candidates, threshold)
        # Return as (smiles, score) for backward compat
        result = []
        for mol_id, score in hits:
            idx = int(mol_id.split("_")[1]) if "_" in mol_id else 0
            if idx < len(database):
                result.append((database[idx], score))
        return result