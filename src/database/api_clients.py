from __future__ import annotations
import requests
from typing import List, Tuple

# PubChem: simple text query -> SMILES
def pubchem_smiles(query: str, limit: int = 100) -> List[Tuple[str,str]]:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(query)}/cids/TXT"
    r = requests.get(url, timeout=30); r.raise_for_status()
    cids = [line.strip() for line in r.text.splitlines() if line.strip()]
    out: List[Tuple[str,str]] = []
    for cid in cids[:limit]:
        s_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES/TXT"
        sr = requests.get(s_url, timeout=30)
        if sr.status_code == 200:
            smi = sr.text.strip()
            out.append((f"PUBCHEM_{cid}", smi))
    return out

# ChEMBL: simple target-based fetch (bioactive molecules -> canonical SMILES)
def chembl_by_target(chembl_target_id: str, limit: int = 100) -> List[Tuple[str,str]]:
    base = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
    params = {"target_chembl_id": chembl_target_id, "limit": str(limit)}
    r = requests.get(base, params=params, timeout=30); r.raise_for_status()
    items = r.json().get("activities", [])
    out: List[Tuple[str,str]] = []
    for a in items:
        smi = a.get("canonical_smiles")
        cid = a.get("molecule_chembl_id") or a.get("record_id") or "CHEMBL_UNK"
        if smi:
            out.append((f"CHEMBL_{cid}", smi))
    return out

# Keep backward compatibility
class CompoundDatabase:
    """Legacy compatibility wrapper"""
    def __init__(self):
        pass
    
    class PubChemClient:
        def search_by_name(self, query: str, limit: int = 100):
            results = pubchem_smiles(query, limit)
            return [{"cid": mol_id.replace("PUBCHEM_", ""), "smiles": smi, "source": "pubchem"} 
                   for mol_id, smi in results]
    
    class ChEMBLClient:
        def search_by_target(self, target: str, limit: int = 100):
            results = chembl_by_target(target, limit)
            return [{"chembl_id": mol_id.replace("CHEMBL_", ""), "smiles": smi, "source": "chembl"}
                   for mol_id, smi in results]
    
    def __init__(self):
        self.pubchem = self.PubChemClient()
        self.chembl = self.ChEMBLClient()
    
    def search(self, query: str, source: str = "pubchem", limit: int = 100):
        if source == "pubchem":
            return self.pubchem.search_by_name(query, limit)
        elif source == "chembl":
            # For general search, this is limited
            return []
        return []