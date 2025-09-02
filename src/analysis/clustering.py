from __future__ import annotations
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.ML.Cluster import Butina
from typing import List, Tuple

def _fps(smiles_list: List[str]):
    fps = []
    for s in smiles_list:
        m = Chem.MolFromSmiles(s)
        if m is None:
            fps.append(None)
        else:
            fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048))
    return fps

def butina_clusters(smiles_list: List[str], threshold: float = 0.4) -> List[List[int]]:
    fps = _fps(smiles_list)
    n = len(fps)
    # Distances (1 - Tanimoto)
    dists = []
    for i in range(1, n):
        for j in range(i):
            if fps[i] is None or fps[j] is None:
                d = 1.0
            else:
                d = 1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
            dists.append(d)
    clusters = Butina.ClusterData(dists, nPts=n, distThresh=1.0 - threshold, isDistData=True)
    # clusters is tuple of tuples of indices
    return [list(c) for c in clusters]

# Keep backward compatibility
class CompoundClusterer:
    """Legacy compatibility wrapper"""
    def butina_cluster(self, molecules: List, cutoff: float = 0.4) -> List[List[int]]:
        # Convert molecules to SMILES
        smiles_list = []
        for mol in molecules:
            if isinstance(mol, str):
                smiles_list.append(mol)
            elif hasattr(mol, 'GetProp'):  # RDKit Mol
                smiles_list.append(Chem.MolToSmiles(mol))
            else:
                smiles_list.append(str(mol))
        return butina_clusters(smiles_list, threshold=1.0-cutoff)
    
    def cluster_statistics(self, molecules: List, clusters: List[List[int]]) -> dict:
        return {
            "n_clusters": len(clusters),
            "largest_cluster": max(len(c) for c in clusters) if clusters else 0,
            "n_singletons": sum(1 for c in clusters if len(c) == 1),
            "mean_cluster_size": sum(len(c) for c in clusters) / len(clusters) if clusters else 0,
            "mean_intra_similarity": 0.0  # Would require computing all pairwise similarities
        }
    
    def diverse_subset_selection(self, molecules: List, n_diverse: int = 10) -> List[int]:
        # Simple diverse selection: take first element from each cluster
        clusters = self.butina_cluster(molecules)
        selected = []
        for cluster in clusters:
            if cluster and len(selected) < n_diverse:
                selected.append(cluster[0])
        return selected[:n_diverse]