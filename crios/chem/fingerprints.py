"""
CriOS Fingerprints and Similarity
High-performance molecular fingerprint generation and similarity calculations
"""

import logging
from typing import List, Optional, Union, Dict, Any, Tuple
import warnings

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem
from rdkit import DataStructs
from scipy.spatial.distance import pdist, squareform

from ..core.exceptions import CriosError

logger = logging.getLogger(__name__)


class FingerprintError(CriosError):
    """Fingerprint-specific errors"""
    pass


class SimilarityError(CriosError):
    """Similarity calculation errors"""
    pass


def generate_morgan_fingerprint(
    mol: Chem.Mol,
    radius: int = 2,
    n_bits: int = 2048,
    use_features: bool = False,
    use_chirality: bool = False
) -> Optional[np.ndarray]:
    """
    Generate Morgan (ECFP-like) fingerprint
    
    Args:
        mol: RDKit molecule
        radius: Fingerprint radius
        n_bits: Number of bits
        use_features: Use feature-based fingerprint
        use_chirality: Include chirality information
    
    Returns:
        Binary fingerprint as numpy array
    """
    if mol is None:
        return None
    
    try:
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol, 
            radius=radius, 
            nBits=n_bits,
            useFeatures=use_features,
            useChirality=use_chirality
        )
        
        # Convert to numpy array
        arr = np.zeros(n_bits, dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        
        return arr
        
    except Exception as e:
        logger.error(f"Morgan fingerprint generation failed: {e}")
        return None


def generate_maccs_fingerprint(mol: Chem.Mol) -> Optional[np.ndarray]:
    """
    Generate MACCS fingerprint
    
    Args:
        mol: RDKit molecule
    
    Returns:
        Binary fingerprint as numpy array
    """
    if mol is None:
        return None
    
    try:
        fp = rdMolDescriptors.GetMACCSKeysFingerprint(mol)
        
        # Convert to numpy array
        arr = np.zeros(167, dtype=np.uint8)  # MACCS has 166 bits + 1 for indexing
        DataStructs.ConvertToNumpyArray(fp, arr)
        
        return arr[1:]  # Remove the first bit (not used)
        
    except Exception as e:
        logger.error(f"MACCS fingerprint generation failed: {e}")
        return None


def generate_rdkit_fingerprint(
    mol: Chem.Mol,
    fp_size: int = 2048,
    min_path: int = 1,
    max_path: int = 7
) -> Optional[np.ndarray]:
    """
    Generate RDKit fingerprint
    
    Args:
        mol: RDKit molecule
        fp_size: Fingerprint size
        min_path: Minimum path length
        max_path: Maximum path length
    
    Returns:
        Binary fingerprint as numpy array
    """
    if mol is None:
        return None
    
    try:
        fp = Chem.RDKFingerprint(
            mol,
            fpSize=fp_size,
            minPath=min_path,
            maxPath=max_path
        )
        
        # Convert to numpy array
        arr = np.zeros(fp_size, dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        
        return arr
        
    except Exception as e:
        logger.error(f"RDKit fingerprint generation failed: {e}")
        return None


def generate_fingerprint(
    mol: Chem.Mol,
    fp_type: str = "morgan",
    **kwargs
) -> Optional[np.ndarray]:
    """
    Generate fingerprint of specified type
    
    Args:
        mol: RDKit molecule
        fp_type: Fingerprint type ("morgan", "maccs", "rdkit")
        **kwargs: Type-specific parameters
    
    Returns:
        Binary fingerprint as numpy array
    """
    fp_type = fp_type.lower()
    
    if fp_type == "morgan":
        return generate_morgan_fingerprint(mol, **kwargs)
    elif fp_type == "maccs":
        return generate_maccs_fingerprint(mol)
    elif fp_type == "rdkit":
        return generate_rdkit_fingerprint(mol, **kwargs)
    else:
        raise FingerprintError(f"Unknown fingerprint type: {fp_type}")


def generate_fingerprints_batch(
    mols: List[Chem.Mol],
    fp_type: str = "morgan",
    **kwargs
) -> List[Optional[np.ndarray]]:
    """
    Generate fingerprints for multiple molecules
    
    Args:
        mols: List of RDKit molecules
        fp_type: Fingerprint type
        **kwargs: Type-specific parameters
    
    Returns:
        List of fingerprints (None for failed)
    """
    fingerprints = []
    
    for mol in mols:
        fp = generate_fingerprint(mol, fp_type, **kwargs)
        fingerprints.append(fp)
    
    return fingerprints


def tanimoto_similarity(
    fp1: np.ndarray,
    fp2: np.ndarray
) -> float:
    """
    Calculate Tanimoto similarity between two fingerprints
    
    Args:
        fp1: First fingerprint
        fp2: Second fingerprint
    
    Returns:
        Tanimoto similarity coefficient (0-1)
    """
    if fp1 is None or fp2 is None:
        return 0.0
    
    if len(fp1) != len(fp2):
        raise SimilarityError(f"Fingerprint size mismatch: {len(fp1)} vs {len(fp2)}")
    
    intersection = np.sum(fp1 & fp2)
    union = np.sum(fp1 | fp2)
    
    if union == 0:
        return 1.0 if intersection == 0 else 0.0
    
    return float(intersection) / float(union)


def dice_similarity(
    fp1: np.ndarray,
    fp2: np.ndarray
) -> float:
    """
    Calculate Dice similarity between two fingerprints
    
    Args:
        fp1: First fingerprint
        fp2: Second fingerprint
    
    Returns:
        Dice similarity coefficient (0-1)
    """
    if fp1 is None or fp2 is None:
        return 0.0
    
    if len(fp1) != len(fp2):
        raise SimilarityError(f"Fingerprint size mismatch: {len(fp1)} vs {len(fp2)}")
    
    intersection = np.sum(fp1 & fp2)
    total = np.sum(fp1) + np.sum(fp2)
    
    if total == 0:
        return 1.0 if intersection == 0 else 0.0
    
    return 2.0 * float(intersection) / float(total)


def cosine_similarity(
    fp1: np.ndarray,
    fp2: np.ndarray
) -> float:
    """
    Calculate cosine similarity between two fingerprints
    
    Args:
        fp1: First fingerprint
        fp2: Second fingerprint
    
    Returns:
        Cosine similarity coefficient (0-1)
    """
    if fp1 is None or fp2 is None:
        return 0.0
    
    if len(fp1) != len(fp2):
        raise SimilarityError(f"Fingerprint size mismatch: {len(fp1)} vs {len(fp2)}")
    
    # Convert to float for calculation
    fp1_f = fp1.astype(np.float64)
    fp2_f = fp2.astype(np.float64)
    
    dot_product = np.dot(fp1_f, fp2_f)
    norm1 = np.linalg.norm(fp1_f)
    norm2 = np.linalg.norm(fp2_f)
    
    if norm1 == 0 or norm2 == 0:
        return 1.0 if np.array_equal(fp1, fp2) else 0.0
    
    return float(dot_product) / (float(norm1) * float(norm2))


def calculate_similarity(
    fp1: np.ndarray,
    fp2: np.ndarray,
    metric: str = "tanimoto"
) -> float:
    """
    Calculate similarity using specified metric
    
    Args:
        fp1: First fingerprint
        fp2: Second fingerprint
        metric: Similarity metric ("tanimoto", "dice", "cosine")
    
    Returns:
        Similarity coefficient (0-1)
    """
    metric = metric.lower()
    
    if metric == "tanimoto":
        return tanimoto_similarity(fp1, fp2)
    elif metric == "dice":
        return dice_similarity(fp1, fp2)
    elif metric == "cosine":
        return cosine_similarity(fp1, fp2)
    else:
        raise SimilarityError(f"Unknown similarity metric: {metric}")


def tanimoto_similarity_matrix(
    fingerprints: List[np.ndarray],
    symmetric: bool = True
) -> np.ndarray:
    """
    Calculate pairwise Tanimoto similarity matrix
    
    Args:
        fingerprints: List of fingerprints
        symmetric: Assume matrix is symmetric (optimization)
    
    Returns:
        Similarity matrix
    """
    n = len(fingerprints)
    if n == 0:
        return np.array([])
    
    # Filter out None fingerprints
    valid_fps = []
    valid_indices = []
    for i, fp in enumerate(fingerprints):
        if fp is not None:
            valid_fps.append(fp)
            valid_indices.append(i)
    
    if not valid_fps:
        return np.zeros((n, n))
    
    # Calculate similarity matrix for valid fingerprints
    n_valid = len(valid_fps)
    sim_matrix = np.zeros((n_valid, n_valid))
    
    for i in range(n_valid):
        for j in range(i if symmetric else 0, n_valid):
            if i == j:
                sim_matrix[i, j] = 1.0
            else:
                sim = tanimoto_similarity(valid_fps[i], valid_fps[j])
                sim_matrix[i, j] = sim
                if symmetric:
                    sim_matrix[j, i] = sim
    
    # Expand to full matrix (including None fingerprints)
    full_matrix = np.zeros((n, n))
    for i, idx_i in enumerate(valid_indices):
        for j, idx_j in enumerate(valid_indices):
            full_matrix[idx_i, idx_j] = sim_matrix[i, j]
    
    return full_matrix


def bulk_tanimoto_similarity(
    query_fp: np.ndarray,
    database_fps: List[np.ndarray],
    threshold: float = 0.0
) -> List[Tuple[int, float]]:
    """
    Fast bulk similarity calculation against database
    
    Args:
        query_fp: Query fingerprint
        database_fps: List of database fingerprints
        threshold: Minimum similarity threshold
    
    Returns:
        List of (index, similarity) tuples above threshold
    """
    if query_fp is None:
        return []
    
    results = []
    
    for i, db_fp in enumerate(database_fps):
        if db_fp is not None:
            sim = tanimoto_similarity(query_fp, db_fp)
            if sim >= threshold:
                results.append((i, sim))
    
    # Sort by similarity (descending)
    results.sort(key=lambda x: x[1], reverse=True)
    
    return results


def nearest_neighbors(
    query_fp: np.ndarray,
    database_fps: List[np.ndarray],
    n_neighbors: int = 10,
    metric: str = "tanimoto"
) -> List[Tuple[int, float]]:
    """
    Find nearest neighbors by similarity
    
    Args:
        query_fp: Query fingerprint
        database_fps: Database fingerprints
        n_neighbors: Number of neighbors to return
        metric: Similarity metric
    
    Returns:
        List of (index, similarity) tuples for top neighbors
    """
    if query_fp is None or not database_fps:
        return []
    
    similarities = []
    
    for i, db_fp in enumerate(database_fps):
        if db_fp is not None:
            sim = calculate_similarity(query_fp, db_fp, metric)
            similarities.append((i, sim))
    
    # Sort by similarity (descending) and take top n
    similarities.sort(key=lambda x: x[1], reverse=True)
    
    return similarities[:n_neighbors]


def diversity_selection(
    fingerprints: List[np.ndarray],
    n_select: int,
    metric: str = "tanimoto",
    seed: Optional[int] = None
) -> List[int]:
    """
    Select diverse subset using MaxMin algorithm
    
    Args:
        fingerprints: List of fingerprints
        n_select: Number to select
        metric: Distance metric
        seed: Random seed for reproducibility
    
    Returns:
        Indices of selected fingerprints
    """
    if seed is not None:
        np.random.seed(seed)
    
    valid_indices = [i for i, fp in enumerate(fingerprints) if fp is not None]
    
    if len(valid_indices) <= n_select:
        return valid_indices
    
    selected = []
    
    # Start with random fingerprint
    first_idx = np.random.choice(valid_indices)
    selected.append(first_idx)
    
    while len(selected) < n_select:
        max_min_dist = -1
        best_candidate = None
        
        for candidate_idx in valid_indices:
            if candidate_idx in selected:
                continue
            
            # Calculate minimum distance to selected fingerprints
            min_dist = float('inf')
            for selected_idx in selected:
                sim = calculate_similarity(
                    fingerprints[candidate_idx],
                    fingerprints[selected_idx],
                    metric
                )
                dist = 1.0 - sim  # Convert similarity to distance
                min_dist = min(min_dist, dist)
            
            # Select candidate with maximum minimum distance
            if min_dist > max_min_dist:
                max_min_dist = min_dist
                best_candidate = candidate_idx
        
        if best_candidate is not None:
            selected.append(best_candidate)
        else:
            break
    
    return selected


# Convenience functions for SMILES
def smiles_to_fingerprint(
    smiles: str,
    fp_type: str = "morgan",
    **kwargs
) -> Optional[np.ndarray]:
    """Convert SMILES to fingerprint"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return generate_fingerprint(mol, fp_type, **kwargs)
    except Exception:
        return None


def smiles_similarity(
    smiles1: str,
    smiles2: str,
    fp_type: str = "morgan",
    metric: str = "tanimoto",
    **kwargs
) -> float:
    """Calculate similarity between two SMILES strings"""
    fp1 = smiles_to_fingerprint(smiles1, fp_type, **kwargs)
    fp2 = smiles_to_fingerprint(smiles2, fp_type, **kwargs)
    
    if fp1 is None or fp2 is None:
        return 0.0
    
    return calculate_similarity(fp1, fp2, metric)