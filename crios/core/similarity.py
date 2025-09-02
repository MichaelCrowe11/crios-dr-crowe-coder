"""
CriOS Similarity Search Engine
High-performance molecular similarity calculations and clustering
"""

import logging
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass
import time

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.metrics.pairwise import cosine_similarity
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina

from .molecule import Molecule, MoleculeCollection
from .compound import Compound, CompoundLibrary
from ..config.schemas import FingerprintType, ClusteringAlgorithm
from ..exceptions import SimilarityError, ClusteringError

logger = logging.getLogger(__name__)


@dataclass
class SimilarityResult:
    """Container for similarity search results"""
    query_id: str
    target_id: str
    similarity: float
    query_molecule: Optional[Molecule] = None
    target_molecule: Optional[Molecule] = None


@dataclass
class SimilarityMatrix:
    """Container for pairwise similarity matrix"""
    molecules: List[str]  # Molecule IDs
    matrix: np.ndarray
    fingerprint_type: str
    calculation_time: float


class SimilaritySearch:
    """
    High-performance molecular similarity search engine
    Supports Tanimoto, Cosine, and custom similarity metrics
    """
    
    def __init__(
        self,
        fingerprint_type: str = "morgan",
        fingerprint_radius: int = 2,
        fingerprint_bits: int = 2048,
        similarity_metric: str = "tanimoto"
    ):
        """
        Initialize similarity search engine
        
        Args:
            fingerprint_type: Type of molecular fingerprint
            fingerprint_radius: Morgan fingerprint radius
            fingerprint_bits: Number of bits in fingerprint
            similarity_metric: Similarity metric to use
        """
        self.fingerprint_type = fingerprint_type.lower()
        self.fingerprint_radius = fingerprint_radius
        self.fingerprint_bits = fingerprint_bits
        self.similarity_metric = similarity_metric.lower()
        
        # Cached fingerprints for efficiency
        self._fingerprint_cache: Dict[str, np.ndarray] = {}
        
        # Supported metrics
        self.supported_metrics = ["tanimoto", "cosine", "dice", "jaccard"]
        
        if self.similarity_metric not in self.supported_metrics:
            raise SimilarityError(f"Unsupported similarity metric: {similarity_metric}")
    
    def calculate_fingerprint(self, molecule: Molecule) -> Optional[np.ndarray]:
        """
        Calculate molecular fingerprint
        
        Args:
            molecule: Molecule object
            
        Returns:
            Fingerprint as numpy array
        """
        # Check cache first
        cache_key = f"{molecule.mol_id}_{self.fingerprint_type}_{self.fingerprint_radius}_{self.fingerprint_bits}"
        if cache_key in self._fingerprint_cache:
            return self._fingerprint_cache[cache_key]
        
        # Calculate fingerprint
        fp = molecule.get_fingerprint(
            fp_type=self.fingerprint_type,
            radius=self.fingerprint_radius,
            n_bits=self.fingerprint_bits
        )
        
        # Cache result
        if fp is not None:
            self._fingerprint_cache[cache_key] = fp
        
        return fp
    
    def calculate_similarity(
        self,
        mol1: Molecule,
        mol2: Molecule,
        metric: Optional[str] = None
    ) -> Optional[float]:
        """
        Calculate similarity between two molecules
        
        Args:
            mol1: First molecule
            mol2: Second molecule
            metric: Similarity metric (overrides default)
            
        Returns:
            Similarity score (0-1)
        """
        fp1 = self.calculate_fingerprint(mol1)
        fp2 = self.calculate_fingerprint(mol2)
        
        if fp1 is None or fp2 is None:
            return None
        
        metric = metric or self.similarity_metric
        
        return self._calculate_fingerprint_similarity(fp1, fp2, metric)
    
    def _calculate_fingerprint_similarity(
        self,
        fp1: np.ndarray,
        fp2: np.ndarray,
        metric: str
    ) -> float:
        """Calculate similarity between fingerprints"""
        
        if metric == "tanimoto":
            return self._tanimoto_similarity(fp1, fp2)
        elif metric == "cosine":
            return cosine_similarity(fp1.reshape(1, -1), fp2.reshape(1, -1))[0, 0]
        elif metric == "dice":
            return self._dice_similarity(fp1, fp2)
        elif metric == "jaccard":
            return self._jaccard_similarity(fp1, fp2)
        else:
            raise SimilarityError(f"Unknown similarity metric: {metric}")
    
    def _tanimoto_similarity(self, fp1: np.ndarray, fp2: np.ndarray) -> float:
        """Calculate Tanimoto similarity coefficient"""
        intersection = np.sum(fp1 & fp2)
        union = np.sum(fp1 | fp2)
        
        if union == 0:
            return 0.0
        
        return intersection / union
    
    def _dice_similarity(self, fp1: np.ndarray, fp2: np.ndarray) -> float:
        """Calculate Dice similarity coefficient"""
        intersection = np.sum(fp1 & fp2)
        total = np.sum(fp1) + np.sum(fp2)
        
        if total == 0:
            return 0.0
        
        return (2 * intersection) / total
    
    def _jaccard_similarity(self, fp1: np.ndarray, fp2: np.ndarray) -> float:
        """Calculate Jaccard similarity coefficient"""
        intersection = np.sum(fp1 & fp2)
        union = np.sum(fp1) + np.sum(fp2) - intersection
        
        if union == 0:
            return 0.0
        
        return intersection / union
    
    def search_single(
        self,
        query_molecule: Molecule,
        target_collection: Union[MoleculeCollection, CompoundLibrary],
        threshold: float = 0.7,
        max_results: Optional[int] = None
    ) -> List[SimilarityResult]:
        """
        Search for similar molecules to a single query
        
        Args:
            query_molecule: Query molecule
            target_collection: Collection to search against
            threshold: Minimum similarity threshold
            max_results: Maximum number of results to return
            
        Returns:
            List of similarity results
        """
        results = []
        query_fp = self.calculate_fingerprint(query_molecule)
        
        if query_fp is None:
            logger.warning(f"Could not generate fingerprint for query {query_molecule.mol_id}")
            return results
        
        # Search through target collection
        for target_mol in target_collection:
            if target_mol.mol_id == query_molecule.mol_id:
                continue  # Skip self-similarity
            
            target_fp = self.calculate_fingerprint(target_mol)
            if target_fp is None:
                continue
            
            similarity = self._calculate_fingerprint_similarity(
                query_fp, target_fp, self.similarity_metric
            )
            
            if similarity >= threshold:
                result = SimilarityResult(
                    query_id=query_molecule.mol_id,
                    target_id=target_mol.mol_id,
                    similarity=similarity,
                    query_molecule=query_molecule,
                    target_molecule=target_mol
                )
                results.append(result)
        
        # Sort by similarity (descending)
        results.sort(key=lambda x: x.similarity, reverse=True)
        
        # Limit results if requested
        if max_results:
            results = results[:max_results]
        
        return results
    
    def search_batch(
        self,
        query_collection: Union[MoleculeCollection, CompoundLibrary],
        target_collection: Union[MoleculeCollection, CompoundLibrary],
        threshold: float = 0.7,
        max_results_per_query: Optional[int] = None
    ) -> Dict[str, List[SimilarityResult]]:
        """
        Batch similarity search for multiple queries
        
        Args:
            query_collection: Collection of query molecules
            target_collection: Collection to search against
            threshold: Minimum similarity threshold
            max_results_per_query: Maximum results per query
            
        Returns:
            Dictionary mapping query IDs to similarity results
        """
        batch_results = {}
        
        for query_mol in query_collection:
            results = self.search_single(
                query_mol, target_collection, threshold, max_results_per_query
            )
            batch_results[query_mol.mol_id] = results
        
        return batch_results
    
    def calculate_pairwise_matrix(
        self,
        molecule_collection: Union[MoleculeCollection, CompoundLibrary],
        symmetric: bool = True
    ) -> SimilarityMatrix:
        """
        Calculate pairwise similarity matrix for a collection
        
        Args:
            molecule_collection: Collection of molecules
            symmetric: Whether to assume symmetric similarity
            
        Returns:
            SimilarityMatrix object
        """
        start_time = time.time()
        
        molecules = list(molecule_collection)
        n_molecules = len(molecules)
        
        if n_molecules == 0:
            raise SimilarityError("Empty molecule collection")
        
        # Calculate all fingerprints
        fingerprints = []
        mol_ids = []
        
        for mol in molecules:
            fp = self.calculate_fingerprint(mol)
            if fp is not None:
                fingerprints.append(fp)
                mol_ids.append(mol.mol_id)
        
        if not fingerprints:
            raise SimilarityError("No valid fingerprints could be calculated")
        
        # Calculate similarity matrix
        n_valid = len(fingerprints)
        similarity_matrix = np.zeros((n_valid, n_valid))
        
        for i in range(n_valid):
            for j in range(i if symmetric else 0, n_valid):
                if i == j:
                    similarity_matrix[i, j] = 1.0
                else:
                    sim = self._calculate_fingerprint_similarity(
                        fingerprints[i], fingerprints[j], self.similarity_metric
                    )
                    similarity_matrix[i, j] = sim
                    
                    if symmetric:
                        similarity_matrix[j, i] = sim
        
        calculation_time = time.time() - start_time
        
        return SimilarityMatrix(
            molecules=mol_ids,
            matrix=similarity_matrix,
            fingerprint_type=self.fingerprint_type,
            calculation_time=calculation_time
        )
    
    def find_diverse_subset(
        self,
        molecule_collection: Union[MoleculeCollection, CompoundLibrary],
        subset_size: int,
        diversity_threshold: float = 0.8
    ) -> List[Molecule]:
        """
        Find diverse subset of molecules using MaxMin algorithm
        
        Args:
            molecule_collection: Collection to select from
            subset_size: Desired subset size
            diversity_threshold: Minimum diversity threshold
            
        Returns:
            List of diverse molecules
        """
        molecules = list(molecule_collection)
        
        if len(molecules) <= subset_size:
            return molecules
        
        # Calculate similarity matrix
        sim_matrix = self.calculate_pairwise_matrix(molecule_collection)
        
        # MaxMin algorithm
        selected_indices = []
        remaining_indices = list(range(len(molecules)))
        
        # Start with random molecule
        selected_indices.append(remaining_indices.pop(0))
        
        while len(selected_indices) < subset_size and remaining_indices:
            best_candidate = -1
            max_min_distance = -1
            
            for candidate_idx in remaining_indices:
                min_distance = float('inf')
                
                for selected_idx in selected_indices:
                    # Convert similarity to distance
                    distance = 1 - sim_matrix.matrix[candidate_idx, selected_idx]
                    min_distance = min(min_distance, distance)
                
                if min_distance > max_min_distance:
                    max_min_distance = min_distance
                    best_candidate = candidate_idx
            
            if best_candidate != -1 and max_min_distance >= (1 - diversity_threshold):
                selected_indices.append(best_candidate)
                remaining_indices.remove(best_candidate)
            else:
                break
        
        return [molecules[i] for i in selected_indices]
    
    def clear_cache(self) -> None:
        """Clear fingerprint cache"""
        self._fingerprint_cache.clear()
    
    def get_cache_stats(self) -> Dict[str, Any]:
        """Get cache statistics"""
        return {
            "cache_size": len(self._fingerprint_cache),
            "fingerprint_type": self.fingerprint_type,
            "fingerprint_bits": self.fingerprint_bits
        }


@dataclass
class ClusterResult:
    """Container for clustering results"""
    cluster_labels: np.ndarray
    n_clusters: int
    silhouette_score: Optional[float] = None
    algorithm: Optional[str] = None
    parameters: Optional[Dict[str, Any]] = None


class ClusterAnalyzer:
    """
    Molecular clustering analysis using various algorithms
    Optimized for chemical space exploration and diversity analysis
    """
    
    def __init__(
        self,
        similarity_search: Optional[SimilaritySearch] = None,
        algorithm: str = "butina"
    ):
        """
        Initialize cluster analyzer
        
        Args:
            similarity_search: SimilaritySearch instance
            algorithm: Clustering algorithm to use
        """
        self.similarity_search = similarity_search or SimilaritySearch()
        self.algorithm = algorithm.lower()
        
        # Supported algorithms
        self.supported_algorithms = ["butina", "ward", "kmeans"]
        
        if self.algorithm not in self.supported_algorithms:
            raise ClusteringError(f"Unsupported clustering algorithm: {algorithm}")
    
    def cluster_butina(
        self,
        molecule_collection: Union[MoleculeCollection, CompoundLibrary],
        distance_threshold: float = 0.6,
        min_cluster_size: int = 2
    ) -> ClusterResult:
        """
        Butina clustering algorithm for molecular data
        
        Args:
            molecule_collection: Collection to cluster
            distance_threshold: Distance threshold for clustering
            min_cluster_size: Minimum cluster size
            
        Returns:
            ClusterResult object
        """
        try:
            molecules = list(molecule_collection)
            n_molecules = len(molecules)
            
            if n_molecules < 2:
                raise ClusteringError("Need at least 2 molecules for clustering")
            
            # Calculate fingerprints
            fps = []
            mol_indices = []
            
            for i, mol in enumerate(molecules):
                fp = self.similarity_search.calculate_fingerprint(mol)
                if fp is not None:
                    fps.append(fp)
                    mol_indices.append(i)
            
            if len(fps) < 2:
                raise ClusteringError("Need at least 2 valid fingerprints for clustering")
            
            # Convert to RDKit bit vectors for Butina clustering
            from rdkit import DataStructs
            bit_vectors = []
            
            for fp in fps:
                bv = DataStructs.ExplicitBitVect(len(fp))
                for i, bit in enumerate(fp):
                    if bit:
                        bv.SetBit(i)
                bit_vectors.append(bv)
            
            # Calculate distance matrix
            distances = []
            for i in range(len(bit_vectors)):
                for j in range(i + 1, len(bit_vectors)):
                    # Tanimoto distance = 1 - Tanimoto similarity
                    tanimoto_sim = DataStructs.TanimotoSimilarity(bit_vectors[i], bit_vectors[j])
                    distance = 1.0 - tanimoto_sim
                    distances.append(distance)
            
            # Butina clustering
            cluster_indices = Butina.ClusterData(
                distances, len(bit_vectors), distance_threshold, isDistData=True
            )
            
            # Create cluster labels array
            cluster_labels = np.full(n_molecules, -1)  # -1 for unclustered
            
            for cluster_id, cluster in enumerate(cluster_indices):
                if len(cluster) >= min_cluster_size:
                    for mol_idx in cluster:
                        original_idx = mol_indices[mol_idx]
                        cluster_labels[original_idx] = cluster_id
            
            n_clusters = len([c for c in cluster_indices if len(c) >= min_cluster_size])
            
            return ClusterResult(
                cluster_labels=cluster_labels,
                n_clusters=n_clusters,
                algorithm="butina",
                parameters={
                    "distance_threshold": distance_threshold,
                    "min_cluster_size": min_cluster_size
                }
            )
            
        except Exception as e:
            raise ClusteringError(f"Butina clustering failed: {e}")
    
    def cluster_hierarchical(
        self,
        molecule_collection: Union[MoleculeCollection, CompoundLibrary],
        n_clusters: Optional[int] = None,
        distance_threshold: Optional[float] = None,
        linkage_method: str = "ward"
    ) -> ClusterResult:
        """
        Hierarchical clustering using Ward linkage
        
        Args:
            molecule_collection: Collection to cluster
            n_clusters: Number of clusters (if specified)
            distance_threshold: Distance threshold (alternative to n_clusters)
            linkage_method: Linkage method for hierarchy
            
        Returns:
            ClusterResult object
        """
        try:
            # Calculate similarity matrix
            sim_matrix = self.similarity_search.calculate_pairwise_matrix(molecule_collection)
            
            # Convert similarity to distance matrix
            distance_matrix = 1 - sim_matrix.matrix
            
            # Hierarchical clustering
            condensed_distances = squareform(distance_matrix)
            linkage_matrix = linkage(condensed_distances, method=linkage_method)
            
            # Get cluster labels
            if n_clusters is not None:
                cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust') - 1
                actual_n_clusters = n_clusters
            elif distance_threshold is not None:
                cluster_labels = fcluster(linkage_matrix, distance_threshold, criterion='distance') - 1
                actual_n_clusters = len(np.unique(cluster_labels))
            else:
                raise ClusteringError("Must specify either n_clusters or distance_threshold")
            
            return ClusterResult(
                cluster_labels=cluster_labels,
                n_clusters=actual_n_clusters,
                algorithm=f"hierarchical_{linkage_method}",
                parameters={
                    "n_clusters": n_clusters,
                    "distance_threshold": distance_threshold,
                    "linkage_method": linkage_method
                }
            )
            
        except Exception as e:
            raise ClusteringError(f"Hierarchical clustering failed: {e}")
    
    def cluster_collection(
        self,
        molecule_collection: Union[MoleculeCollection, CompoundLibrary],
        **kwargs
    ) -> ClusterResult:
        """
        Cluster molecule collection using configured algorithm
        
        Args:
            molecule_collection: Collection to cluster
            **kwargs: Algorithm-specific parameters
            
        Returns:
            ClusterResult object
        """
        if self.algorithm == "butina":
            return self.cluster_butina(molecule_collection, **kwargs)
        elif self.algorithm in ["ward", "hierarchical"]:
            return self.cluster_hierarchical(molecule_collection, **kwargs)
        else:
            raise ClusteringError(f"Unsupported algorithm: {self.algorithm}")
    
    def analyze_clusters(
        self,
        molecule_collection: Union[MoleculeCollection, CompoundLibrary],
        cluster_result: ClusterResult
    ) -> Dict[str, Any]:
        """
        Analyze clustering results and provide statistics
        
        Args:
            molecule_collection: Original molecule collection
            cluster_result: Clustering results
            
        Returns:
            Dictionary of cluster analysis results
        """
        molecules = list(molecule_collection)
        cluster_labels = cluster_result.cluster_labels
        
        # Basic statistics
        unique_clusters = np.unique(cluster_labels)
        n_clusters = len(unique_clusters[unique_clusters >= 0])  # Exclude -1 (noise)
        n_noise = np.sum(cluster_labels == -1)
        
        # Cluster sizes
        cluster_sizes = {}
        for cluster_id in unique_clusters:
            if cluster_id >= 0:
                cluster_sizes[int(cluster_id)] = int(np.sum(cluster_labels == cluster_id))
        
        # Representative molecules (closest to cluster center)
        representatives = {}
        if hasattr(molecule_collection, 'compounds'):
            # CompoundLibrary
            for cluster_id in range(n_clusters):
                cluster_indices = np.where(cluster_labels == cluster_id)[0]
                if len(cluster_indices) > 0:
                    # For now, just take the first compound as representative
                    representatives[cluster_id] = molecules[cluster_indices[0]].mol_id
        
        analysis = {
            "n_clusters": n_clusters,
            "n_noise_points": n_noise,
            "cluster_sizes": cluster_sizes,
            "representatives": representatives,
            "clustering_efficiency": (len(molecules) - n_noise) / len(molecules),
            "average_cluster_size": np.mean(list(cluster_sizes.values())) if cluster_sizes else 0,
            "algorithm": cluster_result.algorithm,
            "parameters": cluster_result.parameters
        }
        
        return analysis
    
    def visualize_clusters(
        self,
        molecule_collection: Union[MoleculeCollection, CompoundLibrary],
        cluster_result: ClusterResult,
        method: str = "pca"
    ) -> Dict[str, Any]:
        """
        Generate data for cluster visualization
        
        Args:
            molecule_collection: Original molecule collection
            cluster_result: Clustering results
            method: Dimensionality reduction method
            
        Returns:
            Visualization data dictionary
        """
        # Calculate fingerprint matrix
        fp_matrix = molecule_collection.get_fingerprint_matrix(
            fp_type=self.similarity_search.fingerprint_type,
            radius=self.similarity_search.fingerprint_radius,
            n_bits=self.similarity_search.fingerprint_bits
        )
        
        if fp_matrix is None:
            raise ClusteringError("Could not generate fingerprint matrix for visualization")
        
        # Dimensionality reduction
        if method.lower() == "pca":
            from sklearn.decomposition import PCA
            reducer = PCA(n_components=2)
            coordinates = reducer.fit_transform(fp_matrix)
        elif method.lower() == "tsne":
            from sklearn.manifold import TSNE
            reducer = TSNE(n_components=2, random_state=42)
            coordinates = reducer.fit_transform(fp_matrix)
        else:
            raise ClusteringError(f"Unsupported visualization method: {method}")
        
        return {
            "coordinates": coordinates.tolist(),
            "cluster_labels": cluster_result.cluster_labels.tolist(),
            "molecule_ids": [mol.mol_id for mol in molecule_collection],
            "method": method
        }