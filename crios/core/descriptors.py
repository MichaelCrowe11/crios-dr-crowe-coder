"""
CriOS Molecular Descriptors Calculator
Comprehensive molecular property calculation and analysis
"""

import logging
from typing import Dict, List, Optional, Union, Any
from dataclasses import dataclass, field
import warnings

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator

from .molecule import Molecule, MoleculeCollection
from .compound import Compound, CompoundLibrary
from ..exceptions import CalculationError

logger = logging.getLogger(__name__)

# Suppress RDKit warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, module="rdkit")


@dataclass
class MolecularDescriptors:
    """Container for comprehensive molecular descriptors"""
    
    # Basic properties
    molecular_weight: Optional[float] = None
    heavy_atom_count: Optional[int] = None
    
    # Lipinski descriptors
    logp: Optional[float] = None
    hbd: Optional[int] = None  # Hydrogen bond donors
    hba: Optional[int] = None  # Hydrogen bond acceptors
    tpsa: Optional[float] = None  # Topological polar surface area
    
    # Structural descriptors
    rotatable_bonds: Optional[int] = None
    aromatic_rings: Optional[int] = None
    aliphatic_rings: Optional[int] = None
    saturated_rings: Optional[int] = None
    
    # Connectivity descriptors
    chi0v: Optional[float] = None  # Valence connectivity index
    chi1v: Optional[float] = None
    chi2v: Optional[float] = None
    chi3v: Optional[float] = None
    chi4v: Optional[float] = None
    
    # Kappa shape indices
    kappa1: Optional[float] = None
    kappa2: Optional[float] = None
    kappa3: Optional[float] = None
    
    # Pharmacophore descriptors
    num_heteroatoms: Optional[int] = None
    num_bridgehead_atoms: Optional[int] = None
    num_spiro_atoms: Optional[int] = None
    
    # Complexity measures
    bertz_complexity: Optional[float] = None
    mol_mr: Optional[float] = None  # Molar refractivity
    
    # Fragment-based descriptors
    num_fragments: Optional[int] = None
    fraction_csp3: Optional[float] = None
    
    # Additional descriptors
    custom_descriptors: Dict[str, float] = field(default_factory=dict)


class DescriptorCalculator:
    """
    High-performance molecular descriptor calculator
    Supports RDKit descriptors plus custom Crowe Discovery Framework descriptors
    """
    
    def __init__(self, include_3d: bool = False, custom_descriptors: Optional[List[str]] = None):
        """
        Initialize descriptor calculator
        
        Args:
            include_3d: Whether to include 3D descriptors (requires 3D conformers)
            custom_descriptors: List of custom descriptor names to calculate
        """
        self.include_3d = include_3d
        self.custom_descriptors = custom_descriptors or []
        
        # Standard RDKit descriptors to calculate
        self.standard_descriptors = [
            'MolWt', 'HeavyAtomCount', 'MolLogP', 'NumHDonors', 'NumHAcceptors',
            'TPSA', 'NumRotatableBonds', 'NumAromaticRings', 'NumAliphaticRings',
            'NumSaturatedRings', 'Chi0v', 'Chi1v', 'Chi2v', 'Chi3v', 'Chi4v',
            'Kappa1', 'Kappa2', 'Kappa3', 'NumHeteroatoms', 'BertzCT',
            'MolMR', 'FractionCsp3', 'BalabanJ', 'HallKierAlpha'
        ]
        
        # Initialize RDKit descriptor calculator
        self.rdkit_calculator = MolecularDescriptorCalculator(self.standard_descriptors)
    
    def calculate_single(self, molecule: Molecule) -> Optional[MolecularDescriptors]:
        """
        Calculate all descriptors for a single molecule
        
        Args:
            molecule: Molecule to calculate descriptors for
            
        Returns:
            MolecularDescriptors object or None if calculation fails
        """
        if not molecule.is_valid():
            logger.warning(f"Invalid molecule {molecule.mol_id}, skipping descriptor calculation")
            return None
        
        try:
            mol = molecule.rdkit_mol
            descriptors = MolecularDescriptors()
            
            # Calculate standard RDKit descriptors
            desc_values = self.rdkit_calculator.CalcDescriptors(mol)
            desc_dict = dict(zip(self.standard_descriptors, desc_values))
            
            # Map to our descriptor container
            descriptors.molecular_weight = desc_dict.get('MolWt')
            descriptors.heavy_atom_count = desc_dict.get('HeavyAtomCount')
            descriptors.logp = desc_dict.get('MolLogP')
            descriptors.hbd = desc_dict.get('NumHDonors')
            descriptors.hba = desc_dict.get('NumHAcceptors')
            descriptors.tpsa = desc_dict.get('TPSA')
            descriptors.rotatable_bonds = desc_dict.get('NumRotatableBonds')
            descriptors.aromatic_rings = desc_dict.get('NumAromaticRings')
            descriptors.aliphatic_rings = desc_dict.get('NumAliphaticRings')
            descriptors.saturated_rings = desc_dict.get('NumSaturatedRings')
            descriptors.chi0v = desc_dict.get('Chi0v')
            descriptors.chi1v = desc_dict.get('Chi1v')
            descriptors.chi2v = desc_dict.get('Chi2v')
            descriptors.chi3v = desc_dict.get('Chi3v')
            descriptors.chi4v = desc_dict.get('Chi4v')
            descriptors.kappa1 = desc_dict.get('Kappa1')
            descriptors.kappa2 = desc_dict.get('Kappa2')
            descriptors.kappa3 = desc_dict.get('Kappa3')
            descriptors.bertz_complexity = desc_dict.get('BertzCT')
            descriptors.mol_mr = desc_dict.get('MolMR')
            descriptors.fraction_csp3 = desc_dict.get('FractionCsp3')
            
            # Additional calculated descriptors
            descriptors.num_heteroatoms = desc_dict.get('NumHeteroatoms')
            descriptors.num_bridgehead_atoms = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
            descriptors.num_spiro_atoms = rdMolDescriptors.CalcNumSpiroAtoms(mol)
            descriptors.num_fragments = len(Chem.GetMolFrags(mol))
            
            # Custom Crowe Discovery descriptors
            self._calculate_custom_descriptors(mol, descriptors)
            
            return descriptors
            
        except Exception as e:
            logger.error(f"Failed to calculate descriptors for {molecule.mol_id}: {e}")
            return None
    
    def _calculate_custom_descriptors(self, mol: Chem.Mol, descriptors: MolecularDescriptors) -> None:
        """Calculate custom descriptors specific to Crowe Discovery Framework"""
        try:
            # Drug-like complexity score
            descriptors.custom_descriptors['drug_complexity'] = self._calculate_drug_complexity(mol)
            
            # Natural product likeness (simplified)
            descriptors.custom_descriptors['natural_product_score'] = self._calculate_np_score(mol)
            
            # Synthetic accessibility approximation
            descriptors.custom_descriptors['synthetic_accessibility'] = self._calculate_sa_score(mol)
            
            # Crowe scaffold diversity metric
            descriptors.custom_descriptors['scaffold_diversity'] = self._calculate_scaffold_diversity(mol)
            
        except Exception as e:
            logger.warning(f"Failed to calculate custom descriptors: {e}")
    
    def _calculate_drug_complexity(self, mol: Chem.Mol) -> float:
        """Calculate drug-like complexity score (0-1)"""
        try:
            # Combination of molecular weight, flexibility, and structural complexity
            mw = Descriptors.MolWt(mol)
            flex = Descriptors.NumRotatableBonds(mol)
            rings = Descriptors.RingCount(mol)
            hetero = rdMolDescriptors.CalcNumHeteroatoms(mol)
            
            # Normalize components
            mw_norm = min(1.0, mw / 500.0)  # Normalize around typical drug MW
            flex_norm = min(1.0, flex / 10.0)  # Normalize around typical flexibility
            ring_norm = min(1.0, rings / 4.0)  # Normalize around typical ring count
            hetero_norm = min(1.0, hetero / 8.0)  # Normalize around typical heteroatom count
            
            # Weighted combination
            complexity = 0.3 * mw_norm + 0.3 * flex_norm + 0.2 * ring_norm + 0.2 * hetero_norm
            
            return min(1.0, complexity)
            
        except Exception:
            return 0.5  # Default moderate complexity
    
    def _calculate_np_score(self, mol: Chem.Mol) -> float:
        """Calculate natural product likeness score (0-1)"""
        try:
            # Simplified natural product score based on structural features
            # Higher scores indicate more natural product-like
            
            score = 0.5  # Base score
            
            # Presence of oxygen increases NP-likeness
            oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
            score += min(0.2, oxygen_count * 0.05)
            
            # Multiple rings increase NP-likeness
            ring_count = Descriptors.RingCount(mol)
            score += min(0.2, ring_count * 0.05)
            
            # Stereogenic centers increase NP-likeness
            stereo_centers = len(Chem.FindMolChiralCenters(mol))
            score += min(0.1, stereo_centers * 0.02)
            
            # Too many aromatic rings decrease NP-likeness
            aromatic_rings = Descriptors.NumAromaticRings(mol)
            if aromatic_rings > 3:
                score -= 0.1
            
            return max(0.0, min(1.0, score))
            
        except Exception:
            return 0.5
    
    def _calculate_sa_score(self, mol: Chem.Mol) -> float:
        """Calculate synthetic accessibility score approximation (0-1, higher = more accessible)"""
        try:
            # Simplified synthetic accessibility based on structural features
            # This is a rough approximation - real SA scores require trained models
            
            score = 0.8  # Start with good accessibility
            
            # Penalize high molecular weight
            mw = Descriptors.MolWt(mol)
            if mw > 500:
                score -= 0.2
            elif mw > 400:
                score -= 0.1
            
            # Penalize high flexibility
            rotatable_bonds = Descriptors.NumRotatableBonds(mol)
            if rotatable_bonds > 10:
                score -= 0.3
            elif rotatable_bonds > 7:
                score -= 0.1
            
            # Penalize complex ring systems
            ring_count = Descriptors.RingCount(mol)
            if ring_count > 4:
                score -= 0.2
            elif ring_count > 2:
                score -= 0.1
            
            # Penalize many stereocenters
            stereo_centers = len(Chem.FindMolChiralCenters(mol))
            if stereo_centers > 4:
                score -= 0.2
            elif stereo_centers > 2:
                score -= 0.1
            
            return max(0.0, min(1.0, score))
            
        except Exception:
            return 0.5
    
    def _calculate_scaffold_diversity(self, mol: Chem.Mol) -> float:
        """Calculate scaffold diversity metric"""
        try:
            # Measure of how diverse the molecular scaffold is
            # Based on ring systems and connectivity
            
            ring_info = mol.GetRingInfo()
            num_rings = ring_info.NumRings()
            
            if num_rings == 0:
                return 0.2  # No rings = low scaffold diversity
            
            # Analyze ring sizes
            ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
            unique_sizes = len(set(ring_sizes))
            
            # Analyze ring connectivity
            fused_rings = sum(1 for ring1 in ring_info.AtomRings() 
                             for ring2 in ring_info.AtomRings()
                             if ring1 != ring2 and len(set(ring1) & set(ring2)) >= 2)
            
            # Calculate diversity score
            size_diversity = unique_sizes / max(1, num_rings)
            fusion_factor = min(1.0, fused_rings / max(1, num_rings))
            
            diversity = 0.5 * size_diversity + 0.3 * fusion_factor + 0.2 * min(1.0, num_rings / 3.0)
            
            return min(1.0, diversity)
            
        except Exception:
            return 0.5
    
    def calculate_batch(
        self,
        molecule_collection: Union[MoleculeCollection, CompoundLibrary],
        parallel: bool = True,
        chunk_size: int = 100
    ) -> Dict[str, MolecularDescriptors]:
        """
        Calculate descriptors for a collection of molecules
        
        Args:
            molecule_collection: Collection to process
            parallel: Whether to use parallel processing
            chunk_size: Size of processing chunks
            
        Returns:
            Dictionary mapping molecule IDs to descriptors
        """
        results = {}
        molecules = list(molecule_collection)
        
        if parallel:
            return self._calculate_batch_parallel(molecules, chunk_size)
        else:
            return self._calculate_batch_sequential(molecules)
    
    def _calculate_batch_sequential(self, molecules: List[Molecule]) -> Dict[str, MolecularDescriptors]:
        """Sequential batch calculation"""
        results = {}
        
        for molecule in molecules:
            descriptors = self.calculate_single(molecule)
            if descriptors is not None:
                results[molecule.mol_id] = descriptors
        
        return results
    
    def _calculate_batch_parallel(
        self, 
        molecules: List[Molecule], 
        chunk_size: int
    ) -> Dict[str, MolecularDescriptors]:
        """Parallel batch calculation using joblib"""
        try:
            from joblib import Parallel, delayed
            import multiprocessing as mp
            
            n_jobs = min(mp.cpu_count(), 8)  # Limit to 8 cores max
            
            # Split molecules into chunks
            chunks = [molecules[i:i + chunk_size] for i in range(0, len(molecules), chunk_size)]
            
            # Process chunks in parallel
            chunk_results = Parallel(n_jobs=n_jobs)(
                delayed(self._calculate_batch_sequential)(chunk) for chunk in chunks
            )
            
            # Combine results
            results = {}
            for chunk_result in chunk_results:
                results.update(chunk_result)
            
            return results
            
        except ImportError:
            logger.warning("joblib not available, falling back to sequential processing")
            return self._calculate_batch_sequential(molecules)
    
    def to_dataframe(
        self,
        descriptor_dict: Dict[str, MolecularDescriptors],
        include_custom: bool = True
    ) -> pd.DataFrame:
        """
        Convert descriptor results to pandas DataFrame
        
        Args:
            descriptor_dict: Dictionary of molecule ID to descriptors
            include_custom: Whether to include custom descriptors
            
        Returns:
            DataFrame with descriptors as columns
        """
        data = []
        
        for mol_id, descriptors in descriptor_dict.items():
            row = {'molecule_id': mol_id}
            
            # Add standard descriptors
            descriptor_fields = [
                'molecular_weight', 'heavy_atom_count', 'logp', 'hbd', 'hba', 'tpsa',
                'rotatable_bonds', 'aromatic_rings', 'aliphatic_rings', 'saturated_rings',
                'chi0v', 'chi1v', 'chi2v', 'chi3v', 'chi4v', 'kappa1', 'kappa2', 'kappa3',
                'num_heteroatoms', 'num_bridgehead_atoms', 'num_spiro_atoms',
                'bertz_complexity', 'mol_mr', 'num_fragments', 'fraction_csp3'
            ]
            
            for field in descriptor_fields:
                row[field] = getattr(descriptors, field, None)
            
            # Add custom descriptors
            if include_custom and descriptors.custom_descriptors:
                for name, value in descriptors.custom_descriptors.items():
                    row[f'custom_{name}'] = value
            
            data.append(row)
        
        return pd.DataFrame(data)
    
    def calculate_lipinski_violations(self, descriptors: MolecularDescriptors) -> List[str]:
        """Calculate Lipinski Rule of Five violations"""
        violations = []
        
        if descriptors.molecular_weight and descriptors.molecular_weight > 500:
            violations.append("Molecular weight > 500 Da")
        
        if descriptors.logp and descriptors.logp > 5:
            violations.append("LogP > 5")
        
        if descriptors.hbd and descriptors.hbd > 5:
            violations.append("Hydrogen bond donors > 5")
        
        if descriptors.hba and descriptors.hba > 10:
            violations.append("Hydrogen bond acceptors > 10")
        
        return violations
    
    def calculate_veber_violations(self, descriptors: MolecularDescriptors) -> List[str]:
        """Calculate Veber rule violations"""
        violations = []
        
        if descriptors.rotatable_bonds and descriptors.rotatable_bonds > 10:
            violations.append("Rotatable bonds > 10")
        
        if descriptors.tpsa and descriptors.tpsa > 140:
            violations.append("TPSA > 140 Ų")
        
        return violations
    
    def get_descriptor_statistics(
        self,
        descriptor_dict: Dict[str, MolecularDescriptors]
    ) -> Dict[str, Dict[str, float]]:
        """Calculate statistics for descriptor values across collection"""
        df = self.to_dataframe(descriptor_dict)
        numeric_columns = df.select_dtypes(include=[np.number]).columns
        
        stats = {}
        for column in numeric_columns:
            if column != 'molecule_id':
                values = df[column].dropna()
                if len(values) > 0:
                    stats[column] = {
                        'mean': float(values.mean()),
                        'std': float(values.std()),
                        'min': float(values.min()),
                        'max': float(values.max()),
                        'median': float(values.median()),
                        'count': len(values)
                    }
        
        return stats