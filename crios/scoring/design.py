"""
CriOS Scoring and Design - Integrated Client Functions
Port and integration of client RDKit functions with Crowe scoring methodology
"""

import logging
import math
from typing import Dict, List, Optional, Union, Any, Tuple
import warnings

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, Lipinski
from rdkit.Chem.rdMolDescriptors import CalcNumHeteroatoms, CalcNumRotatableBonds
from rdkit.Chem.Fragments import fr_benzene, fr_furan, fr_imidazole
from rdkit.Chem.Scaffolds import MurckoScaffold

from ..chem.mol import Molecule
from ..core.exceptions import CriosError

logger = logging.getLogger(__name__)


class ScoringError(CriosError):
    """Scoring-specific errors"""
    pass


# ============================================================================
# CLIENT FUNCTIONS - INTEGRATED VERBATIM WITH ENHANCEMENTS
# ============================================================================

def _passes_drug_like_filters(
    smiles: str,
    config: Optional[Dict[str, Any]] = None
) -> bool:
    """
    Check if molecule passes drug-like filters (Lipinski + extended)
    
    Args:
        smiles: SMILES string
        config: Filter configuration (optional)
    
    Returns:
        True if passes all filters
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        # Default thresholds
        thresholds = {
            "mw_max": 500.0,
            "mw_min": 150.0,
            "logp_max": 5.0,
            "logp_min": -2.0,
            "hbd_max": 5,
            "hba_max": 10,
            "tpsa_max": 140.0,
            "tpsa_min": 20.0,
            "rotbonds_max": 10,
            "rings_max": 4
        }
        
        if config:
            thresholds.update(config)
        
        # Calculate properties
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        rotbonds = CalcNumRotatableBonds(mol)
        rings = Descriptors.RingCount(mol)
        
        # Apply filters
        if not (thresholds["mw_min"] <= mw <= thresholds["mw_max"]):
            return False
        if not (thresholds["logp_min"] <= logp <= thresholds["logp_max"]):
            return False
        if hbd > thresholds["hbd_max"]:
            return False
        if hba > thresholds["hba_max"]:
            return False
        if not (thresholds["tpsa_min"] <= tpsa <= thresholds["tpsa_max"]):
            return False
        if rotbonds > thresholds["rotbonds_max"]:
            return False
        if rings > thresholds["rings_max"]:
            return False
        
        return True
        
    except Exception as e:
        logger.warning(f"Drug-like filter check failed for {smiles}: {e}")
        return False


def _calculate_design_scores(
    smiles: str,
    target_class: str = "GPCR",
    optimize_for: List[str] = None
) -> Dict[str, float]:
    """
    Calculate comprehensive design scores for a molecule
    
    Args:
        smiles: SMILES string
        target_class: Target protein class
        optimize_for: List of aspects to optimize
    
    Returns:
        Dictionary of normalized scores (0-1)
    """
    if optimize_for is None:
        optimize_for = ["potency", "selectivity", "admet"]
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {aspect: 0.0 for aspect in optimize_for}
        
        scores = {}
        
        if "potency" in optimize_for:
            scores["potency"] = _predict_target_potency(mol, target_class)
        
        if "selectivity" in optimize_for:
            scores["selectivity"] = _predict_selectivity(mol, target_class)
        
        if "admet" in optimize_for:
            admet_props = _predict_admet_properties(mol)
            # Composite ADMET score
            admet_components = [
                admet_props.get("oral_absorption", 0.5),
                admet_props.get("bbb_permeation", 0.5),
                admet_props.get("metabolic_stability", 0.5),
                admet_props.get("clearance", 0.5),
                admet_props.get("safety", 0.5)
            ]
            scores["admet"] = np.mean(admet_components)
        
        if "synthesis" in optimize_for:
            scores["synthesis"] = _assess_synthetic_accessibility(mol)
        
        if "novelty" in optimize_for:
            scores["novelty"] = _assess_structural_novelty(mol)
        
        # Ensure all values are in [0,1]
        for key, value in scores.items():
            scores[key] = np.clip(value, 0.0, 1.0)
        
        return scores
        
    except Exception as e:
        logger.error(f"Design score calculation failed for {smiles}: {e}")
        return {aspect: 0.0 for aspect in optimize_for}


def _predict_target_potency(
    mol: Chem.Mol,
    target_class: str = "GPCR"
) -> float:
    """
    Predict target potency based on structural features
    
    Args:
        mol: RDKit molecule
        target_class: Target protein class
    
    Returns:
        Potency score (0-1, higher is better)
    """
    try:
        features = _extract_molecular_features(mol)
        
        # Target class-specific heuristics
        target_class = target_class.upper()
        
        if target_class == "GPCR":
            score = _predict_gpcr_potency(features)
        elif target_class == "KINASE":
            score = _predict_kinase_potency(features)
        elif target_class == "PROTEASE":
            score = _predict_protease_potency(features)
        else:
            # Generic potency heuristic
            score = _predict_generic_potency(features)
        
        return np.clip(score, 0.0, 1.0)
        
    except Exception as e:
        logger.warning(f"Potency prediction failed: {e}")
        return 0.5


def _predict_gpcr_potency(features: Dict[str, Any]) -> float:
    """GPCR-specific potency prediction"""
    score = 0.5  # Base score
    
    # Favorable features for GPCR activity
    if features.get("basic_nitrogen", 0) > 0:
        score += 0.2
    
    if 300 <= features.get("molecular_weight", 0) <= 450:
        score += 0.1
    
    if 2.0 <= features.get("logp", 0) <= 4.0:
        score += 0.1
    
    if features.get("aromatic_rings", 0) >= 1:
        score += 0.1
    
    return min(1.0, score)


def _predict_kinase_potency(features: Dict[str, Any]) -> float:
    """Kinase-specific potency prediction"""
    score = 0.5  # Base score
    
    # Favorable features for kinase activity
    if features.get("hbd", 0) >= 2:
        score += 0.15
    
    if features.get("hba", 0) >= 4:
        score += 0.15
    
    if 250 <= features.get("molecular_weight", 0) <= 400:
        score += 0.1
    
    if features.get("planar_rings", 0) >= 2:
        score += 0.1
    
    return min(1.0, score)


def _predict_protease_potency(features: Dict[str, Any]) -> float:
    """Protease-specific potency prediction"""
    score = 0.5  # Base score
    
    # Favorable features for protease activity
    if features.get("peptide_like", False):
        score += 0.2
    
    if features.get("rotatable_bonds", 0) >= 5:
        score += 0.1
    
    if 400 <= features.get("molecular_weight", 0) <= 600:
        score += 0.1
    
    if features.get("hydroxyl_groups", 0) >= 1:
        score += 0.1
    
    return min(1.0, score)


def _predict_generic_potency(features: Dict[str, Any]) -> float:
    """Generic potency prediction"""
    score = 0.5
    
    # General drug-like features
    if 200 <= features.get("molecular_weight", 0) <= 500:
        score += 0.1
    
    if 1.0 <= features.get("logp", 0) <= 4.0:
        score += 0.1
    
    if features.get("rings", 0) >= 1:
        score += 0.1
    
    return min(1.0, score)


def _predict_selectivity(
    mol: Chem.Mol,
    target_class: str = "GPCR"
) -> float:
    """
    Predict selectivity profile
    
    Args:
        mol: RDKit molecule
        target_class: Target protein class
    
    Returns:
        Selectivity score (0-1, higher is better)
    """
    try:
        features = _extract_molecular_features(mol)
        
        # Selectivity often comes from structural specificity
        score = 0.5
        
        # Moderate complexity suggests selectivity
        complexity = features.get("structural_complexity", 0.5)
        if 0.4 <= complexity <= 0.7:
            score += 0.2
        
        # Specific binding motifs
        if features.get("specific_binding_groups", 0) > 0:
            score += 0.15
        
        # Moderate molecular weight
        mw = features.get("molecular_weight", 300)
        if 300 <= mw <= 500:
            score += 0.1
        
        # Avoid promiscuous substructures
        if features.get("promiscuous_substructures", 0) == 0:
            score += 0.15
        
        return np.clip(score, 0.0, 1.0)
        
    except Exception as e:
        logger.warning(f"Selectivity prediction failed: {e}")
        return 0.5


def _predict_admet_properties(mol: Chem.Mol) -> Dict[str, float]:
    """
    Predict ADMET properties
    
    Args:
        mol: RDKit molecule
    
    Returns:
        Dictionary of ADMET scores (0-1)
    """
    try:
        admet = {}
        
        admet["oral_absorption"] = _predict_oral_absorption(mol)
        admet["bbb_permeation"] = _predict_tissue_distribution(mol, "bbb")
        admet["metabolic_stability"] = _predict_metabolic_stability(mol)
        admet["clearance"] = _predict_clearance(mol)
        admet["safety"] = _predict_safety_profile(mol)
        
        return admet
        
    except Exception as e:
        logger.warning(f"ADMET prediction failed: {e}")
        return {
            "oral_absorption": 0.5,
            "bbb_permeation": 0.5,
            "metabolic_stability": 0.5,
            "clearance": 0.5,
            "safety": 0.5
        }


def _predict_oral_absorption(mol: Chem.Mol) -> float:
    """
    Predict oral absorption based on physicochemical properties
    
    Args:
        mol: RDKit molecule
    
    Returns:
        Oral absorption score (0-1)
    """
    try:
        # Key properties for oral absorption
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Lipinski.NumHDonors(mol)
        
        score = 1.0
        
        # Molecular weight penalty
        if mw > 500:
            score -= 0.3
        elif mw > 400:
            score -= 0.1
        
        # LogP considerations
        if logp < 0 or logp > 5:
            score -= 0.2
        elif 1 <= logp <= 3:
            score += 0.1
        
        # TPSA penalty
        if tpsa > 140:
            score -= 0.3
        elif tpsa > 100:
            score -= 0.1
        
        # HBD penalty
        if hbd > 5:
            score -= 0.2
        
        return np.clip(score, 0.0, 1.0)
        
    except Exception:
        return 0.5


def _predict_tissue_distribution(
    mol: Chem.Mol,
    tissue: str = "bbb"
) -> float:
    """
    Predict tissue distribution (BBB permeation)
    
    Args:
        mol: RDKit molecule
        tissue: Tissue type ("bbb" for blood-brain barrier)
    
    Returns:
        Distribution score (0-1)
    """
    try:
        if tissue.lower() == "bbb":
            # BBB permeation heuristics
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            hbd = Lipinski.NumHDonors(mol)
            
            score = 1.0
            
            # MW cutoff for BBB
            if mw > 450:
                score -= 0.4
            
            # LogP window for BBB
            if not (1.0 <= logp <= 4.0):
                score -= 0.3
            
            # TPSA cutoff
            if tpsa > 90:
                score -= 0.4
            
            # HBD limit
            if hbd > 3:
                score -= 0.2
            
            return np.clip(score, 0.0, 1.0)
        
        # Generic tissue distribution
        return 0.5
        
    except Exception:
        return 0.5


def _predict_metabolic_stability(mol: Chem.Mol) -> float:
    """
    Predict metabolic stability
    
    Args:
        mol: RDKit molecule
    
    Returns:
        Metabolic stability score (0-1)
    """
    try:
        score = 0.8  # Start optimistic
        
        # Count potential metabolic liability sites
        # Simplified - real implementation would use more sophisticated SMARTS
        
        # Aromatic hydroxylation sites
        aromatic_c = len(mol.GetSubstructMatches(Chem.MolFromSmarts("c")))
        if aromatic_c > 6:
            score -= 0.1
        
        # Aliphatic oxidation sites
        aliphatic_c = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CH3]")))
        if aliphatic_c > 3:
            score -= 0.1
        
        # N-dealkylation sites
        n_dealk = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[N]([CH3])")))
        if n_dealk > 1:
            score -= 0.1
        
        # O-dealkylation sites
        o_dealk = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[O][CH3]")))
        if o_dealk > 2:
            score -= 0.1
        
        return np.clip(score, 0.0, 1.0)
        
    except Exception:
        return 0.5


def _predict_clearance(mol: Chem.Mol) -> float:
    """
    Predict clearance (lower clearance = higher score)
    
    Args:
        mol: RDKit molecule
    
    Returns:
        Clearance score (0-1, higher = lower clearance)
    """
    try:
        # Factors affecting clearance
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        flexibility = CalcNumRotatableBonds(mol)
        
        score = 0.7
        
        # Higher MW generally means lower clearance
        if mw > 400:
            score += 0.1
        
        # Very lipophilic compounds may have high clearance
        if logp > 4:
            score -= 0.2
        elif logp < 1:
            score += 0.1
        
        # Flexibility affects clearance
        if flexibility > 8:
            score -= 0.1
        
        return np.clip(score, 0.0, 1.0)
        
    except Exception:
        return 0.5


def _predict_safety_profile(mol: Chem.Mol) -> float:
    """
    Predict safety profile
    
    Args:
        mol: RDKit molecule
    
    Returns:
        Safety score (0-1, higher is safer)
    """
    try:
        score = 0.8  # Start optimistic
        
        # Structural alerts (simplified)
        
        # Aromatic nitro groups
        nitro_aromatic = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[c][N+](=O)[O-]")))
        score -= nitro_aromatic * 0.2
        
        # Aromatic amines
        aromatic_amine = len(mol.GetSubstructMatches(Chem.MolFromSmarts("c[NH2]")))
        score -= aromatic_amine * 0.1
        
        # Epoxides
        epoxides = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C1OC1")))
        score -= epoxides * 0.3
        
        # Aldehydes
        aldehydes = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CH]=O")))
        score -= aldehydes * 0.1
        
        # Very high LogP can indicate toxicity
        logp = Crippen.MolLogP(mol)
        if logp > 5:
            score -= 0.2
        
        return np.clip(score, 0.0, 1.0)
        
    except Exception:
        return 0.5


def _extract_molecular_features(mol: Chem.Mol) -> Dict[str, Any]:
    """
    Extract comprehensive molecular features
    
    Args:
        mol: RDKit molecule
    
    Returns:
        Dictionary of molecular features
    """
    try:
        features = {}
        
        # Basic properties
        features["molecular_weight"] = Descriptors.MolWt(mol)
        features["logp"] = Crippen.MolLogP(mol)
        features["hbd"] = Lipinski.NumHDonors(mol)
        features["hba"] = Lipinski.NumHAcceptors(mol)
        features["tpsa"] = Descriptors.TPSA(mol)
        features["rotatable_bonds"] = CalcNumRotatableBonds(mol)
        features["rings"] = Descriptors.RingCount(mol)
        features["aromatic_rings"] = Descriptors.NumAromaticRings(mol)
        
        # Advanced features
        features["peptide_like"] = _assess_peptide_like_character(mol)
        features["structural_uniqueness"] = _assess_structural_uniqueness(mol)
        features["specific_binding_groups"] = _count_specific_binding_groups(mol)
        features["steric_bulk"] = _assess_steric_bulk(mol)
        features["electronic_properties"] = _assess_electronic_properties(mol)
        
        # Additional derived features
        features["basic_nitrogen"] = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")))
        features["planar_rings"] = _count_planar_rings(mol)
        features["hydroxyl_groups"] = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OH]")))
        features["structural_complexity"] = _calculate_structural_complexity(mol)
        features["promiscuous_substructures"] = _count_promiscuous_substructures(mol)
        
        return features
        
    except Exception as e:
        logger.warning(f"Feature extraction failed: {e}")
        return {}


def _assess_peptide_like_character(mol: Chem.Mol) -> bool:
    """Assess if molecule has peptide-like character"""
    try:
        # Look for amide bonds
        amide_bonds = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)N")))
        
        # Look for amino acid-like features
        amino_groups = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[NX3;H2,H1]")))
        carboxyl_groups = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[OH]")))
        
        # Heuristic: multiple amide bonds or amino acid features
        return amide_bonds >= 2 or (amino_groups >= 1 and carboxyl_groups >= 1)
        
    except Exception:
        return False


def _assess_structural_uniqueness(mol: Chem.Mol) -> float:
    """Assess structural uniqueness (0-1)"""
    try:
        # Based on scaffold diversity and rare structural features
        
        # Get Murcko scaffold
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            scaffold_complexity = scaffold.GetNumHeavyAtoms() / mol.GetNumHeavyAtoms()
        except:
            scaffold_complexity = 0.5
        
        # Count heteroatoms (more heteroatoms = more unique)
        heteroatoms = CalcNumHeteroatoms(mol)
        hetero_ratio = heteroatoms / mol.GetNumHeavyAtoms()
        
        # Spiro and bridgehead atoms are rare
        spiro_atoms = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        bridgehead_atoms = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        
        uniqueness = (
            0.4 * scaffold_complexity +
            0.3 * min(1.0, hetero_ratio * 3) +
            0.3 * min(1.0, (spiro_atoms + bridgehead_atoms) * 0.1)
        )
        
        return np.clip(uniqueness, 0.0, 1.0)
        
    except Exception:
        return 0.5


def _count_specific_binding_groups(mol: Chem.Mol) -> int:
    """Count specific binding motifs"""
    try:
        count = 0
        
        # Hydrogen bond donors/acceptors in specific contexts
        specific_patterns = [
            "[OH]c",  # Phenolic OH
            "C(=O)[NH]",  # Amide
            "[NH]c1ccccc1",  # Aniline-like
            "c1ccc(O)cc1",  # Phenol
            "C(=O)O",  # Carboxylic acid
        ]
        
        for pattern in specific_patterns:
            try:
                matches = len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)))
                count += matches
            except:
                continue
        
        return count
        
    except Exception:
        return 0


def _assess_steric_bulk(mol: Chem.Mol) -> float:
    """Assess steric bulk (0-1)"""
    try:
        # Based on molecular volume and shape
        
        # Use number of rotatable bonds as proxy for flexibility
        rotbonds = CalcNumRotatableBonds(mol)
        
        # Heavy atom count as proxy for size
        heavy_atoms = mol.GetNumHeavyAtoms()
        
        # Ratio gives indication of bulk vs. flexibility
        if rotbonds == 0:
            bulk = min(1.0, heavy_atoms / 20.0)
        else:
            bulk = min(1.0, heavy_atoms / (rotbonds + 5.0) / 5.0)
        
        return bulk
        
    except Exception:
        return 0.5


def _assess_electronic_properties(mol: Chem.Mol) -> Dict[str, float]:
    """Assess electronic properties"""
    try:
        properties = {}
        
        # Electron-withdrawing groups
        ewg_patterns = ["[N+](=O)[O-]", "C(=O)", "C#N", "S(=O)(=O)"]
        ewg_count = 0
        for pattern in ewg_patterns:
            try:
                ewg_count += len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)))
            except:
                continue
        
        # Electron-donating groups
        edg_patterns = ["[OH]", "[NH2]", "[OCH3]"]
        edg_count = 0
        for pattern in edg_patterns:
            try:
                edg_count += len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)))
            except:
                continue
        
        properties["electron_withdrawing"] = min(1.0, ewg_count / 3.0)
        properties["electron_donating"] = min(1.0, edg_count / 3.0)
        
        return properties
        
    except Exception:
        return {"electron_withdrawing": 0.0, "electron_donating": 0.0}


def _count_planar_rings(mol: Chem.Mol) -> int:
    """Count planar ring systems (aromatic rings as proxy)"""
    try:
        return Descriptors.NumAromaticRings(mol)
    except:
        return 0


def _calculate_structural_complexity(mol: Chem.Mol) -> float:
    """Calculate structural complexity (0-1)"""
    try:
        # Bertz complexity index, normalized
        bertz = Descriptors.BertzCT(mol)
        
        # Normalize to 0-1 (typical range 0-100+)
        normalized = min(1.0, bertz / 100.0)
        
        return normalized
        
    except Exception:
        return 0.5


def _count_promiscuous_substructures(mol: Chem.Mol) -> int:
    """Count substructures known to bind promiscuously"""
    try:
        # Simple patterns associated with promiscuity
        promiscuous_patterns = [
            "c1ccc(O)c(O)c1",  # Catechol
            "C1CCCCC1",  # Cyclohexane (simple ring)
            "C(C)(C)C",  # Quaternary carbon
        ]
        
        count = 0
        for pattern in promiscuous_patterns:
            try:
                count += len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)))
            except:
                continue
        
        return count
        
    except Exception:
        return 0


def _assess_synthetic_accessibility(mol: Chem.Mol) -> float:
    """
    Assess synthetic accessibility (0-1, higher = more accessible)
    
    Args:
        mol: RDKit molecule
    
    Returns:
        Synthetic accessibility score
    """
    try:
        score = 0.8  # Start optimistic
        
        # Penalize high molecular weight
        mw = Descriptors.MolWt(mol)
        if mw > 500:
            score -= 0.3
        elif mw > 400:
            score -= 0.1
        
        # Penalize high flexibility
        rotbonds = CalcNumRotatableBonds(mol)
        if rotbonds > 10:
            score -= 0.3
        elif rotbonds > 7:
            score -= 0.1
        
        # Penalize complex ring systems
        rings = Descriptors.RingCount(mol)
        if rings > 4:
            score -= 0.2
        elif rings > 2:
            score -= 0.1
        
        # Penalize many stereocenters
        stereo_centers = len(Chem.FindMolChiralCenters(mol))
        if stereo_centers > 4:
            score -= 0.2
        elif stereo_centers > 2:
            score -= 0.1
        
        # Bonus for common fragments
        common_fragments = [
            "c1ccccc1",  # Benzene
            "C(C)(C)",   # Branched alkyl
            "C(=O)",     # Carbonyl
        ]
        
        for fragment in common_fragments:
            try:
                if len(mol.GetSubstructMatches(Chem.MolFromSmarts(fragment))) > 0:
                    score += 0.05
            except:
                continue
        
        return np.clip(score, 0.0, 1.0)
        
    except Exception:
        return 0.5


def _assess_structural_novelty(mol: Chem.Mol) -> float:
    """
    Assess structural novelty (0-1, higher = more novel)
    
    Args:
        mol: RDKit molecule
    
    Returns:
        Novelty score
    """
    try:
        score = 0.5  # Base novelty
        
        # Novel features increase score
        heteroatoms = CalcNumHeteroatoms(mol)
        hetero_ratio = heteroatoms / mol.GetNumHeavyAtoms()
        score += min(0.3, hetero_ratio * 1.5)
        
        # Rare structural features
        spiro_atoms = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        bridgehead_atoms = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        score += min(0.2, (spiro_atoms + bridgehead_atoms) * 0.05)
        
        # Scaffold complexity
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            scaffold_ratio = scaffold.GetNumHeavyAtoms() / mol.GetNumHeavyAtoms()
            if scaffold_ratio < 0.5:  # Complex substitution pattern
                score += 0.1
        except:
            pass
        
        # Uncommon ring sizes
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            ring_size = len(ring)
            if ring_size not in [5, 6]:  # Common ring sizes
                score += 0.05
        
        return np.clip(score, 0.0, 1.0)
        
    except Exception:
        return 0.5


# ============================================================================
# CROWE SCORING AGGREGATOR
# ============================================================================

def crowe_score(
    smiles: str,
    target_class: str = "GPCR",
    weights: Optional[Dict[str, float]] = None,
    config: Optional[Dict[str, Any]] = None
) -> Dict[str, float]:
    """
    Calculate Crowe Discovery Framework composite score
    
    Args:
        smiles: SMILES string
        target_class: Target protein class
        weights: Component weights (potency, selectivity, admet, synthesis, novelty)
        config: Configuration dictionary
    
    Returns:
        Dictionary with component scores and composite Crowe score
    """
    # Default weights
    if weights is None:
        weights = {
            "potency": 0.25,
            "selectivity": 0.20,
            "admet": 0.25,
            "synthesis": 0.20,
            "novelty": 0.10
        }
    
    # Validate weights sum to 1.0
    weight_sum = sum(weights.values())
    if abs(weight_sum - 1.0) > 0.01:
        raise ScoringError(f"Weights must sum to 1.0, got {weight_sum:.3f}")
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                "potency": 0.0,
                "selectivity": 0.0,
                "admet": 0.0,
                "synthesis": 0.0,
                "novelty": 0.0,
                "crowe_score": 0.0,
                "molecular_weight": 0.0,
                "drug_like": False
            }
        
        # Calculate component scores
        component_scores = _calculate_design_scores(
            smiles, 
            target_class, 
            list(weights.keys())
        )
        
        # Calculate composite Crowe score
        crowe_composite = sum(
            component_scores.get(component, 0.0) * weight
            for component, weight in weights.items()
        )
        
        # Additional molecular properties
        mw = Descriptors.MolWt(mol)
        drug_like = _passes_drug_like_filters(smiles, config)
        
        result = {
            **component_scores,
            "crowe_score": crowe_composite,
            "molecular_weight": mw,
            "drug_like": drug_like
        }
        
        return result
        
    except Exception as e:
        logger.error(f"Crowe scoring failed for {smiles}: {e}")
        return {
            "potency": 0.0,
            "selectivity": 0.0,
            "admet": 0.0,
            "synthesis": 0.0,
            "novelty": 0.0,
            "crowe_score": 0.0,
            "molecular_weight": 0.0,
            "drug_like": False
        }


def design_synthetic_drugs(
    target_class: str = "GPCR",
    n_compounds: int = 100,
    optimize_for: List[str] = None,
    seed_smiles: Optional[List[str]] = None,
    config: Optional[Dict[str, Any]] = None
) -> List[Dict[str, Any]]:
    """
    High-level API for synthetic drug design
    
    Args:
        target_class: Target protein class
        n_compounds: Number of compounds to generate
        optimize_for: Optimization objectives
        seed_smiles: Starting compounds (optional)
        config: Configuration dictionary
    
    Returns:
        List of designed compounds with scores
    """
    if optimize_for is None:
        optimize_for = ["potency", "selectivity", "admet"]
    
    # This is a placeholder for the actual generative algorithm
    # In practice, this would use molecular generation methods
    
    results = []
    
    # For demonstration, use some example SMILES
    example_smiles = [
        "CCO",  # Ethanol (simple)
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
        "CC1=CC2=C(C=C1C)N(C=N2)C3C(C(C(O3)CO)O)O",  # Riboflavin
    ]
    
    # Extend with seed compounds if provided
    if seed_smiles:
        example_smiles.extend(seed_smiles)
    
    # Score each compound
    for i, smiles in enumerate(example_smiles[:n_compounds]):
        try:
            scores = crowe_score(smiles, target_class, config=config)
            
            result = {
                "id": f"compound_{i+1}",
                "smiles": smiles,
                "target_class": target_class,
                **scores
            }
            
            results.append(result)
            
        except Exception as e:
            logger.warning(f"Failed to score {smiles}: {e}")
            continue
    
    # Sort by Crowe score (descending)
    results.sort(key=lambda x: x.get("crowe_score", 0.0), reverse=True)
    
    return results