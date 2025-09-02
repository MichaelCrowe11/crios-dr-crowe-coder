"""
CriOS Crowe Scoring Algorithm
Multi-component scoring system for compound prioritization
Implements the core Crowe Discovery Framework scoring methodology
"""

import logging
from typing import Dict, List, Optional, Union, Any, Tuple
from dataclasses import dataclass
import math

import numpy as np
import pandas as pd

from ..core.molecule import Molecule
from ..core.compound import Compound, CompoundLibrary
from ..core.descriptors import MolecularDescriptors, DescriptorCalculator
from ..config.schemas import SafetyCategory
from ..exceptions import CalculationError

logger = logging.getLogger(__name__)


@dataclass
class CroweScoreComponents:
    """Container for individual Crowe score components"""
    
    # Core components (sum to 1.0)
    novelty: float = 0.0
    drug_likeness: float = 0.0
    synthetic_accessibility: float = 0.0
    safety_profile: float = 0.0
    ethical_alignment: float = 0.0
    
    # Additional metrics
    efficacy_potential: Optional[float] = None
    selectivity_profile: Optional[float] = None
    natural_origin_bonus: Optional[float] = None
    
    # Meta information
    calculation_confidence: float = 1.0
    missing_data_penalty: float = 0.0
    
    def total_score(self, weights: Optional[Dict[str, float]] = None) -> float:
        """Calculate weighted total score"""
        default_weights = {
            "novelty": 0.25,
            "drug_likeness": 0.20, 
            "synthetic_accessibility": 0.15,
            "safety_profile": 0.20,
            "ethical_alignment": 0.20
        }
        
        if weights:
            default_weights.update(weights)
        
        score = (
            self.novelty * default_weights["novelty"] +
            self.drug_likeness * default_weights["drug_likeness"] +
            self.synthetic_accessibility * default_weights["synthetic_accessibility"] +
            self.safety_profile * default_weights["safety_profile"] +
            self.ethical_alignment * default_weights["ethical_alignment"]
        )
        
        # Apply confidence and penalty adjustments
        score *= self.calculation_confidence
        score -= self.missing_data_penalty
        
        return max(0.0, min(1.0, score))


class CroweScorer:
    """
    Crowe Discovery Framework Scoring Engine
    
    Implements the multi-component scoring system for compound prioritization:
    - Novelty (structural uniqueness)
    - Drug-likeness (ADMET and physicochemical properties)
    - Synthetic accessibility (ease of synthesis)
    - Safety profile (toxicity and side effects)
    - Ethical alignment (mission compatibility)
    """
    
    def __init__(
        self,
        reference_library: Optional[CompoundLibrary] = None,
        weights: Optional[Dict[str, float]] = None,
        include_bonuses: bool = True
    ):
        """
        Initialize Crowe scorer
        
        Args:
            reference_library: Reference compound library for novelty assessment
            weights: Custom component weights (must sum to 1.0)
            include_bonuses: Whether to include bonus components
        """
        self.reference_library = reference_library
        self.include_bonuses = include_bonuses
        
        # Set default weights
        self.default_weights = {
            "novelty": 0.25,
            "drug_likeness": 0.20,
            "synthetic_accessibility": 0.15, 
            "safety_profile": 0.20,
            "ethical_alignment": 0.20
        }
        
        if weights:
            self._validate_weights(weights)
            self.default_weights.update(weights)
        
        # Initialize sub-scorers
        self.descriptor_calculator = DescriptorCalculator()
        
        # Caching for performance
        self._novelty_cache: Dict[str, float] = {}
        self._descriptor_cache: Dict[str, MolecularDescriptors] = {}
    
    def _validate_weights(self, weights: Dict[str, float]) -> None:
        """Validate that weights sum to approximately 1.0"""
        total = sum(weights.values())
        if abs(total - 1.0) > 0.01:
            raise ValueError(f"Weights must sum to 1.0, got {total}")
    
    def score_single(
        self,
        compound: Union[Molecule, Compound],
        return_components: bool = False
    ) -> Union[float, Tuple[float, CroweScoreComponents]]:
        """
        Calculate Crowe score for a single compound
        
        Args:
            compound: Compound to score
            return_components: Whether to return component breakdown
            
        Returns:
            Crowe score (0-1) or tuple of (score, components)
        """
        try:
            # Convert Molecule to Compound if needed
            if isinstance(compound, Molecule) and not isinstance(compound, Compound):
                compound = Compound(compound)
            
            # Calculate components
            components = CroweScoreComponents()
            
            # Component 1: Novelty Score
            components.novelty = self._calculate_novelty_score(compound)
            
            # Component 2: Drug-likeness Score
            components.drug_likeness = self._calculate_drug_likeness_score(compound)
            
            # Component 3: Synthetic Accessibility Score
            components.synthetic_accessibility = self._calculate_synthetic_accessibility_score(compound)
            
            # Component 4: Safety Profile Score
            components.safety_profile = self._calculate_safety_score(compound)
            
            # Component 5: Ethical Alignment Score
            components.ethical_alignment = self._calculate_ethical_score(compound)
            
            # Optional bonus components
            if self.include_bonuses:
                components.efficacy_potential = self._calculate_efficacy_potential(compound)
                components.natural_origin_bonus = self._calculate_natural_origin_bonus(compound)
            
            # Meta-scoring adjustments
            components.calculation_confidence = self._calculate_confidence(compound)
            components.missing_data_penalty = self._calculate_missing_data_penalty(compound)
            
            # Calculate final score
            final_score = components.total_score(self.default_weights)
            
            if return_components:
                return final_score, components
            else:
                return final_score
                
        except Exception as e:
            logger.error(f"Failed to calculate Crowe score for {compound.mol_id}: {e}")
            return 0.0 if not return_components else (0.0, CroweScoreComponents())
    
    def _calculate_novelty_score(self, compound: Union[Molecule, Compound]) -> float:
        """Calculate structural novelty score (0-1, higher = more novel)"""
        try:
            # Check cache
            if compound.mol_id in self._novelty_cache:
                return self._novelty_cache[compound.mol_id]
            
            if not self.reference_library:
                # Without reference library, use heuristic based on structural complexity
                return self._heuristic_novelty_score(compound)
            
            # Calculate similarity to reference compounds
            from ..core.similarity import SimilaritySearch
            
            similarity_search = SimilaritySearch()
            similarities = []
            
            for ref_compound in self.reference_library:
                if ref_compound.mol_id != compound.mol_id:
                    sim = similarity_search.calculate_similarity(compound, ref_compound)
                    if sim is not None:
                        similarities.append(sim)
            
            if not similarities:
                novelty_score = 0.7  # Default moderate novelty
            else:
                # Novelty = 1 - max_similarity
                max_similarity = max(similarities)
                novelty_score = 1.0 - max_similarity
                
                # Apply sigmoid transformation for better distribution
                novelty_score = 1.0 / (1.0 + math.exp(-5 * (novelty_score - 0.5)))
            
            # Cache result
            self._novelty_cache[compound.mol_id] = novelty_score
            
            return novelty_score
            
        except Exception as e:
            logger.warning(f"Failed to calculate novelty score: {e}")
            return 0.5  # Default moderate novelty
    
    def _heuristic_novelty_score(self, compound: Union[Molecule, Compound]) -> float:
        """Heuristic novelty based on structural features"""
        try:
            descriptors = self._get_descriptors(compound)
            if not descriptors:
                return 0.5
            
            novelty = 0.5  # Base score
            
            # Rare structural features increase novelty
            if descriptors.num_spiro_atoms and descriptors.num_spiro_atoms > 0:
                novelty += 0.1
            
            if descriptors.num_bridgehead_atoms and descriptors.num_bridgehead_atoms > 2:
                novelty += 0.1
            
            # Unique ring systems
            if descriptors.saturated_rings and descriptors.aromatic_rings:
                if descriptors.saturated_rings > descriptors.aromatic_rings:
                    novelty += 0.05
            
            # High scaffold diversity
            if 'scaffold_diversity' in descriptors.custom_descriptors:
                scaffold_div = descriptors.custom_descriptors['scaffold_diversity']
                novelty += 0.2 * scaffold_div
            
            # Uncommon molecular weight range
            if descriptors.molecular_weight:
                if descriptors.molecular_weight > 600 or descriptors.molecular_weight < 200:
                    novelty += 0.05
            
            return min(1.0, novelty)
            
        except Exception:
            return 0.5
    
    def _calculate_drug_likeness_score(self, compound: Union[Molecule, Compound]) -> float:
        """Calculate drug-likeness score (0-1, higher = more drug-like)"""
        try:
            descriptors = self._get_descriptors(compound)
            if not descriptors:
                return 0.3
            
            score = 1.0  # Start with perfect score
            
            # Lipinski Rule of Five violations
            lipinski_violations = 0
            
            if descriptors.molecular_weight and descriptors.molecular_weight > 500:
                lipinski_violations += 1
            if descriptors.logp and descriptors.logp > 5:
                lipinski_violations += 1
            if descriptors.hbd and descriptors.hbd > 5:
                lipinski_violations += 1
            if descriptors.hba and descriptors.hba > 10:
                lipinski_violations += 1
            
            # Penalty for violations
            score -= lipinski_violations * 0.15
            
            # Veber Rule violations
            veber_violations = 0
            if descriptors.rotatable_bonds and descriptors.rotatable_bonds > 10:
                veber_violations += 1
            if descriptors.tpsa and descriptors.tpsa > 140:
                veber_violations += 1
            
            score -= veber_violations * 0.1
            
            # Additional drug-likeness factors
            if descriptors.aromatic_rings and descriptors.aromatic_rings > 4:
                score -= 0.1
            
            if descriptors.aliphatic_rings and descriptors.aliphatic_rings > 2:
                score -= 0.05
            
            # Bonus for ideal ranges
            if descriptors.logp and 1.5 <= descriptors.logp <= 3.5:
                score += 0.05
            
            if descriptors.tpsa and 40 <= descriptors.tpsa <= 90:
                score += 0.05
            
            return max(0.0, min(1.0, score))
            
        except Exception as e:
            logger.warning(f"Failed to calculate drug-likeness score: {e}")
            return 0.5
    
    def _calculate_synthetic_accessibility_score(self, compound: Union[Molecule, Compound]) -> float:
        """Calculate synthetic accessibility score (0-1, higher = more accessible)"""
        try:
            # Check if compound has synthesis info
            if isinstance(compound, Compound) and compound.synthesis.synthetic_accessibility:
                # SA score is typically 1-10, convert to 0-1 (inverted)
                sa_score = compound.synthesis.synthetic_accessibility
                return max(0.0, (10 - sa_score) / 9)
            
            # Use heuristic based on molecular complexity
            descriptors = self._get_descriptors(compound)
            if not descriptors:
                return 0.5
            
            score = 0.8  # Start optimistic
            
            # Penalize complexity factors
            if descriptors.molecular_weight and descriptors.molecular_weight > 500:
                score -= 0.2
            elif descriptors.molecular_weight and descriptors.molecular_weight > 400:
                score -= 0.1
            
            if descriptors.rotatable_bonds and descriptors.rotatable_bonds > 10:
                score -= 0.3
            elif descriptors.rotatable_bonds and descriptors.rotatable_bonds > 7:
                score -= 0.1
            
            # Ring complexity penalties
            total_rings = (descriptors.aromatic_rings or 0) + (descriptors.aliphatic_rings or 0)
            if total_rings > 4:
                score -= 0.2
            elif total_rings > 2:
                score -= 0.1
            
            # Stereochemistry penalty (if available)
            if hasattr(compound, 'rdkit_mol') and compound.rdkit_mol:
                from rdkit import Chem
                stereo_centers = len(Chem.FindMolChiralCenters(compound.rdkit_mol))
                if stereo_centers > 4:
                    score -= 0.2
                elif stereo_centers > 2:
                    score -= 0.1
            
            # Heteroatom bonus (makes synthesis easier)
            if descriptors.num_heteroatoms and descriptors.num_heteroatoms > 2:
                score += 0.05
            
            return max(0.0, min(1.0, score))
            
        except Exception as e:
            logger.warning(f"Failed to calculate synthetic accessibility score: {e}")
            return 0.5
    
    def _calculate_safety_score(self, compound: Union[Molecule, Compound]) -> float:
        """Calculate safety profile score (0-1, higher = safer)"""
        try:
            score = 1.0  # Start with perfect safety
            
            # Check safety category if compound
            if isinstance(compound, Compound):
                if compound.safety_category == SafetyCategory.YELLOW:
                    score -= 0.3
                elif compound.safety_category == SafetyCategory.RED:
                    score -= 0.7
                
                # Check ADMET properties
                admet = compound.admet
                if admet.hepatotoxicity:
                    score -= 0.4
                if admet.cardiotoxicity:
                    score -= 0.3
                if admet.mutagenicity:
                    score -= 0.2
                if admet.carcinogenicity:
                    score -= 0.5
                
                # Positive ADMET factors
                if admet.bbb_permeability and admet.bbb_permeability < 0.1:
                    # Low BBB permeability can be good for non-CNS drugs
                    pass
                
                if admet.oral_bioavailability and admet.oral_bioavailability > 0.3:
                    score += 0.05
            
            # Structural safety heuristics
            descriptors = self._get_descriptors(compound)
            if descriptors:
                # Reactive functional groups (simplified detection)
                # This is a placeholder - real implementation would use SMARTS patterns
                
                # Very high LogP can indicate toxicity
                if descriptors.logp and descriptors.logp > 6:
                    score -= 0.2
                
                # Very low or high molecular weight concerns
                if descriptors.molecular_weight:
                    if descriptors.molecular_weight < 150 or descriptors.molecular_weight > 800:
                        score -= 0.1
            
            return max(0.0, min(1.0, score))
            
        except Exception as e:
            logger.warning(f"Failed to calculate safety score: {e}")
            return 0.5
    
    def _calculate_ethical_score(self, compound: Union[Molecule, Compound]) -> float:
        """Calculate ethical alignment score (0-1, higher = better aligned)"""
        try:
            score = 0.7  # Base ethical alignment
            
            if isinstance(compound, Compound):
                # Natural origin bonus
                from ..core.compound import CompoundOrigin
                if compound.origin in [CompoundOrigin.NATURAL, CompoundOrigin.HYBRID]:
                    score += 0.1
                
                # Therapeutic area alignment
                ethical_areas = {
                    "neurodegeneration", "alzheimer", "parkinson", "environmental_health",
                    "regenerative_medicine", "aging", "sustainability"
                }
                
                compound_areas = {area.lower() for area in compound.metadata.therapeutic_areas}
                if compound_areas.intersection(ethical_areas):
                    score += 0.15
                
                # Patent and IP considerations
                if compound.patent_info.freedom_to_operate is True:
                    score += 0.05
                elif compound.patent_info.freedom_to_operate is False:
                    score -= 0.1
                
                # Collaboration and open science
                if compound.metadata.publications:
                    score += 0.05
                
                # Quality flags check
                if "unethical" in compound.metadata.quality_flags:
                    score -= 0.3
                if "sustainable" in compound.metadata.quality_flags:
                    score += 0.1
            
            # Molecular characteristics aligned with ethical goals
            descriptors = self._get_descriptors(compound)
            if descriptors and 'natural_product_score' in descriptors.custom_descriptors:
                np_score = descriptors.custom_descriptors['natural_product_score']
                score += 0.1 * np_score  # Bonus for natural product likeness
            
            return max(0.0, min(1.0, score))
            
        except Exception as e:
            logger.warning(f"Failed to calculate ethical score: {e}")
            return 0.7
    
    def _calculate_efficacy_potential(self, compound: Union[Molecule, Compound]) -> float:
        """Calculate potential efficacy based on biological activity data"""
        try:
            if not isinstance(compound, Compound):
                return 0.5
            
            if not compound.biological_activities:
                return 0.5  # No data available
            
            # Find best (lowest) IC50 values
            best_activities = []
            for activity in compound.biological_activities:
                if activity.activity_type.upper() in ["IC50", "KI", "EC50"]:
                    # Convert to nM if needed
                    value_nm = activity.activity_value
                    if activity.activity_unit.upper() == "UM":
                        value_nm *= 1000
                    elif activity.activity_unit.upper() == "MM":
                        value_nm *= 1000000
                    
                    best_activities.append(value_nm)
            
            if not best_activities:
                return 0.5
            
            best_activity = min(best_activities)
            
            # Convert to score (0-1)
            if best_activity <= 1:  # nM
                return 1.0
            elif best_activity <= 10:
                return 0.9
            elif best_activity <= 100:
                return 0.7
            elif best_activity <= 1000:
                return 0.5
            elif best_activity <= 10000:
                return 0.3
            else:
                return 0.1
                
        except Exception as e:
            logger.warning(f"Failed to calculate efficacy potential: {e}")
            return 0.5
    
    def _calculate_natural_origin_bonus(self, compound: Union[Molecule, Compound]) -> float:
        """Calculate bonus for natural product origin"""
        try:
            if not isinstance(compound, Compound):
                return 0.0
            
            from ..core.compound import CompoundOrigin
            
            if compound.origin == CompoundOrigin.NATURAL:
                return 0.1
            elif compound.origin == CompoundOrigin.HYBRID:
                return 0.05
            else:
                # Check natural product likeness score
                descriptors = self._get_descriptors(compound)
                if descriptors and 'natural_product_score' in descriptors.custom_descriptors:
                    np_score = descriptors.custom_descriptors['natural_product_score']
                    return 0.05 * np_score
                return 0.0
                
        except Exception:
            return 0.0
    
    def _calculate_confidence(self, compound: Union[Molecule, Compound]) -> float:
        """Calculate confidence in the scoring based on data availability"""
        try:
            confidence = 1.0
            
            # Penalize if molecule is invalid
            if not compound.is_valid():
                confidence -= 0.5
            
            # Check data availability
            if isinstance(compound, Compound):
                data_points = 0
                total_possible = 5
                
                # Biological activity data
                if compound.biological_activities:
                    data_points += 1
                
                # ADMET data
                admet_fields = [
                    compound.admet.oral_bioavailability,
                    compound.admet.bbb_permeability,
                    compound.admet.hepatotoxicity,
                    compound.admet.cardiotoxicity
                ]
                if any(field is not None for field in admet_fields):
                    data_points += 1
                
                # Synthesis data
                if compound.synthesis.synthetic_accessibility:
                    data_points += 1
                
                # Patent data
                if compound.patent_info.patent_numbers or compound.patent_info.freedom_to_operate is not None:
                    data_points += 1
                
                # Therapeutic area info
                if compound.metadata.therapeutic_areas:
                    data_points += 1
                
                # Adjust confidence based on data availability
                data_completeness = data_points / total_possible
                confidence *= (0.6 + 0.4 * data_completeness)  # Scale from 0.6 to 1.0
            
            return max(0.1, min(1.0, confidence))
            
        except Exception:
            return 0.8
    
    def _calculate_missing_data_penalty(self, compound: Union[Molecule, Compound]) -> float:
        """Calculate penalty for missing critical data"""
        try:
            penalty = 0.0
            
            # Critical missing data penalties
            if not compound.is_valid():
                penalty += 0.3
            
            descriptors = self._get_descriptors(compound)
            if not descriptors:
                penalty += 0.2
            
            if isinstance(compound, Compound):
                # No biological activity data for compounds
                if not compound.biological_activities:
                    penalty += 0.1
                
                # Unknown safety category
                if compound.safety_category == SafetyCategory.GREEN:
                    # This might mean "unknown" rather than "safe"
                    if not any([
                        compound.admet.hepatotoxicity is not None,
                        compound.admet.cardiotoxicity is not None,
                        compound.admet.mutagenicity is not None
                    ]):
                        penalty += 0.05
            
            return min(0.5, penalty)  # Cap penalty at 0.5
            
        except Exception:
            return 0.0
    
    def _get_descriptors(self, compound: Union[Molecule, Compound]) -> Optional[MolecularDescriptors]:
        """Get or calculate molecular descriptors with caching"""
        cache_key = compound.mol_id
        
        if cache_key in self._descriptor_cache:
            return self._descriptor_cache[cache_key]
        
        descriptors = self.descriptor_calculator.calculate_single(compound)
        if descriptors:
            self._descriptor_cache[cache_key] = descriptors
        
        return descriptors
    
    def score_library(
        self,
        compound_library: CompoundLibrary,
        return_components: bool = False,
        parallel: bool = True
    ) -> Union[Dict[str, float], Dict[str, Tuple[float, CroweScoreComponents]]]:
        """
        Score all compounds in a library
        
        Args:
            compound_library: Library to score
            return_components: Whether to return component breakdowns
            parallel: Whether to use parallel processing
            
        Returns:
            Dictionary mapping compound IDs to scores (and components if requested)
        """
        results = {}
        
        if parallel:
            return self._score_library_parallel(compound_library, return_components)
        else:
            return self._score_library_sequential(compound_library, return_components)
    
    def _score_library_sequential(
        self,
        compound_library: CompoundLibrary,
        return_components: bool
    ) -> Union[Dict[str, float], Dict[str, Tuple[float, CroweScoreComponents]]]:
        """Sequential library scoring"""
        results = {}
        
        for compound in compound_library:
            try:
                result = self.score_single(compound, return_components)
                results[compound.mol_id] = result
            except Exception as e:
                logger.error(f"Failed to score compound {compound.mol_id}: {e}")
                if return_components:
                    results[compound.mol_id] = (0.0, CroweScoreComponents())
                else:
                    results[compound.mol_id] = 0.0
        
        return results
    
    def _score_library_parallel(
        self,
        compound_library: CompoundLibrary,
        return_components: bool
    ) -> Union[Dict[str, float], Dict[str, Tuple[float, CroweScoreComponents]]]:
        """Parallel library scoring"""
        try:
            from joblib import Parallel, delayed
            import multiprocessing as mp
            
            compounds = list(compound_library)
            n_jobs = min(mp.cpu_count(), 8)
            
            # Create partial function for scoring
            def score_compound(comp):
                try:
                    return comp.mol_id, self.score_single(comp, return_components)
                except Exception as e:
                    logger.error(f"Failed to score compound {comp.mol_id}: {e}")
                    if return_components:
                        return comp.mol_id, (0.0, CroweScoreComponents())
                    else:
                        return comp.mol_id, 0.0
            
            # Process in parallel
            results_list = Parallel(n_jobs=n_jobs)(
                delayed(score_compound)(compound) for compound in compounds
            )
            
            # Convert to dictionary
            results = dict(results_list)
            
            return results
            
        except ImportError:
            logger.warning("joblib not available, falling back to sequential scoring")
            return self._score_library_sequential(compound_library, return_components)
    
    def rank_compounds(
        self,
        compound_library: CompoundLibrary,
        limit: Optional[int] = None
    ) -> List[Tuple[str, float, CroweScoreComponents]]:
        """
        Rank compounds by Crowe score
        
        Args:
            compound_library: Library to rank
            limit: Maximum number of results
            
        Returns:
            List of (compound_id, score, components) tuples, sorted by score
        """
        # Score all compounds with components
        scores = self.score_library(compound_library, return_components=True)
        
        # Sort by score (descending)
        ranked = sorted(
            [(comp_id, score, components) for comp_id, (score, components) in scores.items()],
            key=lambda x: x[1],
            reverse=True
        )
        
        if limit:
            ranked = ranked[:limit]
        
        return ranked
    
    def clear_cache(self) -> None:
        """Clear all caches"""
        self._novelty_cache.clear()
        self._descriptor_cache.clear()
    
    def get_cache_stats(self) -> Dict[str, int]:
        """Get cache statistics"""
        return {
            "novelty_cache_size": len(self._novelty_cache),
            "descriptor_cache_size": len(self._descriptor_cache)
        }