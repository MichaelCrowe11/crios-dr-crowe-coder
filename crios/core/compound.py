"""
CriOS Compound Classes
Extended molecular representation with biological and discovery metadata
"""

import logging
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional, Union, Any, Set
from pathlib import Path
from enum import Enum
import json

import numpy as np
import pandas as pd
from pydantic import BaseModel, Field, validator

from .molecule import Molecule, MoleculeCollection, MolecularProperties
from ..config.schemas import SafetyCategory
from ..exceptions import CompoundError, ValidationError

logger = logging.getLogger(__name__)


class CompoundOrigin(str, Enum):
    """Origin classification for compounds"""
    NATURAL = "natural"
    SYNTHETIC = "synthetic" 
    SEMI_SYNTHETIC = "semi_synthetic"
    HYBRID = "hybrid"
    UNKNOWN = "unknown"


class DiscoveryStage(str, Enum):
    """Discovery pipeline stage"""
    VIRTUAL = "virtual"
    SCREENING = "screening"
    HIT = "hit"
    LEAD = "lead"
    CANDIDATE = "candidate"
    CLINICAL = "clinical"
    APPROVED = "approved"


class BiologicalActivity(BaseModel):
    """Biological activity data container"""
    target_name: str
    activity_type: str  # IC50, Ki, EC50, etc.
    activity_value: float
    activity_unit: str = "nM"
    assay_description: Optional[str] = None
    reference: Optional[str] = None
    confidence: Optional[float] = Field(None, ge=0.0, le=1.0)


class ADMETProperties(BaseModel):
    """ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) properties"""
    
    # Absorption
    oral_bioavailability: Optional[float] = Field(None, ge=0.0, le=1.0)
    caco2_permeability: Optional[float] = None
    pgp_substrate: Optional[bool] = None
    
    # Distribution
    bbb_permeability: Optional[float] = Field(None, ge=0.0, le=1.0)
    plasma_protein_binding: Optional[float] = Field(None, ge=0.0, le=1.0)
    volume_distribution: Optional[float] = None
    
    # Metabolism
    cyp1a2_inhibitor: Optional[bool] = None
    cyp2c9_inhibitor: Optional[bool] = None
    cyp2c19_inhibitor: Optional[bool] = None
    cyp2d6_inhibitor: Optional[bool] = None
    cyp3a4_inhibitor: Optional[bool] = None
    
    # Excretion
    clearance: Optional[float] = None
    half_life: Optional[float] = None
    renal_clearance: Optional[float] = None
    
    # Toxicity
    hepatotoxicity: Optional[bool] = None
    cardiotoxicity: Optional[bool] = None
    mutagenicity: Optional[bool] = None
    carcinogenicity: Optional[bool] = None
    
    # Custom ADMET properties
    custom_admet: Dict[str, Any] = Field(default_factory=dict)


class SynthesisInfo(BaseModel):
    """Synthesis and accessibility information"""
    synthetic_accessibility: Optional[float] = Field(None, ge=0.0, le=10.0)
    retrosynthetic_complexity: Optional[float] = None
    starting_materials: List[str] = Field(default_factory=list)
    synthesis_steps: Optional[int] = None
    synthesis_route: Optional[str] = None
    estimated_cost: Optional[float] = None
    scalability_score: Optional[float] = Field(None, ge=0.0, le=1.0)


class PatentInfo(BaseModel):
    """Patent and IP information"""
    patent_numbers: List[str] = Field(default_factory=list)
    patent_status: Optional[str] = None
    freedom_to_operate: Optional[bool] = None
    exclusivity_expiry: Optional[datetime] = None
    inventor: Optional[str] = None
    assignee: Optional[str] = None


class CompoundScores(BaseModel):
    """Comprehensive scoring for compound prioritization"""
    
    # Crowe Scoring Methodology
    crowe_score: Optional[float] = Field(None, ge=0.0, le=1.0)
    
    # Individual component scores
    novelty_score: Optional[float] = Field(None, ge=0.0, le=1.0)
    drug_likeness_score: Optional[float] = Field(None, ge=0.0, le=1.0)
    synthetic_accessibility_score: Optional[float] = Field(None, ge=0.0, le=1.0)
    safety_score: Optional[float] = Field(None, ge=0.0, le=1.0)
    ethical_score: Optional[float] = Field(None, ge=0.0, le=1.0)
    
    # Multi-objective scores
    efficacy_score: Optional[float] = Field(None, ge=0.0, le=1.0)
    selectivity_score: Optional[float] = Field(None, ge=0.0, le=1.0)
    
    # Custom scores
    custom_scores: Dict[str, float] = Field(default_factory=dict)


@dataclass
class CompoundMetadata:
    """Extended metadata for compound discovery"""
    
    # Discovery information
    discovery_date: Optional[datetime] = None
    discoverer: Optional[str] = None
    discovery_method: Optional[str] = None
    discovery_stage: DiscoveryStage = DiscoveryStage.VIRTUAL
    
    # Classification
    therapeutic_areas: List[str] = field(default_factory=list)
    mechanism_of_action: Optional[str] = None
    compound_class: Optional[str] = None
    
    # Research status
    publications: List[str] = field(default_factory=list)
    clinical_trials: List[str] = field(default_factory=list)
    
    # Quality flags
    quality_flags: Set[str] = field(default_factory=set)
    validation_status: Optional[str] = None
    
    # Collaboration
    collaborators: List[str] = field(default_factory=list)
    institutions: List[str] = field(default_factory=list)


class Compound(Molecule):
    """
    Enhanced molecular representation for drug discovery
    Extends base Molecule class with biological, synthesis, and discovery metadata
    """
    
    def __init__(
        self,
        identifier: Union[str, Molecule],
        mol_id: Optional[str] = None,
        name: Optional[str] = None,
        source: Optional[str] = None,
        origin: CompoundOrigin = CompoundOrigin.UNKNOWN,
        **kwargs
    ):
        """
        Initialize compound from molecular identifier
        
        Args:
            identifier: SMILES, InChI, or Molecule object
            mol_id: Unique compound identifier
            name: Compound name
            source: Source database or reference
            origin: Compound origin classification
        """
        # Initialize base Molecule
        if isinstance(identifier, Molecule):
            super().__init__(identifier.rdkit_mol, mol_id, name, source)
        else:
            super().__init__(identifier, mol_id, name, source)
        
        # Compound-specific properties
        self.origin = origin
        self.safety_category = SafetyCategory.GREEN
        
        # Biological data
        self.biological_activities: List[BiologicalActivity] = []
        self.admet = ADMETProperties()
        
        # Synthesis and IP
        self.synthesis = SynthesisInfo()
        self.patent_info = PatentInfo()
        
        # Scoring
        self.scores = CompoundScores()
        
        # Discovery metadata
        self.metadata = CompoundMetadata()
        
        # Related compounds
        self.analogs: Set[str] = set()  # IDs of analog compounds
        self.parent_compound: Optional[str] = None
        self.derivative_compounds: Set[str] = set()
    
    def add_biological_activity(
        self,
        target_name: str,
        activity_type: str,
        activity_value: float,
        activity_unit: str = "nM",
        **kwargs
    ) -> None:
        """Add biological activity data"""
        activity = BiologicalActivity(
            target_name=target_name,
            activity_type=activity_type,
            activity_value=activity_value,
            activity_unit=activity_unit,
            **kwargs
        )
        self.biological_activities.append(activity)
    
    def get_activities_for_target(self, target_name: str) -> List[BiologicalActivity]:
        """Get all biological activities for a specific target"""
        return [
            activity for activity in self.biological_activities
            if activity.target_name.lower() == target_name.lower()
        ]
    
    def get_best_activity(self, target_name: str, activity_type: str = "IC50") -> Optional[BiologicalActivity]:
        """Get best (lowest) activity value for target and activity type"""
        activities = [
            activity for activity in self.biological_activities
            if activity.target_name.lower() == target_name.lower() and
               activity.activity_type.upper() == activity_type.upper()
        ]
        
        if not activities:
            return None
        
        return min(activities, key=lambda x: x.activity_value)
    
    def set_safety_category(self, category: SafetyCategory, reason: Optional[str] = None) -> None:
        """Set safety category with optional reason"""
        self.safety_category = category
        if reason:
            self.metadata.quality_flags.add(f"safety_{category.value}_{reason}")
    
    def add_quality_flag(self, flag: str) -> None:
        """Add quality control flag"""
        self.metadata.quality_flags.add(flag)
    
    def is_hit(self, threshold: float = 10.0, activity_type: str = "IC50") -> bool:
        """Check if compound is a hit based on activity threshold"""
        for activity in self.biological_activities:
            if activity.activity_type.upper() == activity_type.upper():
                if activity.activity_unit == "uM":
                    threshold_nm = threshold * 1000
                elif activity.activity_unit == "nM":
                    threshold_nm = threshold
                else:
                    continue  # Skip unknown units
                
                if activity.activity_value <= threshold_nm:
                    return True
        return False
    
    def is_lead(self, potency_threshold: float = 1.0, selectivity: float = 10.0) -> bool:
        """Check if compound qualifies as a lead based on potency and selectivity"""
        # Check potency
        has_potent_activity = any(
            activity.activity_value <= potency_threshold * (1000 if activity.activity_unit == "uM" else 1)
            for activity in self.biological_activities
            if activity.activity_type.upper() == "IC50"
        )
        
        if not has_potent_activity:
            return False
        
        # Additional lead criteria can be added here
        # (selectivity, ADMET, drug-likeness, etc.)
        
        return self.is_drug_like()
    
    def calculate_crowe_score(
        self,
        weights: Optional[Dict[str, float]] = None
    ) -> float:
        """
        Calculate Crowe Score using multi-objective optimization
        
        Args:
            weights: Custom weights for score components
        """
        default_weights = {
            "novelty": 0.25,
            "drug_likeness": 0.20,
            "synthetic_accessibility": 0.15,
            "safety": 0.20,
            "ethical_alignment": 0.20
        }
        
        if weights:
            default_weights.update(weights)
        
        score_components = {}
        
        # Novelty score (based on similarity to known compounds)
        novelty = self.scores.novelty_score or self._calculate_novelty_score()
        score_components["novelty"] = novelty
        
        # Drug-likeness score
        drug_likeness = self.scores.drug_likeness_score or self._calculate_drug_likeness_score()
        score_components["drug_likeness"] = drug_likeness
        
        # Synthetic accessibility score
        synth_acc = self.scores.synthetic_accessibility_score or self._calculate_synthetic_accessibility_score()
        score_components["synthetic_accessibility"] = synth_acc
        
        # Safety score
        safety = self.scores.safety_score or self._calculate_safety_score()
        score_components["safety"] = safety
        
        # Ethical alignment score
        ethical = self.scores.ethical_score or self._calculate_ethical_score()
        score_components["ethical_alignment"] = ethical
        
        # Calculate weighted sum
        crowe_score = sum(
            score_components[component] * default_weights[component]
            for component in score_components
        )
        
        self.scores.crowe_score = crowe_score
        return crowe_score
    
    def _calculate_novelty_score(self) -> float:
        """Calculate novelty score based on structural uniqueness"""
        # Placeholder implementation
        # In practice, this would compare against known compound databases
        return 0.7  # Assume moderate novelty
    
    def _calculate_drug_likeness_score(self) -> float:
        """Calculate drug-likeness score"""
        score = 1.0
        
        # Lipinski violations
        if not self.passes_lipinski():
            score -= 0.3
        
        # Veber violations
        if not self.passes_veber():
            score -= 0.2
        
        # Additional drug-like criteria
        if self.properties.aromatic_rings and self.properties.aromatic_rings > 4:
            score -= 0.1
        
        return max(0.0, score)
    
    def _calculate_synthetic_accessibility_score(self) -> float:
        """Calculate synthetic accessibility score"""
        if self.synthesis.synthetic_accessibility:
            # SA score is typically 0-10, normalize to 0-1 (inverted)
            return max(0.0, (10 - self.synthesis.synthetic_accessibility) / 10)
        return 0.5  # Default moderate accessibility
    
    def _calculate_safety_score(self) -> float:
        """Calculate safety score based on ADMET and toxicity predictions"""
        score = 1.0
        
        # ADMET penalties
        if self.admet.hepatotoxicity:
            score -= 0.4
        if self.admet.cardiotoxicity:
            score -= 0.3
        if self.admet.mutagenicity:
            score -= 0.2
        if self.admet.carcinogenicity:
            score -= 0.5
        
        # Safety category penalties
        if self.safety_category == SafetyCategory.YELLOW:
            score -= 0.2
        elif self.safety_category == SafetyCategory.RED:
            score -= 0.8
        
        return max(0.0, score)
    
    def _calculate_ethical_score(self) -> float:
        """Calculate ethical alignment score"""
        score = 1.0
        
        # Natural product bonus (ethical sourcing assumed)
        if self.origin == CompoundOrigin.NATURAL:
            score += 0.1
        
        # Therapeutic area alignment
        ethical_areas = {
            "neurodegeneration", "alzheimer", "parkinson", 
            "environmental_health", "regenerative_medicine"
        }
        
        if any(area.lower() in ethical_areas for area in self.metadata.therapeutic_areas):
            score += 0.1
        
        return min(1.0, score)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert compound to comprehensive dictionary representation"""
        base_dict = super().to_dict()
        
        compound_dict = {
            **base_dict,
            
            # Compound-specific fields
            "origin": self.origin.value,
            "safety_category": self.safety_category.value,
            "discovery_stage": self.metadata.discovery_stage.value,
            
            # Biological data
            "biological_activities": [
                activity.dict() for activity in self.biological_activities
            ],
            "admet_properties": self.admet.dict(),
            
            # Synthesis and IP
            "synthesis_info": self.synthesis.dict(),
            "patent_info": self.patent_info.dict(),
            
            # Scoring
            "scores": self.scores.dict(),
            
            # Metadata
            "therapeutic_areas": self.metadata.therapeutic_areas,
            "mechanism_of_action": self.metadata.mechanism_of_action,
            "compound_class": self.metadata.compound_class,
            "publications": self.metadata.publications,
            "clinical_trials": self.metadata.clinical_trials,
            "quality_flags": list(self.metadata.quality_flags),
            
            # Relationships
            "analogs": list(self.analogs),
            "parent_compound": self.parent_compound,
            "derivative_compounds": list(self.derivative_compounds)
        }
        
        return compound_dict
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Compound':
        """Create compound from dictionary representation"""
        smiles = data.get("smiles")
        if not smiles:
            raise CompoundError("SMILES required to create compound from dict")
        
        compound = cls(
            identifier=smiles,
            mol_id=data.get("mol_id"),
            name=data.get("name"),
            source=data.get("source"),
            origin=CompoundOrigin(data.get("origin", "unknown"))
        )
        
        # Set safety category
        if "safety_category" in data:
            compound.safety_category = SafetyCategory(data["safety_category"])
        
        # Load biological activities
        if "biological_activities" in data:
            for activity_data in data["biological_activities"]:
                activity = BiologicalActivity(**activity_data)
                compound.biological_activities.append(activity)
        
        # Load ADMET properties
        if "admet_properties" in data:
            compound.admet = ADMETProperties(**data["admet_properties"])
        
        # Load synthesis info
        if "synthesis_info" in data:
            compound.synthesis = SynthesisInfo(**data["synthesis_info"])
        
        # Load patent info
        if "patent_info" in data:
            compound.patent_info = PatentInfo(**data["patent_info"])
        
        # Load scores
        if "scores" in data:
            compound.scores = CompoundScores(**data["scores"])
        
        return compound


class CompoundLibrary(MoleculeCollection):
    """
    Specialized collection for compound discovery and management
    Extended functionality for drug discovery pipelines
    """
    
    def __init__(
        self,
        compounds: Optional[List[Compound]] = None,
        name: Optional[str] = None,
        description: Optional[str] = None
    ):
        """
        Initialize compound library
        
        Args:
            compounds: List of Compound objects
            name: Library name
            description: Library description
        """
        super().__init__()
        self.name = name
        self.description = description
        self.creation_date = datetime.now()
        
        # Replace molecule storage with compound storage
        self.compounds: Dict[str, Compound] = {}
        
        if compounds:
            for compound in compounds:
                self.add_compound(compound)
    
    def add_compound(self, compound: Compound) -> None:
        """Add compound to library"""
        self.compounds[compound.mol_id] = compound
        # Keep molecule dict in sync for compatibility
        self.molecules[compound.mol_id] = compound
    
    def get_compound(self, mol_id: str) -> Optional[Compound]:
        """Get compound by ID"""
        return self.compounds.get(mol_id)
    
    def filter_by_origin(self, origin: CompoundOrigin) -> 'CompoundLibrary':
        """Filter compounds by origin"""
        filtered = [comp for comp in self.compounds.values() if comp.origin == origin]
        return CompoundLibrary(filtered, name=f"{self.name}_filtered_by_{origin.value}")
    
    def filter_by_safety_category(self, category: SafetyCategory) -> 'CompoundLibrary':
        """Filter compounds by safety category"""
        filtered = [comp for comp in self.compounds.values() if comp.safety_category == category]
        return CompoundLibrary(filtered, name=f"{self.name}_safety_{category.value}")
    
    def filter_by_therapeutic_area(self, area: str) -> 'CompoundLibrary':
        """Filter compounds by therapeutic area"""
        filtered = [
            comp for comp in self.compounds.values()
            if area.lower() in [ta.lower() for ta in comp.metadata.therapeutic_areas]
        ]
        return CompoundLibrary(filtered, name=f"{self.name}_therapeutic_{area}")
    
    def get_hits(self, threshold: float = 10.0, activity_type: str = "IC50") -> 'CompoundLibrary':
        """Get all hit compounds"""
        hits = [comp for comp in self.compounds.values() if comp.is_hit(threshold, activity_type)]
        return CompoundLibrary(hits, name=f"{self.name}_hits")
    
    def get_leads(self) -> 'CompoundLibrary':
        """Get all lead compounds"""
        leads = [comp for comp in self.compounds.values() if comp.is_lead()]
        return CompoundLibrary(leads, name=f"{self.name}_leads")
    
    def rank_by_crowe_score(self, limit: Optional[int] = None) -> List[Compound]:
        """Rank compounds by Crowe Score"""
        compounds_with_scores = []
        
        for compound in self.compounds.values():
            if compound.scores.crowe_score is None:
                compound.calculate_crowe_score()
            compounds_with_scores.append(compound)
        
        ranked = sorted(compounds_with_scores, key=lambda x: x.scores.crowe_score or 0, reverse=True)
        
        if limit:
            ranked = ranked[:limit]
        
        return ranked
    
    def export_summary_report(self, filename: Union[str, Path]) -> None:
        """Export comprehensive library summary report"""
        filename = Path(filename)
        
        summary = {
            "library_info": {
                "name": self.name,
                "description": self.description,
                "creation_date": self.creation_date.isoformat(),
                "compound_count": len(self.compounds)
            },
            "statistics": {
                "origin_distribution": self._get_origin_distribution(),
                "safety_distribution": self._get_safety_distribution(),
                "therapeutic_areas": self._get_therapeutic_areas(),
                "discovery_stages": self._get_discovery_stages(),
                "drug_like_percentage": self._get_drug_like_percentage(),
                "hit_count": len(self.get_hits().compounds),
                "lead_count": len(self.get_leads().compounds)
            },
            "top_compounds": [
                comp.to_dict() for comp in self.rank_by_crowe_score(limit=10)
            ]
        }
        
        with open(filename, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
    
    def _get_origin_distribution(self) -> Dict[str, int]:
        """Get compound origin distribution"""
        distribution = {}
        for compound in self.compounds.values():
            origin = compound.origin.value
            distribution[origin] = distribution.get(origin, 0) + 1
        return distribution
    
    def _get_safety_distribution(self) -> Dict[str, int]:
        """Get safety category distribution"""
        distribution = {}
        for compound in self.compounds.values():
            safety = compound.safety_category.value
            distribution[safety] = distribution.get(safety, 0) + 1
        return distribution
    
    def _get_therapeutic_areas(self) -> Dict[str, int]:
        """Get therapeutic area distribution"""
        areas = {}
        for compound in self.compounds.values():
            for area in compound.metadata.therapeutic_areas:
                areas[area] = areas.get(area, 0) + 1
        return areas
    
    def _get_discovery_stages(self) -> Dict[str, int]:
        """Get discovery stage distribution"""
        stages = {}
        for compound in self.compounds.values():
            stage = compound.metadata.discovery_stage.value
            stages[stage] = stages.get(stage, 0) + 1
        return stages
    
    def _get_drug_like_percentage(self) -> float:
        """Get percentage of drug-like compounds"""
        if not self.compounds:
            return 0.0
        
        drug_like_count = sum(1 for comp in self.compounds.values() if comp.is_drug_like())
        return (drug_like_count / len(self.compounds)) * 100
    
    def __str__(self) -> str:
        return f"CompoundLibrary(name={self.name}, compounds={len(self.compounds)})"