"""
CriOS Discovery Engine - Configuration and Policy Schemas
Pydantic models for type-safe configuration and ethical policy validation
"""

from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Union, Any
from pydantic import BaseModel, Field, validator, root_validator


# Enumerations
class EnvironmentType(str, Enum):
    DEVELOPMENT = "development"
    TESTING = "testing"
    PRODUCTION = "production"


class LogLevel(str, Enum):
    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"
    CRITICAL = "CRITICAL"


class FingerprintType(str, Enum):
    MORGAN = "morgan"
    RDKIT = "rdkit"
    MACCS = "maccs"


class ClusteringAlgorithm(str, Enum):
    BUTINA = "butina"
    WARD = "ward"
    KMEANS = "kmeans"


class OutputFormat(str, Enum):
    JSON = "json"
    CSV = "csv"
    SDF = "sdf"
    XLSX = "xlsx"


class SafetyCategory(str, Enum):
    GREEN = "green"
    YELLOW = "yellow"
    RED = "red"


# System Configuration Models
class SystemConfig(BaseModel):
    name: str = "CriOS Discovery Engine"
    version: str = "0.1.0"
    environment: EnvironmentType = EnvironmentType.DEVELOPMENT
    log_level: LogLevel = LogLevel.INFO
    debug: bool = True
    seed: int = 42


class DatabaseConfig(BaseModel):
    engine: str = "sqlite"
    path: str = "data/crios.db"
    backup_enabled: bool = True
    backup_interval_hours: int = 24
    connection_pool_size: int = 20


class APIConfig(BaseModel):
    base_url: str
    timeout: int = 30
    rate_limit_requests: int = 100
    rate_limit_window: int = 3600


class ExternalAPIsConfig(BaseModel):
    chembl: APIConfig
    pubchem: APIConfig
    rdkit: Dict[str, Union[int, bool]]


# Discovery Configuration Models
class NaturalProductsConfig(BaseModel):
    enabled: bool = True
    sources: List[str] = ["fungi", "plants", "marine", "microbial"]
    bioactivity_filters: List[str] = ["neuroprotective", "antioxidant", "anti-inflammatory"]
    scaffold_extraction: bool = True
    pharmacophore_analysis: bool = True


class SyntheticChemistryConfig(BaseModel):
    enabled: bool = True
    target_classes: List[str] = ["GPCR", "kinase", "protease", "ion_channel"]
    drug_like_filters: bool = True
    lipinski_compliant: bool = True
    synthetic_accessibility: bool = True


class HybridDesignConfig(BaseModel):
    enabled: bool = True
    natural_scaffold_optimization: bool = True
    synthetic_modification: bool = True
    bioisosterism: bool = True


class DiscoveryConfig(BaseModel):
    natural_products: NaturalProductsConfig
    synthetic_chemistry: SyntheticChemistryConfig
    hybrid_design: HybridDesignConfig


# Molecular Processing Models
class SimilarityConfig(BaseModel):
    default_threshold: float = Field(0.7, ge=0.0, le=1.0)
    fingerprint_type: FingerprintType = FingerprintType.MORGAN
    fingerprint_radius: int = Field(2, ge=1, le=4)
    fingerprint_bits: int = Field(2048, ge=512, le=8192)
    batch_size: int = Field(1000, ge=1)


class ClusteringConfig(BaseModel):
    algorithm: ClusteringAlgorithm = ClusteringAlgorithm.BUTINA
    distance_threshold: float = Field(0.4, ge=0.0, le=1.0)
    max_clusters: int = Field(100, ge=1)
    min_cluster_size: int = Field(2, ge=1)


class DescriptorsConfig(BaseModel):
    molecular_weight: bool = True
    logp: bool = True
    tpsa: bool = True
    hbd: bool = True
    hba: bool = True
    rotatable_bonds: bool = True
    aromatic_rings: bool = True


class ADMETConfig(BaseModel):
    bbb_permeability: bool = True
    oral_bioavailability: bool = True
    cyp_inhibition: bool = True
    hepatotoxicity: bool = True
    cardiotoxicity: bool = True


class MolecularConfig(BaseModel):
    similarity: SimilarityConfig
    clustering: ClusteringConfig
    descriptors: DescriptorsConfig
    admet: ADMETConfig


# Scoring Configuration Models
class CroweScoreWeights(BaseModel):
    novelty: float = Field(0.25, ge=0.0, le=1.0)
    drug_likeness: float = Field(0.20, ge=0.0, le=1.0)
    synthetic_accessibility: float = Field(0.15, ge=0.0, le=1.0)
    safety_profile: float = Field(0.20, ge=0.0, le=1.0)
    ethical_alignment: float = Field(0.20, ge=0.0, le=1.0)
    
    @root_validator
    def weights_sum_to_one(cls, values):
        total = sum(values.values())
        if abs(total - 1.0) > 0.01:
            raise ValueError(f"Weights must sum to 1.0, got {total}")
        return values


class CroweScoreConfig(BaseModel):
    enabled: bool = True
    weights: CroweScoreWeights
    normalization: str = "min_max"


class OptimizationConfig(BaseModel):
    objectives: List[str] = ["efficacy", "safety", "druggability"]
    pareto_ranking: bool = True
    diversity_penalty: float = Field(0.1, ge=0.0, le=1.0)


class ScoringConfig(BaseModel):
    crowe_score: CroweScoreConfig
    optimization: OptimizationConfig


# Performance Configuration Models
class ParallelConfig(BaseModel):
    enabled: bool = True
    max_workers: int = Field(8, ge=1)
    chunk_size: int = Field(100, ge=1)


class CacheConfig(BaseModel):
    enabled: bool = True
    ttl_seconds: int = Field(3600, ge=0)
    max_size_mb: int = Field(1024, ge=1)


class MemoryConfig(BaseModel):
    max_compounds_in_memory: int = Field(10000, ge=1)
    garbage_collection_threshold: float = Field(0.8, ge=0.0, le=1.0)


class PerformanceConfig(BaseModel):
    parallel: ParallelConfig
    cache: CacheConfig
    memory: MemoryConfig


# Interface Configuration Models
class WebConfig(BaseModel):
    host: str = "0.0.0.0"
    port: int = Field(8000, ge=1024, le=65535)
    reload: bool = True
    cors_origins: List[str] = ["http://localhost:3000", "https://app.crios.ai"]
    max_upload_size_mb: int = Field(100, ge=1)


class CLIConfig(BaseModel):
    progress_bars: bool = True
    colored_output: bool = True
    verbose: bool = False
    output_format: OutputFormat = OutputFormat.JSON


# Quality Control Models
class ValidationRules(BaseModel):
    max_molecular_weight: float = Field(1000.0, gt=0)
    min_molecular_weight: float = Field(100.0, gt=0)
    max_heavy_atoms: int = Field(50, gt=0)
    min_heavy_atoms: int = Field(5, gt=0)
    
    @validator('min_molecular_weight')
    def min_weight_less_than_max(cls, v, values):
        if 'max_molecular_weight' in values and v >= values['max_molecular_weight']:
            raise ValueError('min_molecular_weight must be less than max_molecular_weight')
        return v


class FilterConfig(BaseModel):
    pains: bool = True
    brenk: bool = True
    nih: bool = True
    zinc: bool = True


class QualityConfig(BaseModel):
    validation: ValidationRules
    filters: FilterConfig


class ExportConfig(BaseModel):
    formats: List[OutputFormat] = [OutputFormat.SDF, OutputFormat.CSV, OutputFormat.JSON, OutputFormat.XLSX]
    default_format: OutputFormat = OutputFormat.SDF
    include_descriptors: bool = True
    include_scores: bool = True


# Main Configuration Model
class CriosConfig(BaseModel):
    """Complete CriOS Discovery Engine Configuration"""
    system: SystemConfig
    database: DatabaseConfig
    apis: ExternalAPIsConfig
    discovery: DiscoveryConfig
    molecular: MolecularConfig
    scoring: ScoringConfig
    performance: PerformanceConfig
    web: WebConfig
    cli: CLIConfig
    quality: QualityConfig
    export: ExportConfig
    
    class Config:
        use_enum_values = True
        validate_assignment = True


# Ethical Policy Models
class ResearchPriorities(BaseModel):
    primary: List[str] = [
        "neurodegeneration_therapeutics",
        "environmental_health", 
        "regenerative_medicine",
        "sustainable_materials"
    ]
    secondary: List[str] = [
        "aging_research",
        "metabolic_disorders", 
        "natural_product_optimization"
    ]


class DiscoveryEthics(BaseModel):
    research_priorities: ResearchPriorities
    excluded_applications: List[str] = [
        "biological_weapons",
        "environmental_damage",
        "addictive_substances",
        "predatory_pricing_models"
    ]


class HumanOversightConfig(BaseModel):
    required_for: List[str] = [
        "novel_compound_generation",
        "safety_predictions",
        "therapeutic_recommendations",
        "patent_analysis"
    ]
    validation_threshold: float = Field(0.95, ge=0.0, le=1.0)
    expert_review_required: bool = True


class BiasMitigationConfig(BaseModel):
    training_data_diversity: bool = True
    fairness_metrics: List[str] = ["demographic_parity", "equalized_odds"]
    bias_auditing_frequency: str = "monthly"


class TransparencyConfig(BaseModel):
    model_interpretability: str = "required"
    decision_explanations: str = "detailed"
    uncertainty_quantification: bool = True
    confidence_intervals: bool = True


class AISafetyConfig(BaseModel):
    human_oversight: HumanOversightConfig
    bias_mitigation: BiasMitigationConfig
    transparency: TransparencyConfig


class DataSourcesConfig(BaseModel):
    public_databases: List[str] = ["chembl", "pubchem", "drugbank", "zinc"]
    literature: List[str] = ["pubmed", "patent_databases", "scientific_journals"]


class DataUsagePrinciples(BaseModel):
    attribution_required: bool = True
    open_access_preferred: bool = True
    commercial_restrictions_respected: bool = True
    privacy_preserving: bool = True


class DataEthicsConfig(BaseModel):
    approved_sources: DataSourcesConfig
    prohibited_sources: List[str] = [
        "proprietary_without_permission",
        "personally_identifiable_data",
        "military_classified_data"
    ]
    usage_principles: DataUsagePrinciples


class CompoundSafetyCategories(BaseModel):
    green: str = "Safe for general research and development"
    yellow: str = "Requires additional safety assessment" 
    red: str = "Restricted or prohibited compounds"


class DualUseConfig(BaseModel):
    screening_required: bool = True
    review_committee: str = "CriOS Ethics Board"
    restricted_compounds: List[str] = [
        "chemical_weapons_precursors",
        "controlled_substances",
        "environmental_toxins"
    ]


class NaturalProductEthics(BaseModel):
    indigenous_knowledge_respect: bool = True
    biopiracy_prevention: bool = True
    fair_benefit_sharing: bool = True
    source_attribution: str = "required"


class CompoundEthicsConfig(BaseModel):
    safety_categories: CompoundSafetyCategories
    dual_use: DualUseConfig
    natural_products: NaturalProductEthics


class MissionLockConfig(BaseModel):
    enabled: bool = True
    ethical_licensing_required: bool = True
    public_benefit_clause: bool = True
    profit_cap_percentage: int = Field(15, ge=0, le=100)


class PatentConfig(BaseModel):
    defensive_only: bool = False
    open_licensing_preferred: bool = True
    fair_licensing_terms: bool = True
    research_exemptions: bool = True


class PublicationConfig(BaseModel):
    preprint_encouraged: bool = True
    open_access_preferred: bool = True
    methodology_transparency: str = "complete"
    negative_results_published: bool = True


class IPEthicsConfig(BaseModel):
    mission_lock: MissionLockConfig
    patents: PatentConfig
    publication: PublicationConfig


class EnvironmentalImpact(BaseModel):
    lifecycle_assessment: bool = True
    biodegradability_preferred: bool = True
    green_chemistry_principles: bool = True
    carbon_footprint_minimization: bool = True


class SocialImpact(BaseModel):
    accessibility_prioritized: bool = True
    global_south_considerations: bool = True
    healthcare_equity: bool = True
    community_benefit_focus: bool = True


class ImpactAssessmentConfig(BaseModel):
    environmental: EnvironmentalImpact
    social: SocialImpact


class ReviewBoardConfig(BaseModel):
    composition: List[str] = ["ethicist", "scientist", "community_representative"]
    review_frequency: str = "quarterly"
    decision_authority: str = "binding"
    appeal_process: str = "available"


class ComplianceConfig(BaseModel):
    automated_screening: bool = True
    human_review_required: bool = True
    violation_reporting: str = "mandatory"
    corrective_actions: str = "documented"


class GovernanceConfig(BaseModel):
    review_board: ReviewBoardConfig
    compliance: ComplianceConfig


# Main Ethics Policy Model
class EthicsPolicy(BaseModel):
    """Complete CriOS Discovery Engine Ethics Policy"""
    principles: Dict[str, Union[str, List[str]]]
    discovery_ethics: DiscoveryEthics
    ai_safety: AISafetyConfig
    data_ethics: DataEthicsConfig
    compound_ethics: CompoundEthicsConfig
    ip_ethics: IPEthicsConfig
    impact_assessment: ImpactAssessmentConfig
    governance: GovernanceConfig
    
    class Config:
        use_enum_values = True
        validate_assignment = True


# Configuration Management
class ConfigManager:
    """Configuration and Ethics Policy Manager"""
    
    def __init__(self, config_dir: Path = Path("config")):
        self.config_dir = config_dir
        self._config: Optional[CriosConfig] = None
        self._ethics: Optional[EthicsPolicy] = None
    
    def load_config(self, config_file: str = "config.yaml") -> CriosConfig:
        """Load configuration from YAML file"""
        import yaml
        
        config_path = self.config_dir / config_file
        with open(config_path, 'r') as f:
            config_data = yaml.safe_load(f)
        
        self._config = CriosConfig(**config_data)
        return self._config
    
    def load_ethics(self, ethics_file: str = "ethics.yaml") -> EthicsPolicy:
        """Load ethics policy from YAML file"""
        import yaml
        
        ethics_path = self.config_dir / ethics_file
        with open(ethics_path, 'r') as f:
            ethics_data = yaml.safe_load(f)
        
        self._ethics = EthicsPolicy(**ethics_data)
        return self._ethics
    
    @property
    def config(self) -> CriosConfig:
        """Get current configuration"""
        if self._config is None:
            self._config = self.load_config()
        return self._config
    
    @property
    def ethics(self) -> EthicsPolicy:
        """Get current ethics policy"""
        if self._ethics is None:
            self._ethics = self.load_ethics()
        return self._ethics
    
    def validate_ethical_compliance(self, compound_data: Dict[str, Any]) -> Dict[str, Any]:
        """Validate compound against ethical policies"""
        results = {
            "compliant": True,
            "violations": [],
            "warnings": [],
            "recommendations": []
        }
        
        # Check safety category
        if "safety_category" in compound_data:
            if compound_data["safety_category"] == SafetyCategory.RED:
                results["compliant"] = False
                results["violations"].append("Compound classified as restricted (RED category)")
        
        # Check dual-use concerns
        if "compound_class" in compound_data:
            restricted = self.ethics.compound_ethics.dual_use.restricted_compounds
            if any(cls in compound_data["compound_class"].lower() for cls in restricted):
                results["compliant"] = False
                results["violations"].append("Compound flagged for dual-use restrictions")
        
        # Check environmental impact
        if not self.ethics.impact_assessment.environmental.green_chemistry_principles:
            results["warnings"].append("Green chemistry principles not enforced")
        
        return results


# Export key models
__all__ = [
    "CriosConfig", 
    "EthicsPolicy", 
    "ConfigManager",
    "SafetyCategory",
    "EnvironmentType",
    "LogLevel",
    "FingerprintType",
    "ClusteringAlgorithm", 
    "OutputFormat"
]