"""
CriOS Data Schemas
Pydantic models for type-safe data validation and serialization
"""

from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
from pydantic import BaseModel, Field, validator, root_validator


class LogLevel(str, Enum):
    """Logging levels"""
    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"
    CRITICAL = "CRITICAL"


class OutputFormat(str, Enum):
    """Supported output formats"""
    CSV = "csv"
    SDF = "sdf"
    JSON = "json"
    PARQUET = "parquet"
    XLSX = "xlsx"


class TargetClass(str, Enum):
    """Supported target classes"""
    GPCR = "GPCR"
    KINASE = "kinase"
    PROTEASE = "protease"
    ION_CHANNEL = "ion_channel"
    ENZYME = "enzyme"
    OTHER = "other"


class EngineType(str, Enum):
    """Discovery engine types"""
    SYNTHETIC = "synthetic"
    NEURO = "neuro"
    NATURAL = "natural"


# Configuration Models
class SystemConfig(BaseModel):
    """System configuration"""
    name: str = "CriOS Discovery Engine"
    version: str = "1.0.0"
    environment: str = "production"
    seed: int = 42
    debug: bool = False
    log_level: LogLevel = LogLevel.INFO
    log_json: bool = False


class PerformanceConfig(BaseModel):
    """Performance settings"""
    num_workers: int = Field(4, ge=1, le=64)
    batch_size: int = Field(1000, ge=1)
    max_memory_gb: float = Field(8.0, gt=0)
    enable_ray: bool = False
    cache_size_mb: int = Field(512, ge=0)
    vectorize_ops: bool = True


class DrugLikeFilters(BaseModel):
    """Drug-likeness filter thresholds"""
    molecular_weight: Dict[str, float] = Field(default_factory=lambda: {"min": 150.0, "max": 500.0})
    logp: Dict[str, float] = Field(default_factory=lambda: {"min": -2.0, "max": 5.0})
    hbd: Dict[str, int] = Field(default_factory=lambda: {"max": 5})
    hba: Dict[str, int] = Field(default_factory=lambda: {"max": 10})
    tpsa: Dict[str, float] = Field(default_factory=lambda: {"min": 20.0, "max": 140.0})
    rotatable_bonds: Dict[str, int] = Field(default_factory=lambda: {"max": 10})
    aromatic_rings: Dict[str, int] = Field(default_factory=lambda: {"max": 4})
    heavy_atoms: Dict[str, int] = Field(default_factory=lambda: {"min": 10, "max": 70})


class CroweWeights(BaseModel):
    """Crowe scoring weights"""
    potency: float = Field(0.25, ge=0.0, le=1.0)
    selectivity: float = Field(0.20, ge=0.0, le=1.0)
    admet: float = Field(0.25, ge=0.0, le=1.0)
    synthesis: float = Field(0.20, ge=0.0, le=1.0)
    novelty: float = Field(0.10, ge=0.0, le=1.0)
    
    @root_validator
    def weights_must_sum_to_one(cls, values):
        total = sum(values.values())
        if abs(total - 1.0) > 0.01:
            raise ValueError(f"Weights must sum to 1.0, got {total:.3f}")
        return values


class Config(BaseModel):
    """Main configuration model"""
    system: SystemConfig = SystemConfig()
    performance: PerformanceConfig = PerformanceConfig()
    drug_like_filters: DrugLikeFilters = DrugLikeFilters()
    crowe_weights: CroweWeights = CroweWeights()


# Input/Output Models
class MoleculeInput(BaseModel):
    """Input molecule specification"""
    smiles: str = Field(..., min_length=1, description="SMILES string")
    id: Optional[str] = Field(None, description="Molecule identifier")
    name: Optional[str] = Field(None, description="Molecule name")
    properties: Optional[Dict[str, Any]] = Field(None, description="Additional properties")
    
    @validator('smiles')
    def validate_smiles_not_empty(cls, v):
        if not v.strip():
            raise ValueError("SMILES string cannot be empty")
        return v.strip()


class MoleculeOutput(BaseModel):
    """Output molecule with computed properties"""
    smiles: str
    id: str
    name: Optional[str] = None
    molecular_weight: Optional[float] = None
    logp: Optional[float] = None
    hbd: Optional[int] = None
    hba: Optional[int] = None
    tpsa: Optional[float] = None
    rotatable_bonds: Optional[int] = None
    aromatic_rings: Optional[int] = None
    drug_like: Optional[bool] = None
    properties: Optional[Dict[str, Any]] = None


class ScoringInput(BaseModel):
    """Input for compound scoring"""
    molecules: List[MoleculeInput] = Field(..., min_items=1)
    target_class: Optional[TargetClass] = None
    optimize_for: List[str] = Field(default_factory=lambda: ["potency", "selectivity", "admet"])
    weights: Optional[CroweWeights] = None
    include_explanations: bool = True


class ComponentScores(BaseModel):
    """Individual component scores"""
    potency: float = Field(..., ge=0.0, le=1.0)
    selectivity: float = Field(..., ge=0.0, le=1.0)
    admet: float = Field(..., ge=0.0, le=1.0)
    synthesis: float = Field(..., ge=0.0, le=1.0)
    novelty: float = Field(..., ge=0.0, le=1.0)


class ScoredMolecule(BaseModel):
    """Molecule with computed scores"""
    smiles: str
    id: str
    name: Optional[str] = None
    molecular_weight: Optional[float] = None
    target_class: Optional[TargetClass] = None
    component_scores: ComponentScores
    crowe_score: float = Field(..., ge=0.0, le=1.0)
    drug_like: bool
    passes_filters: bool
    confidence: float = Field(..., ge=0.0, le=1.0)
    explanation: Optional[Dict[str, Any]] = None


class ScoringOutput(BaseModel):
    """Output from scoring operation"""
    molecules: List[ScoredMolecule]
    metadata: Dict[str, Any] = Field(default_factory=dict)
    run_info: Dict[str, Any] = Field(default_factory=dict)


class DesignInput(BaseModel):
    """Input for compound design"""
    engine: EngineType = EngineType.SYNTHETIC
    target_class: TargetClass = TargetClass.GPCR
    optimize_for: List[str] = Field(default_factory=lambda: ["potency", "selectivity", "admet"])
    n_compounds: int = Field(100, ge=1, le=10000)
    top_n: int = Field(20, ge=1)
    seed_molecules: Optional[List[str]] = None
    constraints: Optional[Dict[str, Any]] = None


class DesignOutput(BaseModel):
    """Output from design operation"""
    molecules: List[ScoredMolecule]
    metadata: Dict[str, Any] = Field(default_factory=dict)
    run_info: Dict[str, Any] = Field(default_factory=dict)


class SimilarityInput(BaseModel):
    """Input for similarity calculation"""
    query_molecules: List[str] = Field(..., min_items=1, description="Query SMILES")
    reference_molecules: List[str] = Field(..., min_items=1, description="Reference SMILES")
    metric: str = Field("tanimoto", description="Similarity metric")
    fingerprint_type: str = Field("morgan", description="Fingerprint type")
    threshold: float = Field(0.5, ge=0.0, le=1.0, description="Similarity threshold")


class SimilarityPair(BaseModel):
    """Similarity between two molecules"""
    query_smiles: str
    reference_smiles: str
    similarity: float = Field(..., ge=0.0, le=1.0)
    query_id: Optional[str] = None
    reference_id: Optional[str] = None


class SimilarityOutput(BaseModel):
    """Output from similarity calculation"""
    pairs: List[SimilarityPair]
    metadata: Dict[str, Any] = Field(default_factory=dict)


class EthicsInput(BaseModel):
    """Input for ethics checking"""
    molecules: List[MoleculeInput] = Field(..., min_items=1)
    policy_file: Optional[str] = None
    include_explanations: bool = True


class EthicsViolation(BaseModel):
    """Ethics policy violation"""
    rule_id: str
    rule_name: str
    severity: str
    description: str
    triggered_by: str


class EthicsResult(BaseModel):
    """Ethics check result for a molecule"""
    smiles: str
    id: str
    passed: bool
    violations: List[EthicsViolation] = Field(default_factory=list)
    warnings: List[str] = Field(default_factory=list)
    confidence: float = Field(..., ge=0.0, le=1.0)


class EthicsOutput(BaseModel):
    """Output from ethics checking"""
    results: List[EthicsResult]
    summary: Dict[str, Any] = Field(default_factory=dict)
    metadata: Dict[str, Any] = Field(default_factory=dict)


# Manifest and Provenance Models
class EnvironmentInfo(BaseModel):
    """Runtime environment information"""
    python_version: str
    platform: str
    cpu_count: int
    memory_total_gb: float
    rdkit_version: Optional[str] = None
    numpy_version: Optional[str] = None
    pandas_version: Optional[str] = None


class RunManifest(BaseModel):
    """Complete run provenance information"""
    run_id: str = Field(..., description="Unique run identifier")
    timestamp: datetime = Field(default_factory=datetime.now)
    command: str = Field(..., description="Command executed")
    config_hash: str = Field(..., description="Configuration hash")
    input_hash: str = Field(..., description="Input data hash")
    environment: EnvironmentInfo = Field(..., description="Runtime environment")
    performance_metrics: Dict[str, Any] = Field(default_factory=dict)
    outputs: List[str] = Field(default_factory=list, description="Output file paths")
    success: bool = Field(True, description="Whether run completed successfully")
    error_message: Optional[str] = Field(None, description="Error message if failed")
    wall_time_seconds: Optional[float] = Field(None, description="Total execution time")


# Batch Processing Models
class BatchJob(BaseModel):
    """Batch processing job specification"""
    job_id: str = Field(..., description="Unique job identifier")
    operation: str = Field(..., description="Operation to perform")
    inputs: List[Dict[str, Any]] = Field(..., description="Input specifications")
    config: Dict[str, Any] = Field(default_factory=dict, description="Operation configuration")
    priority: int = Field(0, ge=0, le=10, description="Job priority (higher = more priority)")
    dependencies: List[str] = Field(default_factory=list, description="Dependent job IDs")


class BatchStatus(BaseModel):
    """Batch job status"""
    job_id: str
    status: str = Field(..., description="Job status: pending, running, completed, failed")
    progress: float = Field(0.0, ge=0.0, le=1.0, description="Progress fraction")
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    error_message: Optional[str] = None
    result_files: List[str] = Field(default_factory=list)


# API Models
class APIResponse(BaseModel):
    """Standard API response wrapper"""
    success: bool = True
    data: Optional[Dict[str, Any]] = None
    message: Optional[str] = None
    error: Optional[str] = None
    timestamp: datetime = Field(default_factory=datetime.now)
    request_id: Optional[str] = None


class HealthCheck(BaseModel):
    """Health check response"""
    status: str = "healthy"
    version: str = "1.0.0"
    timestamp: datetime = Field(default_factory=datetime.now)
    uptime_seconds: Optional[float] = None
    memory_usage_mb: Optional[float] = None
    active_jobs: int = 0