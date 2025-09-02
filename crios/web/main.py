"""
CriOS Web API - Main Application
FastAPI backend for compound discovery and analysis
"""

import logging
from contextlib import asynccontextmanager
from typing import Dict, List, Optional

from fastapi import FastAPI, HTTPException, Depends
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from ..core.molecule import Molecule, validate_smiles
from ..core.similarity import SimilaritySearch

logger = logging.getLogger(__name__)

# Global instances
similarity_search: Optional[SimilaritySearch] = None


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Initialize services on startup"""
    global similarity_search
    
    try:
        # Initialize core services
        similarity_search = SimilaritySearch()
        
        logger.info("CriOS web services initialized successfully")
        
    except Exception as e:
        logger.error(f"Failed to initialize services: {e}")
        raise
    
    yield
    
    # Cleanup
    logger.info("CriOS web services shutting down")


# Create FastAPI app
app = FastAPI(
    title="CriOS Discovery Engine API",
    description="Universal Compound Discovery Platform",
    version="0.1.0",
    docs_url="/docs",
    redoc_url="/redoc",
    lifespan=lifespan
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure appropriately for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)


# Pydantic models for API
class MoleculeValidationRequest(BaseModel):
    smiles: str
    standardize: bool = False


class MoleculeValidationResponse(BaseModel):
    valid: bool
    smiles: str
    canonical_smiles: Optional[str] = None
    error: Optional[str] = None
    molecular_weight: Optional[float] = None
    drug_like: Optional[bool] = None


class SimilaritySearchRequest(BaseModel):
    query_smiles: str
    database_smiles: List[str]
    threshold: float = 0.7
    max_results: int = 10
    fingerprint_type: str = "morgan"


class SimilarityResult(BaseModel):
    smiles: str
    similarity: float
    rank: int


class SimilaritySearchResponse(BaseModel):
    query_smiles: str
    results: List[SimilarityResult]
    total_searched: int


class CroweScoreRequest(BaseModel):
    smiles: str
    compound_origin: str = "unknown"
    therapeutic_areas: List[str] = []


class CroweScoreResponse(BaseModel):
    smiles: str
    crowe_score: float
    components: Dict[str, float]
    confidence: float


# Health check endpoint
@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "service": "CriOS Discovery Engine",
        "version": "0.1.0"
    }


# Molecule validation endpoint
@app.post("/validate", response_model=MoleculeValidationResponse)
async def validate_molecule(request: MoleculeValidationRequest):
    """
    Validate and analyze a SMILES string
    """
    try:
        smiles = request.smiles.strip()
        
        if not smiles:
            raise HTTPException(status_code=400, detail="Empty SMILES string")
        
        # Basic validation
        is_valid = validate_smiles(smiles)
        
        if not is_valid:
            return MoleculeValidationResponse(
                valid=False,
                smiles=smiles,
                error="Invalid SMILES syntax"
            )
        
        # Create molecule for detailed analysis
        molecule = Molecule(smiles)
        
        response = MoleculeValidationResponse(
            valid=molecule.is_valid(),
            smiles=smiles,
            molecular_weight=molecule.properties.molecular_weight,
            drug_like=molecule.is_drug_like()
        )
        
        # Add canonical SMILES if requested
        if request.standardize and molecule.is_valid():
            response.canonical_smiles = molecule.smiles
        
        if not molecule.is_valid():
            response.error = "Molecule could not be processed"
        
        return response
        
    except Exception as e:
        logger.error(f"Validation failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# Similarity search endpoint  
@app.post("/similarity", response_model=SimilaritySearchResponse)
async def similarity_search_endpoint(request: SimilaritySearchRequest):
    """
    Perform molecular similarity search
    """
    try:
        if not similarity_search:
            raise HTTPException(status_code=503, detail="Similarity search not initialized")
        
        # Create query molecule
        query_mol = Molecule(request.query_smiles, mol_id="query")
        
        if not query_mol.is_valid():
            raise HTTPException(status_code=400, detail="Invalid query SMILES")
        
        # Create target molecules
        target_molecules = []
        for i, smiles in enumerate(request.database_smiles):
            try:
                mol = Molecule(smiles, mol_id=f"target_{i}")
                if mol.is_valid():
                    target_molecules.append(mol)
            except Exception as e:
                logger.warning(f"Skipping invalid SMILES {smiles}: {e}")
        
        if not target_molecules:
            raise HTTPException(status_code=400, detail="No valid target molecules")
        
        # Create collection and perform search
        from ..core.molecule import MoleculeCollection
        target_collection = MoleculeCollection(target_molecules)
        
        results = similarity_search.search_single(
            query_mol, target_collection, 
            request.threshold, request.max_results
        )
        
        # Format response
        similarity_results = [
            SimilarityResult(
                smiles=result.target_molecule.smiles,
                similarity=result.similarity,
                rank=i + 1
            )
            for i, result in enumerate(results)
        ]
        
        return SimilaritySearchResponse(
            query_smiles=request.query_smiles,
            results=similarity_results,
            total_searched=len(target_molecules)
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Similarity search failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# Crowe score endpoint
@app.post("/crowe-score", response_model=CroweScoreResponse)
async def calculate_crowe_score(request: CroweScoreRequest):
    """
    Calculate Crowe Discovery Framework score for a compound
    """
    try:
        # Create basic molecule for validation
        molecule = Molecule(request.smiles, mol_id="crowe_candidate")
        
        if not molecule.is_valid():
            raise HTTPException(status_code=400, detail="Invalid SMILES")
        
        # Use simplified scoring from design module
        from ..scoring.design import crowe_score
        
        # Map therapeutic areas to target class
        target_class = "GPCR"  # Default
        if "neurotherapeutic" in request.therapeutic_areas:
            target_class = "GPCR"
        elif "synthetic" in request.therapeutic_areas:
            target_class = "Enzyme"
        
        scores = crowe_score(request.smiles, target_class=target_class)
        
        return CroweScoreResponse(
            smiles=request.smiles,
            crowe_score=scores.get('crowe_score', 0.0),
            components={
                "potency": scores.get('potency', 0.0),
                "selectivity": scores.get('selectivity', 0.0), 
                "admet": scores.get('admet', 0.0),
                "synthesis": scores.get('synthesis', 0.0),
                "safety": scores.get('safety', 0.0)
            },
            confidence=0.85  # Default confidence
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Crowe score calculation failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# Configuration endpoint
@app.get("/config")
async def get_config():
    """Get current configuration"""
    return {
        "system": {
            "name": "CriOS Discovery Engine",
            "version": "1.0.0",
            "environment": "development"
        },
        "discovery": {
            "natural_products_enabled": True,
            "synthetic_chemistry_enabled": True,
            "hybrid_design_enabled": True
        },
        "scoring": {
            "crowe_framework_enabled": True,
            "ethics_checking_enabled": True
        }
    }


class EthicsCheckRequest(BaseModel):
    smiles: str

# Ethics check endpoint
@app.post("/ethics-check")
async def ethics_check(request: EthicsCheckRequest):
    """Check compound against ethical policies"""
    try:
        # Basic ethics check using filters module
        from ..ethics.filters import check_ethics_compliance
        
        # Validate molecule first
        molecule = Molecule(request.smiles, mol_id="ethics_check")
        if not molecule.is_valid():
            raise HTTPException(status_code=400, detail="Invalid SMILES")
        
        # Perform ethics check
        result = check_ethics_compliance(request.smiles)
        
        return {
            "compliant": result["passed"],
            "violations": result["violations"],
            "warnings": [],
            "recommendations": [result["explanation"].get("recommendation", "No recommendations")]
        }
        
    except Exception as e:
        logger.error(f"Ethics check failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# Root endpoint
@app.get("/")
async def root():
    """Root endpoint with API information"""
    return {
        "message": "🧬 CriOS Discovery Engine API",
        "tagline": "Science before status. Discovery before profit. Earth and people above extraction.",
        "version": "0.1.0",
        "docs": "/docs",
        "health": "/health"
    }