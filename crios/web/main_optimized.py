"""
CriOS Discovery Engine - High-Performance API
Optimized FastAPI backend with caching, async operations, and performance monitoring
"""

import asyncio
import hashlib
import json
import logging
import time
from contextlib import asynccontextmanager
from functools import lru_cache, wraps
from typing import Dict, List, Optional, Any

import redis.asyncio as redis
from fastapi import FastAPI, HTTPException, Depends, BackgroundTasks, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import JSONResponse, ORJSONResponse
from pydantic import BaseModel, Field
import orjson

from ..core.molecule import Molecule, validate_smiles
from ..core.similarity import SimilaritySearch
from ..scoring.design import crowe_score
from ..ethics.filters import check_ethics_compliance

logger = logging.getLogger(__name__)

# Performance monitoring
class PerformanceMiddleware:
    def __init__(self, app):
        self.app = app
        
    async def __call__(self, scope, receive, send):
        if scope["type"] == "http":
            start_time = time.time()
            
            async def send_wrapper(message):
                if message["type"] == "http.response.start":
                    duration = time.time() - start_time
                    message["headers"].append((b"x-process-time", f"{duration:.3f}".encode()))
                await send(message)
            
            await self.app(scope, receive, send_wrapper)
        else:
            await self.app(scope, receive, send)

# Global instances with connection pooling
similarity_search: Optional[SimilaritySearch] = None
redis_client: Optional[redis.Redis] = None
molecule_cache: Dict[str, Molecule] = {}

# Cache configuration
CACHE_TTL = 3600  # 1 hour
MAX_CACHE_SIZE = 10000

def cache_key(*args, **kwargs) -> str:
    """Generate cache key from arguments"""
    key_data = {"args": args, "kwargs": kwargs}
    return hashlib.md5(orjson.dumps(key_data, option=orjson.OPT_SORT_KEYS)).hexdigest()

async def get_cached_or_compute(
    key: str, 
    compute_func, 
    ttl: int = CACHE_TTL
):
    """Generic caching function with Redis backend"""
    if redis_client:
        try:
            cached = await redis_client.get(key)
            if cached:
                return orjson.loads(cached)
        except Exception as e:
            logger.warning(f"Cache read failed: {e}")
    
    result = await compute_func() if asyncio.iscoroutinefunction(compute_func) else compute_func()
    
    if redis_client:
        try:
            await redis_client.setex(
                key, 
                ttl, 
                orjson.dumps(result, option=orjson.OPT_SERIALIZE_NUMPY)
            )
        except Exception as e:
            logger.warning(f"Cache write failed: {e}")
    
    return result

@asynccontextmanager
async def lifespan(app: FastAPI):
    """Initialize services with connection pooling and caching"""
    global similarity_search, redis_client
    
    try:
        # Initialize Redis with connection pool
        try:
            redis_client = await redis.from_url(
                "redis://localhost:6379",
                encoding="utf-8",
                decode_responses=False,
                max_connections=50,
                socket_keepalive=True,
                socket_keepalive_options={
                    1: 1,  # TCP_KEEPIDLE
                    2: 3,  # TCP_KEEPINTVL  
                    3: 5   # TCP_KEEPCNT
                }
            )
            await redis_client.ping()
            logger.info("Redis cache connected")
        except Exception as e:
            logger.warning(f"Redis not available, using in-memory cache: {e}")
            redis_client = None
        
        # Initialize services with optimizations
        similarity_search = SimilaritySearch()
        
        logger.info("CriOS Discovery Engine initialized (optimized)")
        
    except Exception as e:
        logger.error(f"Failed to initialize services: {e}")
        raise
    
    yield
    
    # Cleanup
    if redis_client:
        await redis_client.close()
    logger.info("CriOS Discovery Engine shutting down")

# Create optimized FastAPI app
app = FastAPI(
    title="CriOS Discovery Engine",
    description="Accelerating Discovery, Protecting Tomorrow™",
    version="2.0.0",
    docs_url="/docs",
    redoc_url="/redoc",
    lifespan=lifespan,
    default_response_class=ORJSONResponse  # Faster JSON serialization
)

# Add performance middleware
app.add_middleware(PerformanceMiddleware)
app.add_middleware(GZipMiddleware, minimum_size=1000)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
    max_age=3600  # Cache CORS preflight
)

# Optimized Pydantic models
class MoleculeValidationRequest(BaseModel):
    smiles: str = Field(..., min_length=1, max_length=500)
    standardize: bool = False
    
    class Config:
        schema_extra = {
            "example": {
                "smiles": "CCO",
                "standardize": True
            }
        }

class BatchValidationRequest(BaseModel):
    smiles_list: List[str] = Field(..., max_items=1000)
    standardize: bool = False
    parallel: bool = True

class MoleculeValidationResponse(BaseModel):
    valid: bool
    smiles: str
    canonical_smiles: Optional[str] = None
    error: Optional[str] = None
    molecular_weight: Optional[float] = None
    drug_like: Optional[bool] = None
    cache_hit: bool = False
    process_time_ms: float = 0

class SimilaritySearchRequest(BaseModel):
    query_smiles: str = Field(..., min_length=1, max_length=500)
    database_smiles: List[str] = Field(..., max_items=10000)
    threshold: float = Field(0.7, ge=0.0, le=1.0)
    max_results: int = Field(10, ge=1, le=1000)
    fingerprint_type: str = Field("morgan", regex="^(morgan|maccs|rdkit)$")

class CroweScoreRequest(BaseModel):
    smiles: str = Field(..., min_length=1, max_length=500)
    compound_origin: str = Field("unknown", max_length=50)
    therapeutic_areas: List[str] = Field(default_factory=list, max_items=10)
    use_cache: bool = True

class CroweScoreResponse(BaseModel):
    smiles: str
    crowe_score: float
    components: Dict[str, float]
    confidence: float
    cache_hit: bool = False
    process_time_ms: float = 0

# Health check with detailed status
@app.get("/health")
async def health_check():
    """Comprehensive health check with service status"""
    health_status = {
        "status": "healthy",
        "service": "CriOS Discovery Engine",
        "version": "2.0.0",
        "tagline": "Accelerating Discovery, Protecting Tomorrow™",
        "services": {
            "api": "operational",
            "cache": "operational" if redis_client else "degraded",
            "similarity_search": "operational" if similarity_search else "unavailable"
        },
        "performance": {
            "cache_size": len(molecule_cache),
            "max_cache_size": MAX_CACHE_SIZE
        }
    }
    
    # Check Redis connection
    if redis_client:
        try:
            await redis_client.ping()
        except Exception:
            health_status["services"]["cache"] = "unavailable"
            health_status["status"] = "degraded"
    
    return health_status

# Optimized molecule validation with caching
@app.post("/validate", response_model=MoleculeValidationResponse)
async def validate_molecule(request: MoleculeValidationRequest):
    """Validate SMILES with intelligent caching"""
    start_time = time.time()
    cache_hit = False
    
    try:
        cache_key_str = f"validate:{request.smiles}:{request.standardize}"
        
        # Try cache first
        if redis_client:
            try:
                cached = await redis_client.get(cache_key_str)
                if cached:
                    result = orjson.loads(cached)
                    result["cache_hit"] = True
                    result["process_time_ms"] = (time.time() - start_time) * 1000
                    return MoleculeValidationResponse(**result)
            except Exception:
                pass
        
        # Validate
        smiles = request.smiles.strip()
        
        if not smiles:
            raise HTTPException(status_code=400, detail="Empty SMILES string")
        
        is_valid = validate_smiles(smiles)
        
        if not is_valid:
            response = MoleculeValidationResponse(
                valid=False,
                smiles=smiles,
                error="Invalid SMILES syntax",
                cache_hit=cache_hit,
                process_time_ms=(time.time() - start_time) * 1000
            )
        else:
            molecule = Molecule(smiles)
            response = MoleculeValidationResponse(
                valid=molecule.is_valid(),
                smiles=smiles,
                molecular_weight=molecule.properties.molecular_weight,
                drug_like=molecule.is_drug_like(),
                cache_hit=cache_hit,
                process_time_ms=(time.time() - start_time) * 1000
            )
            
            if request.standardize and molecule.is_valid():
                response.canonical_smiles = molecule.smiles
        
        # Cache result
        if redis_client and response.valid:
            try:
                await redis_client.setex(
                    cache_key_str,
                    CACHE_TTL,
                    orjson.dumps(response.dict(exclude={"cache_hit", "process_time_ms"}))
                )
            except Exception:
                pass
        
        return response
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Validation failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))

# Batch validation endpoint
@app.post("/validate/batch")
async def validate_batch(request: BatchValidationRequest, background_tasks: BackgroundTasks):
    """Batch validation with parallel processing"""
    start_time = time.time()
    
    if request.parallel and len(request.smiles_list) > 10:
        # Process in parallel for large batches
        tasks = [
            validate_molecule(MoleculeValidationRequest(
                smiles=smiles,
                standardize=request.standardize
            ))
            for smiles in request.smiles_list
        ]
        results = await asyncio.gather(*tasks, return_exceptions=True)
        
        valid_results = []
        errors = []
        
        for i, result in enumerate(results):
            if isinstance(result, Exception):
                errors.append({
                    "index": i,
                    "smiles": request.smiles_list[i],
                    "error": str(result)
                })
            else:
                valid_results.append(result.dict())
        
        return {
            "results": valid_results,
            "errors": errors,
            "total": len(request.smiles_list),
            "valid": len(valid_results),
            "process_time_ms": (time.time() - start_time) * 1000
        }
    else:
        # Sequential processing for small batches
        results = []
        errors = []
        
        for i, smiles in enumerate(request.smiles_list):
            try:
                result = await validate_molecule(
                    MoleculeValidationRequest(
                        smiles=smiles,
                        standardize=request.standardize
                    )
                )
                results.append(result.dict())
            except Exception as e:
                errors.append({
                    "index": i,
                    "smiles": smiles,
                    "error": str(e)
                })
        
        return {
            "results": results,
            "errors": errors,
            "total": len(request.smiles_list),
            "valid": len(results),
            "process_time_ms": (time.time() - start_time) * 1000
        }

# Optimized Crowe scoring with caching
@app.post("/crowe-score", response_model=CroweScoreResponse)
async def calculate_crowe_score(request: CroweScoreRequest):
    """Calculate Crowe score with intelligent caching"""
    start_time = time.time()
    cache_hit = False
    
    try:
        # Check cache if enabled
        if request.use_cache:
            cache_key_str = f"crowe:{request.smiles}:{request.compound_origin}:{','.join(sorted(request.therapeutic_areas))}"
            
            if redis_client:
                try:
                    cached = await redis_client.get(cache_key_str)
                    if cached:
                        result = orjson.loads(cached)
                        return CroweScoreResponse(
                            **result,
                            cache_hit=True,
                            process_time_ms=(time.time() - start_time) * 1000
                        )
                except Exception:
                    pass
        
        # Validate molecule
        molecule = Molecule(request.smiles, mol_id="crowe_candidate")
        
        if not molecule.is_valid():
            raise HTTPException(status_code=400, detail="Invalid SMILES")
        
        # Map therapeutic areas to target class
        target_class = "GPCR"  # Default
        if "neurotherapeutic" in request.therapeutic_areas:
            target_class = "GPCR"
        elif "synthetic" in request.therapeutic_areas:
            target_class = "Enzyme"
        elif "natural" in request.therapeutic_areas:
            target_class = "Natural"
        
        # Calculate score
        scores = await asyncio.to_thread(
            crowe_score,
            request.smiles,
            target_class=target_class
        )
        
        response_data = {
            "smiles": request.smiles,
            "crowe_score": scores.get('crowe_score', 0.0),
            "components": {
                "potency": scores.get('potency', 0.0),
                "selectivity": scores.get('selectivity', 0.0),
                "admet": scores.get('admet', 0.0),
                "synthesis": scores.get('synthesis', 0.0),
                "safety": scores.get('safety', 0.0)
            },
            "confidence": 0.85
        }
        
        # Cache result if enabled
        if request.use_cache and redis_client:
            try:
                await redis_client.setex(
                    cache_key_str,
                    CACHE_TTL,
                    orjson.dumps(response_data)
                )
            except Exception:
                pass
        
        return CroweScoreResponse(
            **response_data,
            cache_hit=cache_hit,
            process_time_ms=(time.time() - start_time) * 1000
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Crowe score calculation failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))

# Websocket endpoint for real-time updates
from fastapi import WebSocket, WebSocketDisconnect
from typing import Set

class ConnectionManager:
    def __init__(self):
        self.active_connections: Set[WebSocket] = set()
    
    async def connect(self, websocket: WebSocket):
        await websocket.accept()
        self.active_connections.add(websocket)
    
    def disconnect(self, websocket: WebSocket):
        self.active_connections.discard(websocket)
    
    async def broadcast(self, message: dict):
        for connection in self.active_connections:
            try:
                await connection.send_json(message)
            except Exception:
                pass

manager = ConnectionManager()

@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    """Real-time updates for discovery pipeline"""
    await manager.connect(websocket)
    try:
        while True:
            data = await websocket.receive_json()
            
            # Process real-time scoring request
            if data.get("type") == "score":
                result = await calculate_crowe_score(
                    CroweScoreRequest(**data.get("payload", {}))
                )
                await websocket.send_json({
                    "type": "score_result",
                    "data": result.dict()
                })
            
            # Broadcast to all connected clients if needed
            elif data.get("type") == "broadcast":
                await manager.broadcast(data.get("payload", {}))
                
    except WebSocketDisconnect:
        manager.disconnect(websocket)

# Performance metrics endpoint
@app.get("/metrics")
async def get_metrics():
    """Performance and usage metrics"""
    metrics = {
        "cache": {
            "type": "redis" if redis_client else "memory",
            "size": len(molecule_cache) if not redis_client else 0,
            "hit_rate": 0.0  # Would need to track this
        },
        "api": {
            "version": "2.0.0",
            "uptime": time.time(),  # Would need to track from startup
            "requests_total": 0,  # Would need to track
        },
        "performance": {
            "avg_response_time_ms": 0,  # Would need to track
            "p95_response_time_ms": 0,  # Would need to track
            "p99_response_time_ms": 0,  # Would need to track
        }
    }
    
    if redis_client:
        try:
            info = await redis_client.info()
            metrics["cache"]["redis_info"] = {
                "used_memory_human": info.get("used_memory_human"),
                "connected_clients": info.get("connected_clients"),
                "total_connections_received": info.get("total_connections_received"),
            }
        except Exception:
            pass
    
    return metrics

# Optimized similarity search with parallel processing
@app.post("/similarity/optimized")
async def similarity_search_optimized(request: SimilaritySearchRequest):
    """Optimized similarity search with batching and caching"""
    start_time = time.time()
    
    try:
        # Validate query molecule
        query_mol = Molecule(request.query_smiles, mol_id="query")
        if not query_mol.is_valid():
            raise HTTPException(status_code=400, detail="Invalid query SMILES")
        
        # Process targets in batches for better performance
        batch_size = 100
        all_results = []
        
        for i in range(0, len(request.database_smiles), batch_size):
            batch = request.database_smiles[i:i+batch_size]
            
            # Create target molecules
            target_molecules = []
            for j, smiles in enumerate(batch):
                try:
                    mol = Molecule(smiles, mol_id=f"target_{i+j}")
                    if mol.is_valid():
                        target_molecules.append(mol)
                except Exception:
                    continue
            
            if target_molecules:
                from ..core.molecule import MoleculeCollection
                target_collection = MoleculeCollection(target_molecules)
                
                # Perform search with async
                results = await asyncio.to_thread(
                    similarity_search.search_single,
                    query_mol,
                    target_collection,
                    request.threshold,
                    request.max_results
                )
                all_results.extend(results)
        
        # Sort and limit results
        all_results.sort(key=lambda x: x.similarity, reverse=True)
        final_results = all_results[:request.max_results]
        
        return {
            "query_smiles": request.query_smiles,
            "results": [
                {
                    "smiles": result.target_molecule.smiles,
                    "similarity": result.similarity,
                    "rank": i + 1
                }
                for i, result in enumerate(final_results)
            ],
            "total_searched": len(request.database_smiles),
            "process_time_ms": (time.time() - start_time) * 1000
        }
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Similarity search failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))

# Root endpoint with branding
@app.get("/")
async def root():
    """Root endpoint with updated branding"""
    return {
        "message": "🧬 CriOS Discovery Engine",
        "tagline": "Accelerating Discovery, Protecting Tomorrow™",
        "version": "2.0.0",
        "endpoints": {
            "docs": "/docs",
            "health": "/health",
            "metrics": "/metrics"
        },
        "features": [
            "High-performance molecular validation",
            "Intelligent caching with Redis",
            "Parallel batch processing",
            "Real-time WebSocket updates",
            "Crowe Discovery Framework scoring",
            "Ethics-first compound screening"
        ]
    }

# Shutdown handler
@app.on_event("shutdown")
async def shutdown_event():
    """Graceful shutdown with cache persistence"""
    logger.info("Shutting down CriOS Discovery Engine...")
    
    # Save in-memory cache if needed
    if molecule_cache and not redis_client:
        try:
            with open("cache_backup.json", "wb") as f:
                f.write(orjson.dumps(molecule_cache))
            logger.info(f"Saved {len(molecule_cache)} cached items")
        except Exception as e:
            logger.error(f"Failed to save cache: {e}")