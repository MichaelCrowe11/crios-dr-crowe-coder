"""
CriOS Dr. Crowe Coder Platform - Backend API
FastAPI backend for 194 PhD Agent Drug Discovery System
"""

from fastapi import FastAPI, HTTPException, WebSocket, WebSocketDisconnect, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, StreamingResponse
from pydantic import BaseModel, Field
from typing import List, Dict, Optional, Any
from datetime import datetime
import asyncio
import json
import uuid
from enum import Enum
import sys
import os

# Add CriOS to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from src.agents.dr_crowe_coder import (
    DrCroweCoderSystem,
    PhDDivision,
    PhDAgent,
    PipelineOrchestrator,
    ClaudeCodeEngine
)
from src.core.molecule import CriOSMolecule
from src.analysis.similarity import tanimoto_search
from src.core.filters import MolecularFilter

# Initialize FastAPI app
app = FastAPI(
    title="CriOS Dr. Crowe Coder Platform",
    description="Revolutionary AI-Driven Drug Discovery with 194 PhD Agents",
    version="1.0.0",
    docs_url="/api/docs",
    redoc_url="/api/redoc"
)

# CORS configuration
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "https://crios.ai"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Initialize Dr. Crowe Coder System
dr_crowe_system = DrCroweCoderSystem()

# WebSocket connection manager
class ConnectionManager:
    def __init__(self):
        self.active_connections: List[WebSocket] = []
        self.agent_statuses: Dict[str, str] = {}
        
    async def connect(self, websocket: WebSocket):
        await websocket.accept()
        self.active_connections.append(websocket)
        
    def disconnect(self, websocket: WebSocket):
        self.active_connections.remove(websocket)
        
    async def broadcast(self, message: dict):
        for connection in self.active_connections:
            try:
                await connection.send_json(message)
            except:
                pass

manager = ConnectionManager()

# ============================================================================
# PYDANTIC MODELS
# ============================================================================

class DiscoveryRequest(BaseModel):
    target: str = Field(..., description="Target disease or protein")
    pipeline: str = Field(default="ai_generative", description="Discovery pipeline")
    num_agents: int = Field(default=24, description="Number of agents to deploy")
    max_compounds: int = Field(default=100, description="Maximum compounds to generate")
    
class MoleculeRequest(BaseModel):
    smiles: str = Field(..., description="SMILES string")
    calculate_descriptors: bool = Field(default=True)
    
class SimilaritySearchRequest(BaseModel):
    query_smiles: str = Field(..., description="Query SMILES")
    database_compounds: List[str] = Field(..., description="Database SMILES")
    threshold: float = Field(default=0.7, description="Similarity threshold")
    
class FilterRequest(BaseModel):
    compounds: List[str] = Field(..., description="SMILES strings")
    rule: str = Field(default="Lipinski", description="Filter rule")
    
class AgentTaskRequest(BaseModel):
    task_description: str = Field(..., description="Task to perform")
    agent_ids: Optional[List[str]] = Field(None, description="Specific agents to use")
    division: Optional[str] = Field(None, description="Division to use")

class PipelineStatus(BaseModel):
    pipeline_id: str
    status: str
    current_stage: int
    total_stages: int
    agents_deployed: List[str]
    compounds_generated: int
    start_time: datetime
    estimated_completion: Optional[datetime]
    
# ============================================================================
# API ENDPOINTS
# ============================================================================

@app.get("/")
async def root():
    """Root endpoint with system info"""
    return {
        "name": "CriOS Dr. Crowe Coder Platform",
        "version": "1.0.0",
        "agents": len(dr_crowe_system.agents),
        "divisions": len(PhDDivision),
        "pipelines": len(dr_crowe_system.pipeline_orchestrator.pipelines),
        "status": "operational",
        "tagline": "Where 194 PhD minds converge to cure disease"
    }

@app.get("/api/agents")
async def get_agents(division: Optional[str] = None):
    """Get all PhD agents or filter by division"""
    agents = dr_crowe_system.agents.values()
    
    if division:
        agents = [a for a in agents if a.division.value == division]
    
    return {
        "total": len(list(agents)),
        "agents": [
            {
                "id": a.agent_id,
                "name": a.name,
                "title": a.title,
                "division": a.division.value,
                "specialization": a.specialization,
                "h_index": a.h_index,
                "publications": a.publications,
                "patents": a.patents
            }
            for a in agents
        ]
    }

@app.get("/api/agents/{agent_id}")
async def get_agent(agent_id: str):
    """Get specific agent details"""
    agent = dr_crowe_system.agents.get(agent_id)
    if not agent:
        raise HTTPException(status_code=404, detail="Agent not found")
    
    return {
        "id": agent.agent_id,
        "name": agent.name,
        "title": agent.title,
        "division": agent.division.value,
        "specialization": agent.specialization,
        "h_index": agent.h_index,
        "publications": agent.publications,
        "patents": agent.patents,
        "breakthrough_discoveries": agent.breakthrough_discoveries,
        "collaboration_style": agent.collaboration_style,
        "thinking_pattern": agent.thinking_pattern
    }

@app.get("/api/pipelines")
async def get_pipelines():
    """Get all discovery pipelines"""
    pipelines = dr_crowe_system.pipeline_orchestrator.pipelines
    
    return {
        "total": len(pipelines),
        "pipelines": [
            {
                "id": name,
                "name": config["name"],
                "stages": len(config["stages"]),
                "stage_details": [
                    {
                        "task": stage["task"],
                        "agents": stage["agents"]
                    }
                    for stage in config["stages"]
                ]
            }
            for name, config in pipelines.items()
        ]
    }

@app.post("/api/discover")
async def discover_compounds(request: DiscoveryRequest, background_tasks: BackgroundTasks):
    """Launch drug discovery pipeline"""
    pipeline_id = str(uuid.uuid4())
    
    # Start discovery in background
    background_tasks.add_task(
        run_discovery_pipeline,
        pipeline_id,
        request.target,
        request.pipeline,
        request.num_agents,
        request.max_compounds
    )
    
    return {
        "pipeline_id": pipeline_id,
        "status": "started",
        "message": f"Discovery pipeline launched for {request.target}",
        "websocket_url": f"/ws/pipeline/{pipeline_id}"
    }

async def run_discovery_pipeline(
    pipeline_id: str,
    target: str,
    pipeline: str,
    num_agents: int,
    max_compounds: int
):
    """Background task to run discovery pipeline"""
    pipeline_config = dr_crowe_system.pipeline_orchestrator.pipelines.get(pipeline)
    if not pipeline_config:
        return
    
    # Notify via WebSocket
    await manager.broadcast({
        "type": "pipeline_started",
        "pipeline_id": pipeline_id,
        "target": target,
        "pipeline": pipeline,
        "stages": len(pipeline_config["stages"])
    })
    
    results = []
    for i, stage in enumerate(pipeline_config["stages"]):
        # Get agents for this stage
        stage_agents = [dr_crowe_system.agents[aid] for aid in stage["agents"]]
        
        # Notify stage start
        await manager.broadcast({
            "type": "stage_started",
            "pipeline_id": pipeline_id,
            "stage": i + 1,
            "task": stage["task"],
            "agents": [a.name for a in stage_agents]
        })
        
        # Execute stage (simulated)
        await asyncio.sleep(2)  # Simulate processing
        
        stage_result = {
            "stage": i + 1,
            "task": stage["task"],
            "compounds_generated": 10 * (i + 1),
            "success": True
        }
        results.append(stage_result)
        
        # Notify stage complete
        await manager.broadcast({
            "type": "stage_completed",
            "pipeline_id": pipeline_id,
            "stage": i + 1,
            "result": stage_result
        })
    
    # Pipeline complete
    await manager.broadcast({
        "type": "pipeline_completed",
        "pipeline_id": pipeline_id,
        "results": results,
        "total_compounds": sum(r["compounds_generated"] for r in results)
    })

@app.post("/api/molecule/analyze")
async def analyze_molecule(request: MoleculeRequest):
    """Analyze a single molecule"""
    try:
        mol = CriOSMolecule(request.smiles)
        
        result = {
            "smiles": mol.canonical_smiles,
            "is_valid": mol.is_valid,
            "mol_id": mol.mol_id
        }
        
        if request.calculate_descriptors and mol.is_valid:
            from src.core.descriptors import DescriptorCalculator
            calc = DescriptorCalculator(["MW", "LogP", "TPSA", "HBA", "HBD", "RotatableBonds"])
            descriptors = calc.calculate(mol.mol)
            result["descriptors"] = descriptors
            
        return result
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/api/similarity/search")
async def similarity_search(request: SimilaritySearchRequest):
    """Perform similarity search"""
    try:
        # Create candidates list
        candidates = [(f"MOL_{i:06d}", smi) for i, smi in enumerate(request.database_compounds)]
        
        # Perform search
        hits = tanimoto_search(
            request.query_smiles,
            candidates,
            threshold=request.threshold
        )
        
        return {
            "query": request.query_smiles,
            "total_hits": len(hits),
            "threshold": request.threshold,
            "hits": [
                {"mol_id": mol_id, "similarity": sim}
                for mol_id, sim in hits[:100]  # Limit to top 100
            ]
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/api/filter")
async def filter_compounds(request: FilterRequest):
    """Filter compounds by rule"""
    try:
        filter_obj = MolecularFilter(request.rule)
        passed = []
        failed = []
        
        for smiles in request.compounds:
            mol = CriOSMolecule(smiles)
            if filter_obj.evaluate(mol):
                passed.append(smiles)
            else:
                failed.append(smiles)
        
        return {
            "rule": request.rule,
            "total": len(request.compounds),
            "passed": len(passed),
            "failed": len(failed),
            "pass_rate": len(passed) / len(request.compounds) if request.compounds else 0,
            "passed_compounds": passed[:100],  # Limit response size
            "failed_compounds": failed[:100]
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/api/orchestrate")
async def orchestrate_task(request: AgentTaskRequest):
    """Orchestrate custom task with specific agents"""
    # Select agents
    if request.agent_ids:
        agents = [dr_crowe_system.agents[aid] for aid in request.agent_ids 
                 if aid in dr_crowe_system.agents]
    elif request.division:
        agents = [a for a in dr_crowe_system.agents.values() 
                 if a.division.value == request.division][:10]
    else:
        agents = list(dr_crowe_system.agents.values())[:10]
    
    # Execute task (simulated)
    result = await dr_crowe_system.claude_code_engine.execute_task(
        {"description": request.task_description, "objective": "solve"},
        agents
    )
    
    return {
        "task": request.task_description,
        "agents_deployed": [a.name for a in agents],
        "result": result
    }

@app.get("/api/stats")
async def get_system_stats():
    """Get system statistics"""
    agents = list(dr_crowe_system.agents.values())
    
    # Calculate division stats
    division_stats = {}
    for division in PhDDivision:
        division_agents = [a for a in agents if a.division == division]
        if division_agents:
            division_stats[division.value] = {
                "count": len(division_agents),
                "avg_h_index": sum(a.h_index for a in division_agents) / len(division_agents),
                "total_publications": sum(a.publications for a in division_agents),
                "total_patents": sum(a.patents for a in division_agents)
            }
    
    return {
        "total_agents": len(agents),
        "total_divisions": len(PhDDivision),
        "total_pipelines": len(dr_crowe_system.pipeline_orchestrator.pipelines),
        "total_publications": sum(a.publications for a in agents),
        "total_patents": sum(a.patents for a in agents),
        "avg_h_index": sum(a.h_index for a in agents) / len(agents),
        "division_stats": division_stats,
        "top_agents": [
            {"name": a.name, "h_index": a.h_index}
            for a in sorted(agents, key=lambda x: x.h_index, reverse=True)[:10]
        ]
    }

# ============================================================================
# WEBSOCKET ENDPOINTS
# ============================================================================

@app.websocket("/ws/pipeline/{pipeline_id}")
async def websocket_pipeline(websocket: WebSocket, pipeline_id: str):
    """WebSocket for real-time pipeline updates"""
    await manager.connect(websocket)
    try:
        await websocket.send_json({
            "type": "connected",
            "message": f"Connected to pipeline {pipeline_id}"
        })
        
        while True:
            # Keep connection alive
            await asyncio.sleep(1)
            
    except WebSocketDisconnect:
        manager.disconnect(websocket)

@app.websocket("/ws/agents")
async def websocket_agents(websocket: WebSocket):
    """WebSocket for real-time agent status updates"""
    await manager.connect(websocket)
    try:
        # Send initial agent statuses
        await websocket.send_json({
            "type": "agent_statuses",
            "agents": {
                agent_id: "idle" 
                for agent_id in dr_crowe_system.agents.keys()
            }
        })
        
        while True:
            # Simulate agent activity
            await asyncio.sleep(5)
            
            # Random agent becomes active
            import random
            agent_id = random.choice(list(dr_crowe_system.agents.keys()))
            
            await websocket.send_json({
                "type": "agent_active",
                "agent_id": agent_id,
                "status": "working",
                "task": "Processing molecular structure"
            })
            
    except WebSocketDisconnect:
        manager.disconnect(websocket)

# ============================================================================
# HEALTH CHECK
# ============================================================================

@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "timestamp": datetime.now().isoformat(),
        "agents_operational": len(dr_crowe_system.agents),
        "pipelines_available": len(dr_crowe_system.pipeline_orchestrator.pipelines)
    }

# ============================================================================
# STARTUP/SHUTDOWN EVENTS
# ============================================================================

@app.on_event("startup")
async def startup_event():
    """Initialize system on startup"""
    print("=" * 60)
    print("CriOS Dr. Crowe Coder Platform Starting...")
    print(f"Loaded {len(dr_crowe_system.agents)} PhD agents")
    print(f"Available pipelines: {len(dr_crowe_system.pipeline_orchestrator.pipelines)}")
    print("API Documentation: http://localhost:8000/api/docs")
    print("=" * 60)

@app.on_event("shutdown")
async def shutdown_event():
    """Cleanup on shutdown"""
    print("CriOS Platform shutting down...")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        "app:app",
        host="0.0.0.0",
        port=8000,
        reload=True,
        log_level="info"
    )