"""
CriOS CLI - Cluster Commands
Molecular clustering and chemical space analysis
"""

import logging
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from ...core.molecule import MoleculeCollection
from ...core.similarity import ClusterAnalyzer
from ..main import app

logger = logging.getLogger(__name__)
console = Console()


@app.command("cluster")
def cluster_molecules(
    input_file: Path = typer.Option(..., "--input", "-i", help="Input SDF file"),
    algorithm: str = typer.Option("butina", "--algorithm", "-a", help="Clustering algorithm"),
    threshold: float = typer.Option(0.6, "--threshold", "-t", help="Distance threshold"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output results")
):
    """
    🎯 Cluster molecules by structural similarity
    
    Group compounds into clusters based on structural similarity
    for chemical space analysis and diversity assessment.
    """
    try:
        console.print(f"📂 Loading molecules from {input_file}...")
        collection = MoleculeCollection.from_sdf(input_file)
        console.print(f"✅ Loaded {len(collection)} molecules")
        
        console.print("🎯 Clustering molecules...")
        analyzer = ClusterAnalyzer(algorithm=algorithm)
        
        if algorithm == "butina":
            result = analyzer.cluster_butina(collection, distance_threshold=threshold)
        else:
            result = analyzer.cluster_hierarchical(collection, distance_threshold=threshold)
        
        analysis = analyzer.analyze_clusters(collection, result)
        
        console.print(f"📊 Clustering Results:")
        console.print(f"  Clusters found: {analysis['n_clusters']}")
        console.print(f"  Average cluster size: {analysis['average_cluster_size']:.1f}")
        console.print(f"  Clustering efficiency: {analysis['clustering_efficiency']:.1%}")
        
    except Exception as e:
        logger.error(f"Clustering failed: {e}")
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)