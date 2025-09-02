"""
CriOS CLI - Search Commands
Molecular similarity search and compound retrieval
"""

import logging
from pathlib import Path
from typing import Optional, List

import typer
from rich.console import Console
from rich.table import Table

from ...core.molecule import Molecule, MoleculeCollection
from ...core.similarity import SimilaritySearch
from ..main import app

logger = logging.getLogger(__name__)
console = Console()


@app.command("search")
def similarity_search(
    query: str = typer.Option(..., "--query", "-q", help="Query SMILES string"),
    database: Path = typer.Option(..., "--database", "-d", help="Database file (SDF)"),
    threshold: float = typer.Option(0.7, "--threshold", "-t", help="Similarity threshold"),
    max_results: int = typer.Option(10, "--max-results", "-n", help="Maximum results"),
    fingerprint: str = typer.Option("morgan", "--fingerprint", "-f", help="Fingerprint type"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output file")
):
    """
    🔍 Perform molecular similarity search
    
    Search for structurally similar compounds in a database using
    molecular fingerprints and similarity metrics.
    """
    try:
        # Create query molecule
        console.print(f"🎯 Query: {query}")
        query_mol = Molecule(query, mol_id="query")
        
        # Load database
        console.print(f"📚 Loading database: {database}")
        db_collection = MoleculeCollection.from_sdf(database)
        console.print(f"✅ Loaded {len(db_collection)} compounds")
        
        # Perform search
        console.print("🔍 Searching for similar compounds...")
        search_engine = SimilaritySearch(fingerprint_type=fingerprint)
        
        results = search_engine.search_single(
            query_mol, db_collection, threshold, max_results
        )
        
        # Display results
        console.print(f"\n📊 Found {len(results)} similar compounds:")
        
        table = Table()
        table.add_column("Rank", style="cyan")
        table.add_column("ID", style="white") 
        table.add_column("Similarity", style="green")
        table.add_column("SMILES", style="yellow")
        
        for i, result in enumerate(results, 1):
            table.add_row(
                str(i),
                result.target_id,
                f"{result.similarity:.3f}",
                result.target_molecule.smiles[:50] + "..." if len(result.target_molecule.smiles) > 50 else result.target_molecule.smiles
            )
        
        console.print(table)
        
        if output:
            # Save results
            import csv
            with open(output, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['rank', 'id', 'similarity', 'smiles'])
                for i, result in enumerate(results, 1):
                    writer.writerow([i, result.target_id, result.similarity, result.target_molecule.smiles])
            console.print(f"✅ Results saved to {output}")
        
    except Exception as e:
        logger.error(f"Search failed: {e}")
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)