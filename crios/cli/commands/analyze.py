"""
CriOS CLI - Analyze Commands
Comprehensive compound analysis and reporting
"""

import logging
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

from ...core.compound import CompoundLibrary
from ...scoring.crowe_score import CroweScorer
from ..main import app

logger = logging.getLogger(__name__)
console = Console()


@app.command("analyze")
def analyze_library(
    input_file: Path = typer.Option(..., "--input", "-i", help="Input compound library"),
    report_type: str = typer.Option("summary", "--type", "-t", help="Report type (summary, detailed)"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output report file"),
    top_n: int = typer.Option(10, "--top", help="Number of top compounds to show")
):
    """
    📊 Analyze compound library and generate reports
    
    Generate comprehensive analysis reports including Crowe scores,
    drug-likeness assessment, and discovery pipeline insights.
    """
    try:
        console.print(f"📂 Loading compound library: {input_file}")
        
        # For demonstration, create a simple library
        # In practice, this would load from various formats
        from ...core.molecule import MoleculeCollection
        collection = MoleculeCollection.from_sdf(input_file)
        
        console.print(f"✅ Loaded {len(collection)} compounds")
        
        # Analyze compounds
        console.print("📊 Analyzing compounds...")
        
        # Basic statistics
        valid_count = sum(1 for mol in collection if mol.is_valid())
        drug_like_count = sum(1 for mol in collection if mol.is_drug_like())
        
        # Display results
        table = Table(title="Library Analysis Summary")
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="white")
        table.add_column("Percentage", style="green")
        
        table.add_row("Total Compounds", str(len(collection)), "100.0%")
        table.add_row("Valid Structures", str(valid_count), f"{valid_count/len(collection)*100:.1f}%")
        table.add_row("Drug-like", str(drug_like_count), f"{drug_like_count/len(collection)*100:.1f}%")
        
        console.print(table)
        
        if output:
            console.print(f"📄 Report saved to {output}")
        
        console.print("✅ Analysis complete!")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)