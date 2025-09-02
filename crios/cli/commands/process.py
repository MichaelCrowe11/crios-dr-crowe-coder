"""
CriOS CLI - Process Commands  
Molecular descriptor calculation and property analysis
"""

import logging
from pathlib import Path
from typing import Optional, List

import typer
from rich.console import Console
from rich.progress import Progress

from ...core.molecule import MoleculeCollection
from ...core.compound import CompoundLibrary
from ...core.descriptors import DescriptorCalculator
from ...scoring.crowe_score import CroweScorer
from ..main import app

logger = logging.getLogger(__name__)
console = Console()


@app.command("process")
def process_molecules(
    input_file: Path = typer.Option(..., "--input", "-i", help="Input file (SDF, CSV)"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output file"),
    descriptors: bool = typer.Option(True, "--descriptors", help="Calculate molecular descriptors"),
    crowe_score: bool = typer.Option(False, "--crowe-score", help="Calculate Crowe scores"),
    drug_filters: bool = typer.Option(False, "--drug-filters", help="Apply drug-like filters"),
    parallel: bool = typer.Option(True, "--parallel", help="Use parallel processing")
):
    """
    🧮 Calculate molecular descriptors and properties
    
    Process compounds to calculate comprehensive molecular descriptors,
    drug-likeness scores, and Crowe Discovery Framework scores.
    """
    try:
        # Load molecules
        console.print(f"📂 Loading molecules from {input_file}...")
        
        if input_file.suffix.lower() == '.sdf':
            collection = MoleculeCollection.from_sdf(input_file)
        else:
            raise typer.BadParameter(f"Unsupported file format: {input_file.suffix}")
        
        console.print(f"✅ Loaded {len(collection)} molecules")
        
        # Process molecules
        with Progress() as progress:
            task = progress.add_task("Processing molecules...", total=len(collection))
            
            if descriptors:
                calc = DescriptorCalculator()
                descriptor_results = calc.calculate_batch(collection, parallel=parallel)
                progress.advance(task, 50)
                
                if output:
                    df = calc.to_dataframe(descriptor_results)
                    df.to_csv(output, index=False)
                    console.print(f"✅ Results saved to {output}")
            
            if crowe_score:
                scorer = CroweScorer()
                # Convert to CompoundLibrary if needed
                if isinstance(collection, MoleculeCollection):
                    compounds = [compound for compound in collection if hasattr(compound, 'biological_activities')]
                    if compounds:
                        compound_lib = CompoundLibrary(compounds)
                        scores = scorer.score_library(compound_lib, parallel=parallel)
                        console.print(f"📊 Calculated Crowe scores for {len(scores)} compounds")
                
                progress.advance(task, 50)
        
        console.print("🎉 Processing complete!")
        
    except Exception as e:
        logger.error(f"Processing failed: {e}")
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)