"""
CriOS CLI - Main Application
Typer-based command line interface for CriOS Discovery Engine
"""

import logging
import sys
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.logging import RichHandler
from rich.traceback import install as install_rich_traceback

# Install rich traceback for better error display
install_rich_traceback()

# Initialize Typer app
app = typer.Typer(
    name="crios",
    help="🧬 CriOS Discovery Engine - Universal Compound Discovery Platform",
    add_completion=False,
    rich_markup_mode="rich"
)

# Global console for rich output
console = Console()

# Global configuration
global_config = {
    "seed": None,
    "num_workers": None,
    "config_file": None,
    "policy_file": None,
    "log_json": False,
    "debug": False
}


def setup_logging(debug: bool = False, log_json: bool = False):
    """Setup logging configuration"""
    level = logging.DEBUG if debug else logging.INFO
    
    if log_json:
        # JSON logging for production
        logging.basicConfig(
            level=level,
            format='{"timestamp": "%(asctime)s", "level": "%(levelname)s", "message": "%(message)s"}',
            handlers=[logging.StreamHandler()]
        )
    else:
        # Rich logging for development
        logging.basicConfig(
            level=level,
            format="%(message)s",
            datefmt="[%X]",
            handlers=[RichHandler(console=console, rich_tracebacks=True)]
        )


@app.callback()
def main(
    seed: Optional[int] = typer.Option(None, "--seed", help="Random seed for reproducibility"),
    num_workers: Optional[int] = typer.Option(None, "--num-workers", help="Number of parallel workers"),
    config: Optional[Path] = typer.Option(None, "--config", help="Configuration file path"),
    policy: Optional[Path] = typer.Option(None, "--policy", help="Ethics policy file path"),
    log_json: bool = typer.Option(False, "--log-json", help="Use JSON logging format"),
    debug: bool = typer.Option(False, "--debug", help="Enable debug logging"),
):
    """
    🧬 CriOS Discovery Engine
    
    Universal compound discovery platform with ethical AI constraints.
    
    Examples:
        crios featurize --in molecules.smi --out features.parquet
        crios score --in molecules.smi --engine synthetic --out scored.csv
        crios design --engine synthetic --target-class GPCR --n 100 --top 20
    """
    # Setup logging first
    setup_logging(debug, log_json)
    
    # Store global configuration
    global_config.update({
        "seed": seed,
        "num_workers": num_workers,
        "config_file": config,
        "policy_file": policy,
        "log_json": log_json,
        "debug": debug
    })
    
    # Set random seed if provided
    if seed is not None:
        import random
        import numpy as np
        random.seed(seed)
        np.random.seed(seed)
        
        # Try to set RDKit seed (if available)
        try:
            from rdkit.rdBase import RandomSeed
            RandomSeed(seed)
        except ImportError:
            pass


# Import and register commands
from .commands import (
    init,
    featurize, 
    score,
    design,
    similarity,
    explain
)

# Register commands with the app
app.command("init")(init.init_project)
app.command("featurize")(featurize.featurize_molecules)
app.command("score")(score.score_molecules)
app.command("design")(design.design_molecules)
app.command("similarity")(similarity.calculate_similarity)
app.command("explain")(explain.explain_molecules)


if __name__ == "__main__":
    app()