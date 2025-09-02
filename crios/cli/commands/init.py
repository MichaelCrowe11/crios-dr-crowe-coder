"""
CriOS CLI - Init Command
Project initialization and scaffolding
"""

import logging
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn

logger = logging.getLogger(__name__)
console = Console()


def init_project(
    directory: Optional[Path] = typer.Argument(None, help="Project directory (default: current)"),
    config_only: bool = typer.Option(False, "--config-only", help="Only create config files"),
    overwrite: bool = typer.Option(False, "--overwrite", help="Overwrite existing files")
):
    """
    🚀 Initialize CriOS project structure
    
    Creates configuration files, example data, and project scaffolding.
    
    Examples:
        crios init                    # Initialize in current directory
        crios init my_project         # Initialize in new directory
        crios init --config-only      # Only create config files
    """
    
    # Determine project directory
    if directory is None:
        project_dir = Path.cwd()
    else:
        project_dir = Path(directory)
        project_dir.mkdir(parents=True, exist_ok=True)
    
    console.print(f"🚀 Initializing CriOS project in [cyan]{project_dir}[/cyan]")
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:
        
        # Create directory structure
        task1 = progress.add_task("Creating directory structure...", total=None)
        _create_directory_structure(project_dir, config_only)
        progress.advance(task1, advance=1)
        
        # Create configuration files
        task2 = progress.add_task("Creating configuration files...", total=None)
        _create_config_files(project_dir, overwrite)
        progress.advance(task2, advance=1)
        
        if not config_only:
            # Create example data
            task3 = progress.add_task("Creating example data...", total=None)
            _create_example_data(project_dir, overwrite)
            progress.advance(task3, advance=1)
            
            # Create documentation
            task4 = progress.add_task("Creating documentation...", total=None)
            _create_documentation(project_dir, overwrite)
            progress.advance(task4, advance=1)
    
    console.print("✅ Project initialized successfully!")
    console.print("\n📋 Next steps:")
    console.print("  1. Review and customize [cyan]config/config.yaml[/cyan]")
    console.print("  2. Review safety policies in [cyan]config/ethics.yaml[/cyan]")
    console.print("  3. Try the demo: [green]crios featurize --in data/examples/molecules.smi[/green]")
    console.print("  4. Read the documentation: [blue]docs/README.md[/blue]")


def _create_directory_structure(project_dir: Path, config_only: bool):
    """Create project directory structure"""
    
    # Always create these directories
    directories = [
        "configs",
        "manifests",
        "logs"
    ]
    
    if not config_only:
        directories.extend([
            "data/inputs",
            "data/outputs", 
            "data/examples",
            "docs",
            "notebooks",
            "scripts"
        ])
    
    for directory in directories:
        dir_path = project_dir / directory
        dir_path.mkdir(parents=True, exist_ok=True)
        
        # Create .gitkeep files in empty directories
        gitkeep = dir_path / ".gitkeep"
        if not gitkeep.exists():
            gitkeep.touch()


def _create_config_files(project_dir: Path, overwrite: bool):
    """Create configuration files"""
    
    # Main configuration
    config_file = project_dir / "configs" / "config.yaml"
    if overwrite or not config_file.exists():
        config_content = """# CriOS Discovery Engine Configuration
# Customize these settings for your project

system:
  name: "CriOS Discovery Engine"
  version: "1.0.0"
  environment: "development"
  seed: 42
  debug: true
  log_level: "INFO"

performance:
  num_workers: 4
  batch_size: 1000
  vectorize_ops: true

drug_like_filters:
  molecular_weight:
    min: 150.0
    max: 500.0
  logp:
    min: -2.0
    max: 5.0
  hbd:
    max: 5
  hba:
    max: 10
  tpsa:
    min: 20.0
    max: 140.0

crowe_weights:
  potency: 0.25
  selectivity: 0.20
  admet: 0.25
  synthesis: 0.20
  novelty: 0.10

target_classes:
  GPCR:
    mw_preference: [300, 450]
    logp_preference: [2.0, 4.0]
    tpsa_preference: [40, 90]
"""
        config_file.write_text(config_content)
    
    # Ethics policy
    ethics_file = project_dir / "configs" / "ethics.yaml"
    if overwrite or not ethics_file.exists():
        ethics_content = """# CriOS Ethics & Safety Policy
# IMPORTANT: Customize this for your research context

metadata:
  version: "1.0.0"
  policy_name: "Research Safety Framework"
  description: "Ethics and safety constraints for compound discovery"

principles:
  - "Science before status"
  - "Discovery before profit"
  - "Research only - no medical advice"
  - "Harm reduction and safety first"

enforcement:
  mode: "strict"  # strict, warning, permissive
  log_all_checks: true
  auto_reject_violations: true

# WARNING: Populate these lists based on your jurisdiction and research needs
structural_alerts:
  examples:
    - name: "epoxide"
      smarts: "C1OC1"
      severity: "high"
      description: "Epoxide group - potential DNA reactivity"

controlled_substances:
  # IMPORTANT: These lists are intentionally empty
  # Users must populate based on their research context
  blocked_hashes: []
  blocked_patterns: []
"""
        ethics_file.write_text(ethics_content)


def _create_example_data(project_dir: Path, overwrite: bool):
    """Create example data files"""
    
    # Example molecules SMILES file
    molecules_file = project_dir / "data" / "examples" / "molecules.smi"
    if overwrite or not molecules_file.exists():
        molecules_content = """# Example molecules for CriOS demonstration
# SMILES format: one molecule per line
CCO	ethanol
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O	ibuprofen
CN1C=NC2=C1C(=O)N(C(=O)N2C)C	caffeine
CC1=CC2=C(C=C1C)N(C=N2)C3C(C(C(O3)CO)O)O	riboflavin
C1=CC=C(C=C1)C(=O)O	benzoic_acid
CC(C)(C)C1=CC=C(C=C1)O	BHT
C1CCC(CC1)N	cyclohexylamine
C1=CC=C2C(=C1)C=CC=C2	naphthalene
CC1=CC=CC=C1C	toluene
CCCCCCCCCCCCCCC(=O)O	palmitic_acid
"""
        molecules_file.write_text(molecules_content)
    
    # Example CSV file
    csv_file = project_dir / "data" / "examples" / "compounds.csv"
    if overwrite or not csv_file.exists():
        csv_content = """id,smiles,name,activity
mol_001,CCO,ethanol,inactive
mol_002,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O,ibuprofen,active
mol_003,CN1C=NC2=C1C(=O)N(C(=O)N2C)C,caffeine,active
mol_004,CC1=CC2=C(C=C1C)N(C=N2)C3C(C(C(O3)CO)O)O,riboflavin,active
mol_005,C1=CC=C(C=C1)C(=O)O,benzoic_acid,inactive
"""
        csv_file.write_text(csv_content)


def _create_documentation(project_dir: Path, overwrite: bool):
    """Create project documentation"""
    
    readme_file = project_dir / "docs" / "README.md"
    if overwrite or not readme_file.exists():
        readme_content = """# CriOS Project Documentation

## Getting Started

This CriOS project contains configuration files, example data, and documentation
to help you get started with compound discovery and analysis.

## Project Structure

```
├── configs/           # Configuration files
│   ├── config.yaml   # Main system configuration
│   └── ethics.yaml   # Ethics and safety policies
├── data/             # Data files
│   ├── inputs/       # Input data files
│   ├── outputs/      # Generated output files
│   └── examples/     # Example molecules and compounds
├── docs/             # Documentation
├── logs/             # Log files
├── manifests/        # Run provenance files
├── notebooks/        # Jupyter notebooks
└── scripts/          # Custom scripts
```

## Quick Start

1. **Validate the example molecules:**
   ```bash
   crios featurize --in data/examples/molecules.smi --out data/outputs/features.csv
   ```

2. **Score compounds using Crowe methodology:**
   ```bash
   crios score --in data/examples/molecules.smi --engine synthetic --out data/outputs/scored.csv
   ```

3. **Check ethics compliance:**
   ```bash
   crios explain --in data/examples/molecules.smi --policy configs/ethics.yaml
   ```

## Configuration

### Main Configuration (`configs/config.yaml`)

Customize system settings, performance parameters, drug-like filters, and scoring weights.

### Ethics Policy (`configs/ethics.yaml`)

**IMPORTANT**: Review and customize the ethics policy for your research context:

- Structural alerts and toxicophore patterns
- Controlled substance filters (populate based on your jurisdiction)
- Contextual guardrails and safety checks

## Safety and Ethics

CriOS includes comprehensive safety and ethics frameworks:

- **Structural Alerts**: Automated detection of potentially harmful substructures
- **Controlled Substance Screening**: Hash and pattern-based detection (user-configured)
- **Contextual Guardrails**: Protection against misuse scenarios
- **Audit Trails**: Complete provenance tracking for all operations

## Example Workflows

### Natural Product Discovery
```bash
crios score --in natural_products.smi --engine natural --optimize-for potency novelty
```

### Neurotherapeutic Design
```bash
crios design --engine neuro --target-class GPCR --n 100 --top 20
```

### Similarity Analysis
```bash
crios similarity --in query.smi --ref database.smi --metric tanimoto --out similarity.csv
```

## Support

- Documentation: See additional files in this `docs/` directory
- Issues: Report problems via your project management system
- Ethics Questions: Consult with your institutional ethics board

## Disclaimer

**RESEARCH USE ONLY**: This software is for research purposes only and should not be used 
for medical diagnosis, treatment, or drug development without appropriate regulatory oversight.
"""
        readme_file.write_text(readme_content)