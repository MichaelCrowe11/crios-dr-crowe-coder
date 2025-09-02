"""
CriOS CLI - Score Command
Molecular scoring using Crowe methodology and multi-objective optimization
"""

import logging
from pathlib import Path
from typing import List, Optional
import csv
import json

import typer
from rich.console import Console
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn
from rich.table import Table

from ...scoring.design import crowe_score, _passes_drug_like_filters
from ...ethics.filters import check_ethics_compliance
from ...io.schemas import TargetClass, EngineType
from ..main import global_config

logger = logging.getLogger(__name__)
console = Console()


def score_molecules(
    input_file: Path = typer.Option(..., "--in", "-i", help="Input file (SMILES, CSV, or SDF)"),
    output_file: Optional[Path] = typer.Option(None, "--out", "-o", help="Output file path"),
    engine: EngineType = typer.Option(EngineType.SYNTHETIC, "--engine", "-e", help="Scoring engine"),
    target_class: TargetClass = typer.Option(TargetClass.GPCR, "--target-class", "-t", help="Target protein class"),
    optimize_for: List[str] = typer.Option(["potency", "selectivity", "admet"], "--opt", help="Optimization objectives"),
    include_ethics: bool = typer.Option(True, "--ethics/--no-ethics", help="Include ethics screening"),
    format: str = typer.Option("csv", "--format", "-f", help="Output format (csv, json, sdf)"),
    top_n: Optional[int] = typer.Option(None, "--top", help="Only output top N compounds"),
    min_score: Optional[float] = typer.Option(None, "--min-score", help="Minimum Crowe score threshold")
):
    """
    📊 Score molecules using Crowe methodology
    
    Calculate comprehensive scores including potency, selectivity, ADMET, 
    synthetic accessibility, and novelty using the Crowe Discovery Framework.
    
    Examples:
        crios score --in molecules.smi --out scored.csv
        crios score --in compounds.csv --engine neuro --target-class kinase
        crios score --in drugs.sdf --opt potency admet --top 50 --min-score 0.5
    """
    
    console.print(f"📊 Scoring molecules from [cyan]{input_file}[/cyan]")
    console.print(f"🎯 Engine: [green]{engine.value}[/green], Target: [green]{target_class.value}[/green]")
    console.print(f"🔧 Optimizing for: [yellow]{', '.join(optimize_for)}[/yellow]")
    
    try:
        # Load input molecules
        molecules = _load_molecules(input_file)
        console.print(f"✅ Loaded {len(molecules)} molecules")
        
        if not molecules:
            console.print("[red]❌ No valid molecules found in input file[/red]")
            raise typer.Exit(1)
        
        # Score molecules with progress bar
        scored_molecules = []
        
        with Progress(
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=console
        ) as progress:
            
            scoring_task = progress.add_task("Scoring molecules...", total=len(molecules))
            
            for mol_data in molecules:
                smiles = mol_data["smiles"]
                mol_id = mol_data.get("id", smiles[:20])
                
                try:
                    # Calculate Crowe score
                    scores = crowe_score(
                        smiles=smiles,
                        target_class=target_class.value,
                        config=_get_scoring_config()
                    )
                    
                    # Ethics screening
                    ethics_result = None
                    if include_ethics:
                        policy_file = global_config.get("policy_file")
                        ethics_result = check_ethics_compliance(smiles, policy_file)
                    
                    # Compile result
                    result = {
                        "id": mol_id,
                        "smiles": smiles,
                        "target_class": target_class.value,
                        "engine": engine.value,
                        **scores,
                        "ethics_passed": ethics_result["passed"] if ethics_result else True,
                        "ethics_violations": len(ethics_result["violations"]) if ethics_result else 0
                    }
                    
                    # Add original properties if available
                    for key, value in mol_data.items():
                        if key not in ["smiles", "id"]:
                            result[f"original_{key}"] = value
                    
                    scored_molecules.append(result)
                    
                except Exception as e:
                    logger.warning(f"Failed to score {smiles}: {e}")
                    continue
                
                progress.advance(scoring_task)
        
        console.print(f"✅ Scored {len(scored_molecules)} molecules successfully")
        
        # Filter results
        if min_score is not None:
            before_count = len(scored_molecules)
            scored_molecules = [m for m in scored_molecules if m.get("crowe_score", 0) >= min_score]
            console.print(f"🔍 Filtered to {len(scored_molecules)} molecules (score >= {min_score})")
        
        if include_ethics:
            before_count = len(scored_molecules)
            scored_molecules = [m for m in scored_molecules if m.get("ethics_passed", False)]
            console.print(f"⚖️ Ethics screening passed: {len(scored_molecules)}/{before_count} molecules")
        
        # Sort by Crowe score
        scored_molecules.sort(key=lambda x: x.get("crowe_score", 0), reverse=True)
        
        # Limit results
        if top_n is not None:
            scored_molecules = scored_molecules[:top_n]
            console.print(f"📋 Limited to top {len(scored_molecules)} compounds")
        
        # Display summary
        _display_scoring_summary(scored_molecules)
        
        # Save results
        if output_file:
            _save_results(scored_molecules, output_file, format)
            console.print(f"💾 Results saved to [cyan]{output_file}[/cyan]")
        else:
            _display_top_results(scored_molecules[:10])
        
        console.print("🎉 Scoring complete!")
        
    except Exception as e:
        logger.error(f"Scoring failed: {e}")
        console.print(f"[red]❌ Scoring failed: {e}[/red]")
        raise typer.Exit(1)


def _load_molecules(input_file: Path) -> List[dict]:
    """Load molecules from various input formats"""
    
    molecules = []
    
    try:
        if input_file.suffix.lower() == '.csv':
            # CSV format
            with open(input_file, 'r') as f:
                reader = csv.DictReader(f)
                for i, row in enumerate(reader):
                    smiles = row.get('smiles', '').strip()
                    if smiles:
                        mol_data = {
                            "id": row.get('id', f'mol_{i+1}'),
                            "smiles": smiles
                        }
                        # Add other columns as properties
                        for key, value in row.items():
                            if key not in ['id', 'smiles']:
                                mol_data[key] = value
                        molecules.append(mol_data)
        
        elif input_file.suffix.lower() in ['.smi', '.smiles']:
            # SMILES format
            with open(input_file, 'r') as f:
                for i, line in enumerate(f):
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split('\t')
                        smiles = parts[0].strip()
                        name = parts[1].strip() if len(parts) > 1 else f'mol_{i+1}'
                        
                        molecules.append({
                            "id": name,
                            "smiles": smiles
                        })
        
        elif input_file.suffix.lower() == '.sdf':
            # SDF format (requires RDKit)
            try:
                from rdkit import Chem
                supplier = Chem.SDMolSupplier(str(input_file))
                
                for i, mol in enumerate(supplier):
                    if mol is not None:
                        smiles = Chem.MolToSmiles(mol)
                        mol_id = mol.GetProp('_Name') if mol.HasProp('_Name') else f'mol_{i+1}'
                        
                        mol_data = {
                            "id": mol_id,
                            "smiles": smiles
                        }
                        
                        # Extract properties
                        for prop_name in mol.GetPropNames():
                            mol_data[prop_name] = mol.GetProp(prop_name)
                        
                        molecules.append(mol_data)
            except ImportError:
                console.print("[red]❌ RDKit required for SDF format[/red]")
                raise typer.Exit(1)
        
        else:
            console.print(f"[red]❌ Unsupported file format: {input_file.suffix}[/red]")
            raise typer.Exit(1)
    
    except FileNotFoundError:
        console.print(f"[red]❌ Input file not found: {input_file}[/red]")
        raise typer.Exit(1)
    except Exception as e:
        console.print(f"[red]❌ Failed to load molecules: {e}[/red]")
        raise typer.Exit(1)
    
    return molecules


def _get_scoring_config():
    """Get scoring configuration"""
    config = {}
    
    # Load from config file if specified
    config_file = global_config.get("config_file")
    if config_file and config_file.exists():
        try:
            import yaml
            with open(config_file, 'r') as f:
                full_config = yaml.safe_load(f)
                config = full_config.get("drug_like_filters", {})
        except Exception as e:
            logger.warning(f"Failed to load config: {e}")
    
    return config


def _display_scoring_summary(molecules: List[dict]):
    """Display scoring summary statistics"""
    
    if not molecules:
        return
    
    # Calculate statistics
    crowe_scores = [m.get("crowe_score", 0) for m in molecules]
    drug_like_count = sum(1 for m in molecules if m.get("drug_like", False))
    ethics_passed_count = sum(1 for m in molecules if m.get("ethics_passed", True))
    
    # Component score averages
    components = ["potency", "selectivity", "admet", "synthesis", "novelty"]
    avg_scores = {}
    for comp in components:
        scores = [m.get(comp, 0) for m in molecules if comp in m]
        avg_scores[comp] = sum(scores) / len(scores) if scores else 0
    
    console.print("\n📊 Scoring Summary:")
    console.print(f"  Total molecules: {len(molecules)}")
    console.print(f"  Average Crowe score: {sum(crowe_scores) / len(crowe_scores):.3f}")
    console.print(f"  Drug-like compounds: {drug_like_count} ({drug_like_count/len(molecules)*100:.1f}%)")
    console.print(f"  Ethics compliant: {ethics_passed_count} ({ethics_passed_count/len(molecules)*100:.1f}%)")
    
    console.print("\n🎯 Average Component Scores:")
    for comp, score in avg_scores.items():
        console.print(f"  {comp.capitalize()}: {score:.3f}")


def _display_top_results(molecules: List[dict]):
    """Display top scoring molecules in a table"""
    
    if not molecules:
        return
    
    table = Table(title="Top Scoring Molecules")
    table.add_column("Rank", style="cyan", width=4)
    table.add_column("ID", style="white", width=12)
    table.add_column("SMILES", style="yellow", width=30)
    table.add_column("Crowe Score", style="green", width=10)
    table.add_column("Drug-like", style="blue", width=8)
    table.add_column("Ethics", style="red", width=6)
    
    for i, mol in enumerate(molecules, 1):
        table.add_row(
            str(i),
            mol.get("id", "")[:12],
            mol.get("smiles", "")[:30],
            f"{mol.get('crowe_score', 0):.3f}",
            "✓" if mol.get("drug_like", False) else "✗",
            "✓" if mol.get("ethics_passed", True) else "✗"
        )
    
    console.print("\n")
    console.print(table)


def _save_results(molecules: List[dict], output_file: Path, format: str):
    """Save results to file"""
    
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    if format.lower() == "csv":
        # CSV format
        if molecules:
            fieldnames = list(molecules[0].keys())
            
            with open(output_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(molecules)
    
    elif format.lower() == "json":
        # JSON format
        with open(output_file, 'w') as f:
            json.dump(molecules, f, indent=2, default=str)
    
    elif format.lower() == "sdf":
        # SDF format (requires RDKit)
        try:
            from rdkit import Chem
            
            with Chem.SDWriter(str(output_file)) as writer:
                for mol_data in molecules:
                    smiles = mol_data.get("smiles")
                    if not smiles:
                        continue
                    
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is None:
                        continue
                    
                    # Add properties
                    for key, value in mol_data.items():
                        if key != "smiles" and value is not None:
                            mol.SetProp(key, str(value))
                    
                    writer.write(mol)
        except ImportError:
            console.print("[red]❌ RDKit required for SDF output format[/red]")
            raise typer.Exit(1)
    
    else:
        console.print(f"[red]❌ Unsupported output format: {format}[/red]")
        raise typer.Exit(1)