"""
CriOS CLI - Validate Commands
SMILES validation and molecular structure verification
"""

import logging
from pathlib import Path
from typing import List, Optional, Dict, Any
import csv

import typer
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn

from ...core.molecule import Molecule, validate_smiles, standardize_smiles
from ...exceptions import ValidationError
from ..main import app

logger = logging.getLogger(__name__)
console = Console()


@app.command("validate")
def validate_molecules(
    input_file: Optional[Path] = typer.Option(
        None, "--input", "-i", help="Input file (CSV, SDF, or text file with SMILES)"
    ),
    smiles_column: str = typer.Option(
        "smiles", "--smiles-col", help="Column name for SMILES (CSV input)"
    ),
    id_column: str = typer.Option(
        "id", "--id-col", help="Column name for molecule IDs (CSV input)"
    ),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Output file for validation results"
    ),
    smiles: Optional[List[str]] = typer.Option(
        None, "--smiles", "-s", help="Individual SMILES strings to validate"
    ),
    standardize: bool = typer.Option(
        False, "--standardize", help="Standardize SMILES to canonical form"
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Show detailed validation results"
    )
):
    """
    🧪 Validate molecular structures from SMILES strings
    
    Validates SMILES strings and optionally standardizes them to canonical form.
    Supports input from files (CSV, SDF) or individual SMILES strings.
    
    Examples:
        crios validate -s "CCO" -s "c1ccccc1"
        crios validate -i compounds.csv --standardize
        crios validate -i molecules.sdf -o validation_results.csv
    """
    try:
        validation_results = []
        
        if input_file:
            # Process file input
            validation_results = _validate_from_file(
                input_file, smiles_column, id_column, standardize, verbose
            )
        elif smiles:
            # Process individual SMILES
            validation_results = _validate_smiles_list(
                smiles, standardize, verbose
            )
        else:
            console.print("[red]Error:[/red] No input provided. Use --input or --smiles")
            raise typer.Exit(1)
        
        # Display results
        _display_validation_results(validation_results, verbose)
        
        # Save results if output specified
        if output:
            _save_validation_results(validation_results, output)
            console.print(f"✅ Results saved to {output}")
        
        # Exit with non-zero code if any validation failed
        failed_count = sum(1 for result in validation_results if not result["valid"])
        if failed_count > 0:
            console.print(f"\n[yellow]Warning:[/yellow] {failed_count} molecule(s) failed validation")
            raise typer.Exit(1)
        
        console.print("\n✅ All molecules validated successfully!")
        
    except Exception as e:
        logger.error(f"Validation failed: {e}")
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)


def _validate_from_file(
    input_file: Path,
    smiles_column: str,
    id_column: str,
    standardize: bool,
    verbose: bool
) -> List[Dict[str, Any]]:
    """Validate molecules from file input"""
    
    if not input_file.exists():
        raise ValidationError(f"Input file not found: {input_file}")
    
    results = []
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:
        task = progress.add_task("Validating molecules...", total=None)
        
        if input_file.suffix.lower() == '.csv':
            results = _validate_from_csv(
                input_file, smiles_column, id_column, standardize, verbose
            )
        elif input_file.suffix.lower() == '.sdf':
            results = _validate_from_sdf(input_file, standardize, verbose)
        elif input_file.suffix.lower() in ['.txt', '.smi']:
            results = _validate_from_text(input_file, standardize, verbose)
        else:
            raise ValidationError(f"Unsupported file format: {input_file.suffix}")
        
        progress.update(task, completed=True)
    
    return results


def _validate_from_csv(
    csv_file: Path,
    smiles_column: str,
    id_column: str,
    standardize: bool,
    verbose: bool
) -> List[Dict[str, Any]]:
    """Validate molecules from CSV file"""
    
    results = []
    
    with open(csv_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        
        if smiles_column not in reader.fieldnames:
            raise ValidationError(f"SMILES column '{smiles_column}' not found in CSV")
        
        for i, row in enumerate(reader, 1):
            mol_id = row.get(id_column, f"mol_{i}")
            smiles_str = row.get(smiles_column, "").strip()
            
            if not smiles_str:
                results.append({
                    "id": mol_id,
                    "smiles": smiles_str,
                    "valid": False,
                    "error": "Empty SMILES string",
                    "canonical_smiles": None
                })
                continue
            
            result = _validate_single_smiles(mol_id, smiles_str, standardize)
            results.append(result)
    
    return results


def _validate_from_sdf(
    sdf_file: Path,
    standardize: bool,
    verbose: bool
) -> List[Dict[str, Any]]:
    """Validate molecules from SDF file"""
    
    from rdkit import Chem
    from ...core.molecule import MoleculeCollection
    
    results = []
    collection = MoleculeCollection.from_sdf(sdf_file)
    
    for molecule in collection:
        result = {
            "id": molecule.mol_id,
            "smiles": molecule.smiles,
            "valid": molecule.is_valid(),
            "error": None if molecule.is_valid() else "Invalid molecular structure",
            "canonical_smiles": None
        }
        
        if standardize and molecule.is_valid():
            result["canonical_smiles"] = molecule.smiles
        
        results.append(result)
    
    return results


def _validate_from_text(
    text_file: Path,
    standardize: bool,
    verbose: bool
) -> List[Dict[str, Any]]:
    """Validate molecules from text file (one SMILES per line)"""
    
    results = []
    
    with open(text_file, 'r', encoding='utf-8') as f:
        for i, line in enumerate(f, 1):
            smiles_str = line.strip()
            
            if not smiles_str or smiles_str.startswith('#'):
                continue  # Skip empty lines and comments
            
            mol_id = f"mol_{i}"
            result = _validate_single_smiles(mol_id, smiles_str, standardize)
            results.append(result)
    
    return results


def _validate_smiles_list(
    smiles_list: List[str],
    standardize: bool,
    verbose: bool
) -> List[Dict[str, Any]]:
    """Validate list of SMILES strings"""
    
    results = []
    
    for i, smiles_str in enumerate(smiles_list, 1):
        mol_id = f"mol_{i}"
        result = _validate_single_smiles(mol_id, smiles_str, standardize)
        results.append(result)
    
    return results


def _validate_single_smiles(
    mol_id: str,
    smiles_str: str,
    standardize: bool
) -> Dict[str, Any]:
    """Validate a single SMILES string"""
    
    result = {
        "id": mol_id,
        "smiles": smiles_str,
        "valid": False,
        "error": None,
        "canonical_smiles": None
    }
    
    try:
        # Basic validation
        is_valid = validate_smiles(smiles_str)
        result["valid"] = is_valid
        
        if not is_valid:
            result["error"] = "Invalid SMILES syntax"
        else:
            # Try to create molecule for additional validation
            molecule = Molecule(smiles_str, mol_id=mol_id)
            result["valid"] = molecule.is_valid()
            
            if not molecule.is_valid():
                result["error"] = "Molecule could not be processed by RDKit"
            
            # Standardize if requested
            if standardize and molecule.is_valid():
                canonical = standardize_smiles(smiles_str)
                result["canonical_smiles"] = canonical
                
    except Exception as e:
        result["valid"] = False
        result["error"] = str(e)
    
    return result


def _display_validation_results(
    results: List[Dict[str, Any]],
    verbose: bool
):
    """Display validation results in a table"""
    
    if not results:
        console.print("[yellow]No molecules to validate[/yellow]")
        return
    
    # Summary statistics
    total = len(results)
    valid_count = sum(1 for r in results if r["valid"])
    invalid_count = total - valid_count
    
    console.print(f"\n📊 Validation Summary:")
    console.print(f"  Total molecules: {total}")
    console.print(f"  Valid: [green]{valid_count}[/green]")
    console.print(f"  Invalid: [red]{invalid_count}[/red]")
    console.print(f"  Success rate: {valid_count/total*100:.1f}%")
    
    if verbose or invalid_count > 0:
        # Detailed results table
        table = Table(title="Validation Results")
        table.add_column("ID", style="cyan")
        table.add_column("SMILES", style="white")
        table.add_column("Status", style="bold")
        table.add_column("Error", style="red")
        
        if any(r.get("canonical_smiles") for r in results):
            table.add_column("Canonical SMILES", style="green")
        
        for result in results:
            status = "✅ Valid" if result["valid"] else "❌ Invalid"
            error = result.get("error", "") or ""
            
            row = [
                result["id"][:20] + "..." if len(result["id"]) > 20 else result["id"],
                result["smiles"][:40] + "..." if len(result["smiles"]) > 40 else result["smiles"],
                status,
                error[:50] + "..." if len(error) > 50 else error
            ]
            
            if any(r.get("canonical_smiles") for r in results):
                canonical = result.get("canonical_smiles", "") or ""
                canonical_display = canonical[:40] + "..." if len(canonical) > 40 else canonical
                row.append(canonical_display)
            
            table.add_row(*row)
        
        console.print("\n")
        console.print(table)


def _save_validation_results(
    results: List[Dict[str, Any]],
    output_file: Path
):
    """Save validation results to file"""
    
    if output_file.suffix.lower() == '.csv':
        _save_results_csv(results, output_file)
    elif output_file.suffix.lower() == '.json':
        _save_results_json(results, output_file)
    else:
        # Default to CSV
        _save_results_csv(results, output_file)


def _save_results_csv(
    results: List[Dict[str, Any]],
    output_file: Path
):
    """Save results as CSV"""
    
    if not results:
        return
    
    fieldnames = ["id", "smiles", "valid", "error"]
    if any(r.get("canonical_smiles") for r in results):
        fieldnames.append("canonical_smiles")
    
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for result in results:
            # Only write fields that exist
            filtered_result = {k: v for k, v in result.items() if k in fieldnames}
            writer.writerow(filtered_result)


def _save_results_json(
    results: List[Dict[str, Any]],
    output_file: Path
):
    """Save results as JSON"""
    
    import json
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)