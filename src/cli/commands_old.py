import click
import sys
from pathlib import Path
from typing import Optional, List
import logging
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeRemainingColumn
from rich.logging import RichHandler
import time

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.core import CriOSMolecule, DescriptorCalculator, MolecularFilter, MoleculeValidator
from src.io import ReaderFactory, WriterFactory, FormatConverter

console = Console()

def setup_logging(verbose: bool):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(message)s",
        handlers=[RichHandler(console=console, rich_tracebacks=True)]
    )

@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
def cli(verbose):
    setup_logging(verbose)

@cli.command()
@click.option('--input', '-i', 'input_file', required=True, type=click.Path(exists=True),
              help='Input file with molecules')
@click.option('--output', '-o', 'output_file', required=True, type=click.Path(),
              help='Output file for results')
@click.option('--descriptors', '-d', default='MW,LogP,TPSA,HBA,HBD',
              help='Comma-separated list of descriptors to calculate')
@click.option('--parallel', '-p', default=1, type=int,
              help='Number of parallel workers')
@click.option('--validate', '-val', default='strict', type=click.Choice(['strict', 'loose', 'none']),
              help='Validation level for molecules')
def process(input_file, output_file, descriptors, parallel, validate):
    console.print(f"[bold blue]Processing molecules from {input_file}[/bold blue]")
    
    descriptor_list = [d.strip() for d in descriptors.split(',')]
    console.print(f"Calculating descriptors: {', '.join(descriptor_list)}")
    
    try:
        reader = ReaderFactory.create_reader(input_file, validate=(validate != 'none'))
        calc = DescriptorCalculator()
        
        processed_molecules = []
        errors = 0
        
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeRemainingColumn(),
            console=console
        ) as progress:
            
            molecules = list(reader.read())
            task = progress.add_task("Processing molecules...", total=len(molecules))
            
            for mol_data in molecules:
                try:
                    mol = CriOSMolecule(mol_data['smiles'], mol_data['mol_id'])
                    mol_descriptors = mol.calculate_descriptors(descriptor_list)
                    
                    mol_data['descriptors'] = mol_descriptors
                    processed_molecules.append(mol_data)
                    
                except Exception as e:
                    logging.error(f"Error processing {mol_data.get('mol_id', 'unknown')}: {e}")
                    errors += 1
                
                progress.update(task, advance=1)
        
        writer = WriterFactory.create_writer(output_file)
        written = writer.write(processed_molecules)
        
        console.print(f"[green]✓[/green] Processed {len(processed_molecules)} molecules")
        console.print(f"[green]✓[/green] Written {written} molecules to {output_file}")
        if errors > 0:
            console.print(f"[yellow]⚠[/yellow] {errors} molecules had errors")
        
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        sys.exit(1)

@cli.command()
@click.option('--input', '-i', 'input_file', required=True, type=click.Path(exists=True),
              help='Input file with molecules')
@click.option('--output', '-o', 'output_file', required=True, type=click.Path(),
              help='Output file for filtered molecules')
@click.option('--rules', '-r', help='Filter rules (e.g., "MW<500 AND LogP<5")')
@click.option('--filter-type', '-f', default='lipinski',
              type=click.Choice(['lipinski', 'lead_like', 'fragment', 'drug_like', 'custom']),
              help='Predefined filter type')
@click.option('--save-rejected', '-rej', type=click.Path(),
              help='Save rejected molecules to this file')
def filter(input_file, output_file, rules, filter_type, save_rejected):
    console.print(f"[bold blue]Filtering molecules from {input_file}[/bold blue]")
    
    try:
        if filter_type == 'custom' and not rules:
            console.print("[red]Error:[/red] Custom filter requires --rules parameter")
            sys.exit(1)
        
        mol_filter = MolecularFilter(filter_type)
        if rules:
            mol_filter.parse_rule_string(rules)
        
        console.print(f"Using filter: {filter_type}")
        if rules:
            console.print(f"Custom rules: {rules}")
        
        reader = ReaderFactory.create_reader(input_file)
        calc = DescriptorCalculator()
        
        passed_molecules = []
        rejected_molecules = []
        
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            
            molecules = list(reader.read())
            task = progress.add_task("Filtering molecules...", total=len(molecules))
            
            for mol_data in molecules:
                try:
                    mol = CriOSMolecule(mol_data['smiles'], mol_data['mol_id'])
                    
                    if mol_filter.apply(mol):
                        passed_molecules.append(mol_data)
                    else:
                        rejected_molecules.append(mol_data)
                    
                except Exception as e:
                    logging.error(f"Error filtering {mol_data.get('mol_id', 'unknown')}: {e}")
                    rejected_molecules.append(mol_data)
                
                progress.update(task, advance=1)
        
        writer = WriterFactory.create_writer(output_file)
        written = writer.write(passed_molecules)
        
        console.print(f"[green]✓[/green] {len(passed_molecules)} molecules passed filter")
        console.print(f"[yellow]✗[/yellow] {len(rejected_molecules)} molecules rejected")
        console.print(f"[green]✓[/green] Written {written} molecules to {output_file}")
        
        if save_rejected and rejected_molecules:
            rej_writer = WriterFactory.create_writer(save_rejected)
            rej_written = rej_writer.write(rejected_molecules)
            console.print(f"[green]✓[/green] Written {rej_written} rejected molecules to {save_rejected}")
        
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        sys.exit(1)

@cli.command()
@click.option('--input', '-i', 'input_file', required=True, type=click.Path(exists=True),
              help='Input file with SMILES')
@click.option('--output', '-o', 'output_file', type=click.Path(),
              help='Output file for validation results')
@click.option('--strict', is_flag=True, help='Use strict validation')
@click.option('--standardize', is_flag=True, help='Standardize molecules')
@click.option('--remove-salts', is_flag=True, help='Remove salts from molecules')
def validate(input_file, output_file, strict, standardize, remove_salts):
    console.print(f"[bold blue]Validating molecules from {input_file}[/bold blue]")
    
    try:
        validator = MoleculeValidator(
            strict_mode=strict,
            standardize=standardize,
            remove_salts=remove_salts
        )
        
        reader = ReaderFactory.create_reader(input_file, validate=False)
        
        results = []
        valid_count = 0
        invalid_count = 0
        
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console
        ) as progress:
            
            task = progress.add_task("Validating molecules...")
            
            for mol_data in reader.read():
                smiles = mol_data['smiles']
                mol_id = mol_data['mol_id']
                
                is_valid, error = validator.validate_smiles(smiles, mol_id)
                
                if is_valid:
                    valid_count += 1
                    if standardize:
                        processed_smiles, _ = validator.process_smiles(smiles, mol_id)
                        mol_data['smiles'] = processed_smiles
                else:
                    invalid_count += 1
                
                results.append({
                    'mol_id': mol_id,
                    'smiles': smiles,
                    'valid': is_valid,
                    'error': error
                })
        
        table = Table(title="Validation Results")
        table.add_column("Status", style="cyan")
        table.add_column("Count", style="magenta")
        table.add_column("Percentage", style="green")
        
        total = valid_count + invalid_count
        table.add_row("Valid", str(valid_count), f"{valid_count/total*100:.1f}%")
        table.add_row("Invalid", str(invalid_count), f"{invalid_count/total*100:.1f}%")
        table.add_row("Total", str(total), "100.0%")
        
        console.print(table)
        
        if output_file:
            import pandas as pd
            df = pd.DataFrame(results)
            df.to_csv(output_file, index=False)
            console.print(f"[green]✓[/green] Results saved to {output_file}")
        
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        sys.exit(1)

@cli.command()
@click.option('--from', '-f', 'input_format', type=click.Choice(['smi', 'sdf', 'csv']),
              help='Input format (auto-detected if not specified)')
@click.option('--to', '-t', 'output_format', type=click.Choice(['smi', 'sdf', 'csv']),
              help='Output format (auto-detected if not specified)')
@click.option('--input', '-i', 'input_file', required=True, type=click.Path(exists=True),
              help='Input file')
@click.option('--output', '-o', 'output_file', required=True, type=click.Path(),
              help='Output file')
@click.option('--canonical', is_flag=True, default=True, help='Canonicalize SMILES')
@click.option('--remove-salts', is_flag=True, help='Remove salts from molecules')
@click.option('--remove-duplicates', is_flag=True, help='Remove duplicate molecules')
def convert(input_format, output_format, input_file, output_file, canonical, remove_salts, remove_duplicates):
    console.print(f"[bold blue]Converting {input_file} to {output_file}[/bold blue]")
    
    try:
        converter = FormatConverter(
            canonical=canonical,
            remove_duplicates=remove_duplicates
        )
        
        result = converter.convert(
            input_file,
            output_file,
            input_format=input_format,
            output_format=output_format,
            show_progress=True
        )
        
        console.print(f"[green]✓[/green] Processed {result['molecules_processed']} molecules")
        console.print(f"[green]✓[/green] Written {result['molecules_written']} molecules")
        
        if result['errors'] > 0:
            console.print(f"[yellow]⚠[/yellow] {result['errors']} errors occurred")
        
        if result['duplicates_removed'] > 0:
            console.print(f"[yellow]ℹ[/yellow] {result['duplicates_removed']} duplicates removed")
        
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        sys.exit(1)

@cli.command()
@click.option('--source', '-s', type=click.Choice(['pubchem', 'chembl', 'all']), 
              default='pubchem', help='Database source')
@click.option('--query', '-q', required=True, help='Search query (compound name or disease)')
@click.option('--limit', '-l', default=100, type=int, help='Maximum results')
@click.option('--output', '-o', type=click.Path(), help='Output file for results')
@click.option('--target', '-t', help='ChEMBL target ID (e.g., CHEMBL203)')
def fetch(source, query, limit, output, target):
    """Fetch compounds from public databases"""
    console.print(f"[bold blue]Fetching compounds from {source}...[/bold blue]")
    
    try:
        from src.database.api_clients import CompoundDatabase
        
        db = CompoundDatabase()
        
        if target and source == 'chembl':
            console.print(f"Searching ChEMBL for target: {target}")
            compounds = db.chembl.search_by_target(target, limit)
        else:
            console.print(f"Searching for: {query}")
            compounds = db.search(query, source, limit)
        
        console.print(f"[green]✓[/green] Found {len(compounds)} compounds")
        
        if output:
            import pandas as pd
            df = pd.DataFrame(compounds)
            if output.endswith('.csv'):
                df.to_csv(output, index=False)
            else:
                # Save as SMILES file
                with open(output, 'w') as f:
                    for comp in compounds:
                        if comp.get('smiles'):
                            f.write(f"{comp['smiles']}\t{comp.get('chembl_id', comp.get('cid', 'Unknown'))}\n")
            
            console.print(f"[green]✓[/green] Results saved to {output}")
        else:
            # Display results in table
            table = Table(title=f"Search Results from {source}")
            table.add_column("ID", style="cyan")
            table.add_column("Name", style="magenta")
            table.add_column("MW", style="green")
            table.add_column("Source", style="yellow")
            
            for comp in compounds[:10]:  # Show first 10
                table.add_row(
                    str(comp.get('chembl_id', comp.get('cid', ''))),
                    str(comp.get('name', comp.get('iupac', '')))[:30],
                    str(comp.get('mw', comp.get('molecular_weight', '')))[:10],
                    comp.get('source', '')
                )
            
            console.print(table)
            if len(compounds) > 10:
                console.print(f"[dim]... and {len(compounds) - 10} more compounds[/dim]")
        
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        sys.exit(1)

@cli.command()
@click.option('--query', '-q', required=True, help='Query SMILES or file')
@click.option('--database', '-d', required=True, type=click.Path(exists=True),
              help='Database file to search')
@click.option('--threshold', '-t', default=0.7, type=float,
              help='Similarity threshold (0-1)')
@click.option('--method', '-m', default='morgan',
              type=click.Choice(['morgan', 'rdkit', 'maccs']),
              help='Fingerprint method')
@click.option('--top-k', '-k', default=10, type=int, help='Number of results')
@click.option('--output', '-o', type=click.Path(), help='Output file')
def search(query, database, threshold, method, top_k, output):
    """Search for similar compounds"""
    console.print(f"[bold blue]Searching for similar compounds...[/bold blue]")
    
    try:
        from src.analysis.similarity import CompoundDiscovery
        from src.io.readers import ReaderFactory
        
        # Read query
        if Path(query).exists():
            reader = ReaderFactory.create_reader(query)
            query_mols = list(reader.read())
            query_smiles = query_mols[0]['smiles']
        else:
            query_smiles = query
        
        console.print(f"Query: {query_smiles[:50]}...")
        
        # Read database
        reader = ReaderFactory.create_reader(database)
        db_compounds = list(reader.read())
        db_smiles = [c['smiles'] for c in db_compounds]
        
        console.print(f"Searching {len(db_smiles)} compounds...")
        
        # Perform similarity search
        discovery = CompoundDiscovery()
        similar = discovery.find_similar_compounds(query_smiles, db_smiles, threshold)
        
        # Take top K results
        similar = similar[:top_k]
        
        console.print(f"[green]✓[/green] Found {len(similar)} similar compounds")
        
        if output:
            # Save results
            with open(output, 'w') as f:
                f.write("smiles\tsimilarity\tid\n")
                for smiles, score in similar:
                    f.write(f"{smiles}\t{score:.3f}\tSIM_{score:.0f}\n")
            console.print(f"[green]✓[/green] Results saved to {output}")
        else:
            # Display results
            table = Table(title="Similar Compounds")
            table.add_column("Rank", style="cyan")
            table.add_column("Similarity", style="magenta")
            table.add_column("SMILES", style="green")
            
            for i, (smiles, score) in enumerate(similar, 1):
                table.add_row(str(i), f"{score:.3f}", smiles[:50] + "..." if len(smiles) > 50 else smiles)
            
            console.print(table)
        
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        sys.exit(1)

@cli.command()
@click.option('--input', '-i', 'input_file', required=True, type=click.Path(exists=True),
              help='Input file with compounds')
@click.option('--method', '-m', default='butina',
              type=click.Choice(['butina', 'hierarchical']),
              help='Clustering method')
@click.option('--threshold', '-t', default=0.4, type=float,
              help='Distance threshold for Butina clustering')
@click.option('--n-clusters', '-n', default=10, type=int,
              help='Number of clusters for hierarchical')
@click.option('--output', '-o', type=click.Path(), help='Output file')
def cluster(input_file, method, threshold, n_clusters, output):
    """Cluster compound library"""
    console.print(f"[bold blue]Clustering compounds using {method} method...[/bold blue]")
    
    try:
        from src.analysis.clustering import CompoundClusterer
        from src.io.readers import ReaderFactory
        from rdkit import Chem
        
        # Read compounds
        reader = ReaderFactory.create_reader(input_file)
        compounds = list(reader.read())
        mols = [Chem.MolFromSmiles(c['smiles']) for c in compounds if c.get('smiles')]
        mols = [m for m in mols if m is not None]
        
        console.print(f"Clustering {len(mols)} compounds...")
        
        # Perform clustering
        clusterer = CompoundClusterer()
        
        if method == 'butina':
            clusters = clusterer.butina_cluster(mols, cutoff=threshold)
        else:
            clusters = clusterer.hierarchical_cluster(mols, n_clusters=n_clusters)
        
        # Get statistics
        stats = clusterer.cluster_statistics(mols, clusters)
        
        console.print(f"[green]✓[/green] Created {stats['n_clusters']} clusters")
        console.print(f"  Largest cluster: {stats['largest_cluster']} compounds")
        console.print(f"  Singletons: {stats['n_singletons']}")
        console.print(f"  Mean cluster size: {stats['mean_cluster_size']:.1f}")
        
        if output:
            # Save clustering results
            import json
            results = {
                'method': method,
                'parameters': {'threshold': threshold} if method == 'butina' else {'n_clusters': n_clusters},
                'statistics': stats,
                'clusters': [[compounds[i]['mol_id'] for i in cluster] for cluster in clusters]
            }
            
            with open(output, 'w') as f:
                json.dump(results, f, indent=2)
            
            console.print(f"[green]✓[/green] Results saved to {output}")
        
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        sys.exit(1)

if __name__ == '__main__':
    cli()