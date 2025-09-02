from __future__ import annotations
import click, csv, sys
from pathlib import Path
from typing import Optional, List, Iterable, Tuple
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeRemainingColumn

from rdkit import Chem

from src.core.molecule import CriOSMolecule
from src.core.descriptors import DescriptorCalculator
from src.core.filters import MolecularFilter
from src.database import models as dbm
from src.database import api_clients as apis
from src.analysis.similarity import tanimoto_search
from src.analysis.clustering import butina_clusters

console = Console()

def _iter_smi_file(path: Path) -> Iterable[Tuple[str,str]]:
    with path.open("r", encoding="utf-8") as fh:
        for i, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"): continue
            parts = line.split()
            if len(parts) == 1:
                yield (f"MOL_{i:06d}", parts[0])
            else:
                yield (parts[1], parts[0])

@click.group(help="CriOS – Compound Discovery CLI")
def cli(): pass

@cli.command(help="Validate SMILES file")
@click.option("--input", "inp", type=click.Path(exists=True), required=True)
def validate(inp):
    inp = Path(inp)
    ok = 0; bad = 0
    for mol_id, smi in _iter_smi_file(inp):
        if Chem.MolFromSmiles(smi) is None:
            bad += 1
        else:
            ok += 1
    console.print(f"[bold]Valid:[/bold] {ok}   [bold red]Invalid:[/bold red] {bad}")
    sys.exit(0 if bad == 0 else 2)

@cli.command(help="Process SMILES -> descriptors CSV")
@click.option("--input", "inp", type=click.Path(exists=True), required=True)
@click.option("--output", "outp", type=click.Path(), required=True)
@click.option("--descriptors", default="MW,LogP,TPSA,HBA,HBD,RotatableBonds")
def process(inp, outp, descriptors):
    inp = Path(inp); outp = Path(outp)
    names = [s.strip() for s in descriptors.split(",") if s.strip()]
    calc = DescriptorCalculator(names)
    rows = []
    
    with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"),
                  BarColumn(), TimeRemainingColumn()) as prog:
        task = prog.add_task("Calculating descriptors", total=None)
        for mol_id, smi in _iter_smi_file(inp):
            try:
                mol = CriOSMolecule(smi, mol_id)
                d = calc.calculate(mol.mol)
                drow = {"mol_id": mol_id, "smiles": mol.canonical_smiles}
                drow.update(d)
                rows.append(drow)
            except Exception:
                continue
        prog.update(task, completed=True)
    
    outp.parent.mkdir(parents=True, exist_ok=True)
    with outp.open("w", newline="", encoding="utf-8") as fh:
        if rows:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    console.print(f"[green]Wrote[/green] {len(rows)} rows -> {outp}")

@cli.command(help="Filter SMILES by rule")
@click.option("--input", "inp", type=click.Path(exists=True), required=True)
@click.option("--rule", required=True, help="e.g., Lipinski OR 'MW<500 AND LogP<5'")
@click.option("--output", "outp", type=click.Path(), required=True)
@click.option("--save-rejected", type=click.Path(), default=None)
def filter(inp, rule, outp, save_rejected):
    filt = MolecularFilter(rule)
    kept = []; rej = []
    
    for mol_id, smi in _iter_smi_file(Path(inp)):
        try:
            m = CriOSMolecule(smi, mol_id)
            if filt.evaluate(m): 
                kept.append((mol_id, m.canonical_smiles))
            else: 
                rej.append((mol_id, m.canonical_smiles))
        except Exception: 
            pass
    
    outp = Path(outp); outp.parent.mkdir(parents=True, exist_ok=True)
    with outp.open("w", encoding="utf-8") as fh:
        for mol_id, smi in kept: 
            fh.write(f"{smi}\t{mol_id}\n")
    
    if save_rejected:
        with Path(save_rejected).open("w", encoding="utf-8") as fh:
            for mol_id, smi in rej: 
                fh.write(f"{smi}\t{mol_id}\n")
    
    console.print(f"[green]Kept[/green] {len(kept)}  [red]Rejected[/red] {len(rej)}")

@cli.command(name="db-import", help="Import SMILES to SQLite")
@click.option("--database", "db", type=click.Path(), required=True)
@click.option("--input", "inp", type=click.Path(exists=True), required=True)
def db_import(db, inp):
    con = dbm.connect(db)
    n = dbm.import_smiles(con, _iter_smi_file(Path(inp)))
    console.print(f"[green]Imported[/green] {n} molecules -> {db}; total={dbm.count(con)}")

@cli.command(help="Tanimoto similarity search")
@click.option("--database", "db", type=click.Path(exists=True), required=True)
@click.option("--query", "query_smiles", required=True)
@click.option("--threshold", type=float, default=0.7)
@click.option("--top-k", type=int, default=50)
def search(db, query_smiles, threshold, top_k):
    con = dbm.connect(db)
    hits = tanimoto_search(query_smiles, dbm.iter_all(con), threshold=threshold)
    hits = hits[:top_k]
    
    if hits:
        table = Table(title=f"Top {len(hits)} hits (Tanimoto >= {threshold})")
        table.add_column("mol_id", style="cyan")
        table.add_column("similarity", style="magenta")
        for mol_id, sim in hits:
            table.add_row(mol_id, f"{sim:.3f}")
        console.print(table)
    else:
        console.print(f"[yellow]No hits found with threshold >= {threshold}[/yellow]")

@cli.command(help="Cluster compounds with Butina")
@click.option("--input", "inp", type=click.Path(exists=True), required=True)
@click.option("--threshold", type=float, default=0.4)
@click.option("--output", "outp", type=click.Path(), required=True)
def cluster(inp, threshold, outp):
    smiles = [s for _, s in _iter_smi_file(Path(inp))]
    clusters = butina_clusters(smiles, threshold=threshold)
    
    outp = Path(outp); outp.parent.mkdir(parents=True, exist_ok=True)
    with outp.open("w", encoding="utf-8") as fh:
        for i, idxs in enumerate(clusters, 1):
            fh.write(f"# Cluster {i} (n={len(idxs)})\n")
            for j in idxs:
                fh.write(f"{smiles[j]}\n")
    console.print(f"[green]Wrote[/green] {len(clusters)} clusters -> {outp}")

@cli.command(help="Fetch from PubChem/ChEMBL")
@click.option("--source", type=click.Choice(["pubchem","chembl"]), required=True)
@click.option("--query", help="PubChem name query")
@click.option("--target", help="ChEMBL target ID")
@click.option("--limit", type=int, default=200)
@click.option("--output", "outp", type=click.Path(), required=True)
def fetch(source, query, target, limit, outp):
    rows = []
    if source == "pubchem":
        if not query: 
            raise click.ClickException("--query required for pubchem")
        rows = apis.pubchem_smiles(query, limit=limit)
    else:
        if not target: 
            raise click.ClickException("--target required for chembl")
        rows = apis.chembl_by_target(target, limit=limit)
    
    outp = Path(outp); outp.parent.mkdir(parents=True, exist_ok=True)
    with outp.open("w", encoding="utf-8") as fh:
        for mol_id, smi in rows:
            fh.write(f"{smi}\t{mol_id}\n")
    console.print(f"[green]Fetched[/green] {len(rows)} molecules -> {outp}")

if __name__ == '__main__':
    cli()