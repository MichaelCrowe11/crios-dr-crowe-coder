from __future__ import annotations
import sqlite3
from typing import Iterable, Tuple, List, Optional
from rdkit import Chem

SCHEMA = """
CREATE TABLE IF NOT EXISTS compounds (
  id INTEGER PRIMARY KEY,
  mol_id TEXT UNIQUE,
  smiles TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_smiles ON compounds(smiles);
"""

def connect(db_path: str) -> sqlite3.Connection:
    con = sqlite3.connect(db_path)
    con.execute("PRAGMA journal_mode=WAL;")
    con.executescript(SCHEMA)
    return con

def import_smiles(con: sqlite3.Connection, rows: Iterable[Tuple[str,str]], batch: int = 1000) -> int:
    cur = con.cursor()
    count = 0; buf: List[Tuple[str,str]] = []
    for mol_id, smiles in rows:
        # ensure valid
        if Chem.MolFromSmiles(smiles) is None: continue
        buf.append((mol_id, smiles))
        if len(buf) >= batch:
            cur.executemany("INSERT OR IGNORE INTO compounds(mol_id,smiles) VALUES (?,?)", buf)
            con.commit(); count += len(buf); buf.clear()
    if buf:
        cur.executemany("INSERT OR IGNORE INTO compounds(mol_id,smiles) VALUES (?,?)", buf)
        con.commit(); count += len(buf)
    return count

def iter_all(con: sqlite3.Connection) -> Iterable[Tuple[str,str]]:
    for mol_id, smiles in con.execute("SELECT mol_id, smiles FROM compounds"):
        yield mol_id, smiles

def count(con: sqlite3.Connection) -> int:
    return int(con.execute("SELECT COUNT(*) FROM compounds").fetchone()[0])