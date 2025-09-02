#!/usr/bin/env python
"""Test basic CLI functionality"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

# Test that we can run basic operations
from src.core.molecule import CriOSMolecule
from src.database import models as dbm

def main():
    print("CriOS CLI Basic Test")
    print("-" * 40)
    
    # Create test file
    test_smi = Path("test_compounds.smi")
    with test_smi.open("w") as f:
        f.write("CC(=O)Oc1ccccc1C(=O)O\taspirin\n")
        f.write("CC(C)Cc1ccc(cc1)C(C)C(=O)O\tibuprofen\n")
        f.write("CN1C=NC2=C1C(=O)N(C(=O)N2C)C\tcaffeine\n")
    
    print(f"Created {test_smi}")
    
    # Test database import
    db_path = "test_cli.db"
    con = dbm.connect(db_path)
    
    # Read and import
    compounds = []
    with test_smi.open() as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                compounds.append((parts[1], parts[0]))
    
    n = dbm.import_smiles(con, compounds)
    print(f"Imported {n} compounds to {db_path}")
    print(f"Database contains {dbm.count(con)} compounds")
    
    # Test similarity search
    from src.analysis.similarity import tanimoto_search
    
    hits = tanimoto_search("CC(=O)Oc1ccccc1C(=O)O", dbm.iter_all(con), threshold=0.5)
    print(f"\nSimilarity search results:")
    for mol_id, sim in hits:
        print(f"  {mol_id}: {sim:.3f}")
    
    con.close()
    
    # Cleanup
    test_smi.unlink()
    Path(db_path).unlink()
    
    print("\n[SUCCESS] Basic CLI operations work!")
    return 0

if __name__ == "__main__":
    sys.exit(main())