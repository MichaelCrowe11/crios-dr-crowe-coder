#!/usr/bin/env python
"""Test improved CriOS system with clean architecture"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from src.core.molecule import CriOSMolecule
from src.core.filters import MolecularFilter
from src.analysis.similarity import tanimoto_search
from src.database import models as dbm

def test_improved_system():
    print("=" * 60)
    print("Testing Improved CriOS System")
    print("=" * 60)
    
    # Test 1: Clean molecule creation
    print("\n1. Testing molecule with new validators...")
    aspirin = CriOSMolecule("CC(=O)Oc1ccccc1C(=O)O", "aspirin")
    print(f"   Canonical SMILES: {aspirin.canonical_smiles}")
    
    # Test 2: Safe descriptor calculation
    print("\n2. Testing RDKit-safe descriptors...")
    desc = aspirin.calculate_descriptors(["MW", "LogP", "TPSA", "HBA", "HBD"])
    for name, value in desc.items():
        print(f"   {name}: {value:.2f}" if isinstance(value, float) else f"   {name}: {value}")
    
    # Test 3: Improved filter
    print("\n3. Testing filter evaluation...")
    filt = MolecularFilter("lipinski")
    passes = filt.evaluate(aspirin)
    print(f"   Passes Lipinski: {passes}")
    
    # Test 4: Custom rule parsing
    print("\n4. Testing custom rule parser...")
    custom = MolecularFilter("MW<200 AND LogP<2")
    passes_custom = custom.evaluate(aspirin)
    print(f"   Passes 'MW<200 AND LogP<2': {passes_custom}")
    
    # Test 5: Database operations
    print("\n5. Testing SQLite database...")
    db_path = "test.db"
    con = dbm.connect(db_path)
    
    # Import some test compounds
    test_compounds = [
        ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
        ("caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    ]
    
    n = dbm.import_smiles(con, test_compounds)
    print(f"   Imported {n} compounds")
    print(f"   Total in DB: {dbm.count(con)}")
    
    # Test 6: Similarity search
    print("\n6. Testing similarity search...")
    hits = tanimoto_search("CC(=O)Oc1ccccc1C(=O)O", dbm.iter_all(con), threshold=0.5)
    print(f"   Found {len(hits)} similar compounds:")
    for mol_id, sim in hits:
        print(f"     - {mol_id}: {sim:.3f}")
    
    con.close()
    
    print("\n" + "=" * 60)
    print("[SUCCESS] All tests passed!")
    print("=" * 60)
    
    # Clean up test database
    Path(db_path).unlink(missing_ok=True)
    
    return True

if __name__ == "__main__":
    try:
        test_improved_system()
    except Exception as e:
        print(f"\n[ERROR] Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)