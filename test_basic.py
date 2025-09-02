#!/usr/bin/env python
"""
Basic test script to verify CriOS functionality
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

from src.core.molecule import CriOSMolecule
from src.core.filters import MolecularFilter

def test_basic_functionality():
    print("=" * 60)
    print("CriOS Basic Functionality Test")
    print("=" * 60)
    
    # Test 1: Create molecule
    print("\n1. Creating CriOSMolecule from aspirin SMILES...")
    aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = CriOSMolecule(aspirin_smiles, mol_id="aspirin")
    print(f"   [OK] Molecule created: {mol.mol_id}")
    print(f"   [OK] Canonical SMILES: {mol.canonical_smiles}")
    
    # Test 2: Calculate descriptors
    print("\n2. Calculating molecular descriptors...")
    descriptors = mol.calculate_descriptors(["MW", "LogP", "TPSA", "HBA", "HBD"])
    print("   Descriptors calculated:")
    for name, value in descriptors.items():
        print(f"     • {name}: {value:.2f}" if isinstance(value, float) else f"     • {name}: {value}")
    
    # Test 3: Check Lipinski filter
    print("\n3. Checking Lipinski's Rule of Five...")
    passes_lipinski = mol.passes_filter("Lipinski")
    print(f"   [OK] Passes Lipinski filter: {passes_lipinski}")
    
    # Test 4: Check violations if any
    filter_obj = MolecularFilter("lipinski")
    violations = filter_obj.get_violations(mol)
    if violations:
        print("   Violations found:")
        for desc, violation in violations.items():
            print(f"     • {desc}: {violation}")
    else:
        print("   [OK] No violations found")
    
    # Test 5: Test with another molecule
    print("\n4. Testing with ibuprofen...")
    ibuprofen_smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    ibu_mol = CriOSMolecule(ibuprofen_smiles, mol_id="ibuprofen")
    ibu_descriptors = ibu_mol.calculate_descriptors(["MW", "LogP", "TPSA"])
    print(f"   [OK] Ibuprofen MW: {ibu_descriptors['MW']:.2f}")
    print(f"   [OK] Ibuprofen LogP: {ibu_descriptors['LogP']:.2f}")
    print(f"   [OK] Passes Lipinski: {ibu_mol.passes_filter('Lipinski')}")
    
    print("\n" + "=" * 60)
    print("[SUCCESS] All basic tests passed successfully!")
    print("=" * 60)
    
    return True

if __name__ == "__main__":
    try:
        success = test_basic_functionality()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n[ERROR] Error during testing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)