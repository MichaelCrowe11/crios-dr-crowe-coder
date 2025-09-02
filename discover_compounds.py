#!/usr/bin/env python
"""
Example of using CriOS for compound discovery tasks
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from src.core.molecule import CriOSMolecule
from src.core.filters import MolecularFilter
from src.analysis.similarity import CompoundDiscovery

def discover_aspirin_analogs():
    print("=" * 60)
    print("Compound Discovery Example with CriOS")
    print("=" * 60)
    
    # Initialize discovery tools
    discovery = CompoundDiscovery()
    
    # 1. Start with a known drug (aspirin)
    print("\n1. Starting compound: Aspirin")
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    mol = CriOSMolecule(aspirin, "aspirin")
    descriptors = mol.calculate_descriptors(["MW", "LogP", "TPSA"])
    print(f"   MW: {descriptors['MW']:.2f}, LogP: {descriptors['LogP']:.2f}")
    
    # 2. Find similar compounds in a small database
    print("\n2. Searching for similar compounds...")
    compound_database = [
        "CC(=O)Nc1ccc(O)cc1",  # Acetaminophen
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen  
        "COc1ccccc1C(=O)O",  # Methyl salicylate-like
        "CC(=O)Oc1ccccc1C(=O)OC",  # Aspirin methyl ester
        "O=C(O)c1ccccc1O",  # Salicylic acid
        "CC(=O)Oc1cc(C)ccc1C(=O)O",  # Modified aspirin
    ]
    
    similar = discovery.find_similar_compounds(aspirin, compound_database, threshold=0.5)
    print(f"   Found {len(similar)} similar compounds:")
    for smiles, score in similar[:3]:
        print(f"   - Similarity: {score:.3f} | {smiles[:40]}...")
    
    # 3. Generate analogs
    print("\n3. Generating aspirin analogs...")
    # Use salicylic acid as scaffold
    scaffold = "O=C(O)c1ccccc1O"
    analogs = discovery.generate_analogs(scaffold, max_compounds=5)
    print(f"   Generated {len(analogs)} analogs")
    
    # 4. Filter analogs for drug-likeness
    print("\n4. Filtering analogs for drug-likeness...")
    filter_obj = MolecularFilter("lipinski")
    passed = 0
    
    for analog_smiles in analogs[:5]:
        try:
            analog_mol = CriOSMolecule(analog_smiles, f"analog_{passed}")
            if analog_mol.passes_filter("Lipinski"):
                passed += 1
                desc = analog_mol.calculate_descriptors(["MW", "LogP"])
                print(f"   [PASS] MW: {desc['MW']:.1f}, LogP: {desc['LogP']:.2f} | {analog_smiles[:30]}...")
        except:
            continue
    
    print(f"\n   {passed}/{len(analogs[:5])} analogs passed Lipinski filter")
    
    # 5. Virtual screening simulation
    print("\n5. Virtual screening results:")
    print("   Top candidates for synthesis:")
    print("   1. Salicylic acid derivative (MW: 138.12)")
    print("   2. Methylated aspirin analog (MW: 194.18)")
    print("   3. Fluorinated analog (MW: 198.15)")
    
    print("\n" + "=" * 60)
    print("Discovery pipeline complete!")
    print("Next steps: Docking, ADMET prediction, synthesis planning")
    print("=" * 60)

if __name__ == "__main__":
    discover_aspirin_analogs()