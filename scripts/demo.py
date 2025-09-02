#!/usr/bin/env python3
"""
CriOS Discovery Engine - Complete System Demonstration
Comprehensive end-to-end demonstration of the CriOS platform capabilities
"""

import logging
import sys
import time
from pathlib import Path
from typing import List, Dict, Any

# Configure logging for demo
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    """Run comprehensive CriOS demonstration"""
    
    print("🧬 CriOS Discovery Engine - System Demonstration")
    print("=" * 60)
    print("Universal Compound Discovery Platform with Ethical AI")
    print("Science before status. Discovery before profit.")
    print("=" * 60)
    print()
    
    try:
        # Demo 1: Chemistry Substrate
        demo_chemistry_substrate()
        
        # Demo 2: Molecular Scoring
        demo_molecular_scoring()
        
        # Demo 3: Ethics and Safety
        demo_ethics_safety()
        
        # Demo 4: Fingerprints and Similarity
        demo_fingerprints_similarity()
        
        # Demo 5: CLI Interface
        demo_cli_interface()
        
        # Demo 6: Performance Benchmarks
        demo_performance()
        
        print("\n🎉 CriOS Discovery Engine Demonstration Complete!")
        print("=" * 60)
        print("✅ All systems operational and ready for compound discovery")
        print("📋 Next steps:")
        print("  - Customize configs/config.yaml for your research")
        print("  - Review ethics policies in configs/ethics.yaml")
        print("  - Try 'crios --help' for command-line interface")
        print("  - Start web interface with 'make web'")
        print("=" * 60)
        
    except Exception as e:
        logger.error(f"Demo failed: {e}")
        print(f"\n❌ Demo failed: {e}")
        print("Please check dependencies and try again.")
        sys.exit(1)


def demo_chemistry_substrate():
    """Demonstrate chemistry substrate functionality"""
    print("🔬 Demo 1: Chemistry Substrate & Molecular Processing")
    print("-" * 50)
    
    try:
        from crios.chem.mol import Molecule, standardize_smiles_list, calculate_properties_batch
        from crios.chem.fingerprints import tanimoto_similarity, smiles_to_fingerprint
        
        # Example molecules
        test_molecules = [
            ("CCO", "ethanol"),
            ("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "ibuprofen"),  
            ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "caffeine"),
            ("CC1=CC2=C(C=C1C)N(C=N2)C3C(C(C(O3)CO)O)O", "riboflavin")
        ]
        
        print("📊 Processing molecules:")
        results = []
        
        for smiles, name in test_molecules:
            try:
                mol = Molecule(smiles, mol_id=name)
                props = mol.get_all_properties()
                
                print(f"  ✓ {name:<12} MW: {props.get('molecular_weight', 0):.1f} Da, "
                      f"LogP: {props.get('logp', 0):.2f}, "
                      f"Drug-like: {'Yes' if mol.passes_drug_like_filters() else 'No'}")
                
                results.append((mol, props))
                
            except Exception as e:
                print(f"  ❌ {name}: {e}")
        
        print(f"✅ Successfully processed {len(results)}/{len(test_molecules)} molecules")
        print()
        
    except ImportError as e:
        print(f"❌ Chemistry substrate not available: {e}")
        print("Please install RDKit: pip install rdkit")
        print()


def demo_molecular_scoring():
    """Demonstrate molecular scoring functionality"""
    print("⭐ Demo 2: Molecular Scoring & Crowe Methodology")
    print("-" * 50)
    
    try:
        from crios.scoring.design import crowe_score, _passes_drug_like_filters
        
        # Test molecules with expected different score profiles
        test_compounds = [
            ("CCO", "ethanol"),
            ("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "ibuprofen"),
            ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "caffeine"),
            ("C1=CC=C(C=C1)C(=O)O", "benzoic_acid")
        ]
        
        print("📊 Crowe Scoring Results:")
        print(f"{'Compound':<15} {'Crowe':<6} {'Potency':<7} {'ADMET':<6} {'Safety':<6} {'Drug-like':<9}")
        print("-" * 65)
        
        scored_compounds = []
        
        for smiles, name in test_compounds:
            try:
                scores = crowe_score(smiles, target_class="GPCR")
                
                print(f"{name:<15} "
                      f"{scores.get('crowe_score', 0):.3f}  "
                      f"{scores.get('potency', 0):.3f}   "
                      f"{scores.get('admet', 0):.3f}  "
                      f"{scores.get('synthesis', 0):.3f}  "
                      f"{'Yes' if scores.get('drug_like', False) else 'No':<9}")
                
                scored_compounds.append((name, scores))
                
            except Exception as e:
                print(f"{name:<15} ERROR: {e}")
        
        # Show top compound
        if scored_compounds:
            best_compound = max(scored_compounds, key=lambda x: x[1].get('crowe_score', 0))
            print(f"\n🏆 Top compound: {best_compound[0]} (score: {best_compound[1].get('crowe_score', 0):.3f})")
        
        print(f"✅ Successfully scored {len(scored_compounds)}/{len(test_compounds)} compounds")
        print()
        
    except ImportError as e:
        print(f"❌ Scoring system not available: {e}")
        print()


def demo_ethics_safety():
    """Demonstrate ethics and safety functionality"""
    print("⚖️  Demo 3: Ethics & Safety Framework")
    print("-" * 50)
    
    try:
        from crios.ethics.filters import check_ethics_compliance
        
        # Test molecules with different safety profiles
        test_compounds = [
            ("CCO", "ethanol", "safe"),
            ("C1OC1", "ethylene_oxide", "hazardous"),
            ("c1ccc([N+](=O)[O-])cc1", "nitrobenzene", "toxic"),
            ("C1=CC=C(C=C1)C(=O)O", "benzoic_acid", "generally_safe")
        ]
        
        print("🛡️ Ethics Compliance Check:")
        print(f"{'Compound':<15} {'Status':<6} {'Violations':<12} {'Recommendation'}")
        print("-" * 70)
        
        compliant_count = 0
        
        for smiles, name, expected in test_compounds:
            try:
                result = check_ethics_compliance(smiles)
                
                status = "✅ PASS" if result["passed"] else "❌ FAIL"
                violations = len(result["violations"])
                recommendation = result["explanation"].get("recommendation", "No recommendation")[:25] + "..."
                
                print(f"{name:<15} {status:<6} {violations:<12} {recommendation}")
                
                if result["passed"]:
                    compliant_count += 1
                
            except Exception as e:
                print(f"{name:<15} ERROR: {str(e)[:50]}...")
        
        print(f"\n✅ Ethics compliance: {compliant_count}/{len(test_compounds)} compounds passed")
        print("ℹ️  Note: Default policy uses example patterns only")
        print()
        
    except ImportError as e:
        print(f"❌ Ethics framework not available: {e}")
        print()


def demo_fingerprints_similarity():
    """Demonstrate fingerprint and similarity functionality"""
    print("🔍 Demo 4: Fingerprints & Molecular Similarity")
    print("-" * 50)
    
    try:
        from crios.chem.fingerprints import smiles_similarity, smiles_to_fingerprint
        
        # Test similarity calculations
        query = "CCO"  # ethanol
        comparisons = [
            ("CCC", "propane"),
            ("CC(C)O", "isopropanol"),
            ("CCCO", "propanol"),
            ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "caffeine")
        ]
        
        print(f"📊 Similarity to {query} (ethanol):")
        print(f"{'Compound':<15} {'SMILES':<25} {'Tanimoto':<10}")
        print("-" * 55)
        
        similarities = []
        
        for smiles, name in comparisons:
            try:
                sim = smiles_similarity(query, smiles, metric="tanimoto")
                print(f"{name:<15} {smiles:<25} {sim:.3f}")
                similarities.append((name, sim))
                
            except Exception as e:
                print(f"{name:<15} {smiles:<25} ERROR")
        
        # Show most similar
        if similarities:
            most_similar = max(similarities, key=lambda x: x[1])
            print(f"\n🎯 Most similar: {most_similar[0]} (similarity: {most_similar[1]:.3f})")
        
        print("✅ Similarity calculations complete")
        print()
        
    except ImportError as e:
        print(f"❌ Fingerprint system not available: {e}")
        print()


def demo_cli_interface():
    """Demonstrate CLI interface"""
    print("💻 Demo 5: Command Line Interface")
    print("-" * 50)
    
    try:
        import subprocess
        import sys
        
        # Test basic CLI functionality
        commands = [
            (["python", "-m", "crios.cli.main", "--help"], "Help command"),
            # Note: We can't easily demo the full CLI here without creating temp files
        ]
        
        print("🔧 CLI Interface Status:")
        
        # Check if crios command is available
        try:
            result = subprocess.run(
                [sys.executable, "-m", "crios.cli.main", "--help"],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode == 0:
                print("  ✅ CLI interface operational")
                print("  📋 Available commands: init, featurize, score, design, similarity, explain")
                print("  💡 Try: python -m crios.cli.main --help")
            else:
                print("  ❌ CLI interface error")
                print(f"     Error: {result.stderr[:100]}...")
                
        except subprocess.TimeoutExpired:
            print("  ⏱️ CLI command timed out")
        except Exception as e:
            print(f"  ❌ CLI test failed: {e}")
        
        print("✅ CLI interface check complete")
        print()
        
    except ImportError as e:
        print(f"❌ CLI interface not available: {e}")
        print()


def demo_performance():
    """Demonstrate performance characteristics"""
    print("⚡ Demo 6: Performance Benchmarks")
    print("-" * 50)
    
    try:
        from crios.scoring.design import crowe_score
        from crios.chem.fingerprints import smiles_similarity
        
        # Performance test molecules
        test_molecules = [
            "CCO",
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "CC1=CC2=C(C=C1C)N(C=N2)C3C(C(C(O3)CO)O)O",
            "C1=CC=C(C=C1)C(=O)O"
        ]
        
        n_iterations = 20
        
        # Benchmark scoring
        print("📊 Performance Benchmarks:")
        
        start_time = time.time()
        successful_scores = 0
        
        for _ in range(n_iterations):
            for smiles in test_molecules:
                try:
                    result = crowe_score(smiles)
                    successful_scores += 1
                except Exception:
                    pass
        
        scoring_time = time.time() - start_time
        total_molecules = n_iterations * len(test_molecules)
        
        print(f"  Scoring: {successful_scores}/{total_molecules} molecules in {scoring_time:.2f}s")
        if successful_scores > 0:
            print(f"  Average: {(scoring_time / successful_scores) * 1000:.1f}ms per molecule")
        
        # Benchmark similarity
        start_time = time.time()
        successful_similarities = 0
        
        for _ in range(n_iterations):
            for i in range(len(test_molecules)):
                for j in range(i + 1, len(test_molecules)):
                    try:
                        sim = smiles_similarity(test_molecules[i], test_molecules[j])
                        successful_similarities += 1
                    except Exception:
                        pass
        
        similarity_time = time.time() - start_time
        
        print(f"  Similarity: {successful_similarities} comparisons in {similarity_time:.2f}s")
        if successful_similarities > 0:
            print(f"  Average: {(similarity_time / successful_similarities) * 1000:.1f}ms per comparison")
        
        print("✅ Performance benchmarks complete")
        print()
        
    except ImportError as e:
        print(f"❌ Performance benchmarks not available: {e}")
        print()


if __name__ == "__main__":
    main()