import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.core.molecule import CriOSMolecule
from src.core.descriptors import DescriptorCalculator
from src.core.filters import MolecularFilter
from src.core.validators import MoleculeValidator


class TestMoleculeProcessing:
    
    @pytest.fixture
    def valid_molecules(self):
        return [
            ("CC(=O)Oc1ccccc1C(=O)O", "aspirin"),
            ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "ibuprofen"),
            ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "caffeine"),
            ("CC(C)(C)NCC(O)c1ccc(O)c(O)c1", "epinephrine"),
            ("CC(=O)Nc1ccc(O)cc1", "acetaminophen"),
        ]
    
    @pytest.fixture
    def invalid_molecules(self):
        return [
            ("C1CCC", "unclosed_ring"),
            ("C(C)(C)(C)(C)C", "invalid_valence"),
            ("", "empty_string"),
            ("INVALID", "invalid_syntax"),
        ]
    
    def test_valid_molecule_creation(self, valid_molecules):
        for smiles, name in valid_molecules:
            mol = CriOSMolecule(smiles, mol_id=name)
            assert mol.mol is not None
            assert mol.mol_id == name
            assert mol.canonical_smiles is not None
    
    def test_descriptor_calculation(self, valid_molecules):
        for smiles, name in valid_molecules:
            mol = CriOSMolecule(smiles, mol_id=name)
            descriptors = mol.calculate_descriptors(["MW", "LogP", "TPSA", "HBA", "HBD"])
            
            assert "MW" in descriptors
            assert descriptors["MW"] > 0
            assert "LogP" in descriptors
            assert "TPSA" in descriptors
            assert "HBA" in descriptors
            assert isinstance(descriptors["HBA"], (int, float))
            assert "HBD" in descriptors
            assert isinstance(descriptors["HBD"], (int, float))
    
    def test_invalid_molecule_handling(self, invalid_molecules):
        for smiles, error_type in invalid_molecules:
            if smiles == "":
                with pytest.raises(ValueError):
                    mol = CriOSMolecule(smiles)
                    _ = mol.mol
            else:
                try:
                    mol = CriOSMolecule(smiles)
                    _ = mol.mol
                except ValueError:
                    pass
    
    def test_fingerprint_generation(self):
        mol = CriOSMolecule("CC(=O)Oc1ccccc1C(=O)O", "aspirin")
        
        fp_morgan = mol.get_fingerprint(fp_type="morgan", radius=2)
        assert fp_morgan is not None
        
        fp_rdkit = mol.get_fingerprint(fp_type="rdkit")
        assert fp_rdkit is not None
    
    def test_lipinski_filter(self):
        aspirin = CriOSMolecule("CC(=O)Oc1ccccc1C(=O)O", "aspirin")
        assert aspirin.passes_filter("Lipinski") == True
        
        caffeine = CriOSMolecule("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "caffeine")
        assert caffeine.passes_filter("Lipinski") == True
    
    def test_lead_like_filter(self):
        acetaminophen = CriOSMolecule("CC(=O)Nc1ccc(O)cc1", "acetaminophen")
        filter_obj = MolecularFilter("lead_like")
        assert filter_obj.apply(acetaminophen) == True
        
        testosterone = CriOSMolecule("CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C", "testosterone")
        assert filter_obj.apply(testosterone) == False


class TestDescriptorCalculator:
    
    def test_basic_descriptors(self):
        from rdkit import Chem
        
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        calc = DescriptorCalculator()
        
        descriptors = calc.calculate(mol, ["MW", "LogP", "TPSA"])
        
        assert "MW" in descriptors
        assert 180 < descriptors["MW"] < 181
        assert "LogP" in descriptors
        assert "TPSA" in descriptors
    
    def test_descriptor_sets(self):
        from rdkit import Chem
        
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        calc = DescriptorCalculator()
        
        lipinski_descs = calc.calculate(mol, calc.DESCRIPTOR_SETS["lipinski"])
        assert len(lipinski_descs) == 4
        assert all(k in lipinski_descs for k in ["MW", "LogP", "HBA", "HBD"])
    
    def test_batch_calculation(self):
        from rdkit import Chem
        
        smiles_list = [
            "CC(=O)Oc1ccccc1C(=O)O",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        ]
        
        mols = [Chem.MolFromSmiles(s) for s in smiles_list]
        calc = DescriptorCalculator()
        
        results = calc.calculate_batch(mols, ["MW", "LogP"])
        
        assert len(results) == 3
        for result in results:
            assert "MW" in result
            assert "LogP" in result


class TestMolecularFilter:
    
    def test_predefined_filters(self):
        filter_lipinski = MolecularFilter("lipinski")
        assert len(filter_lipinski.range_rules) == 4
        
        filter_lead = MolecularFilter("lead_like")
        assert len(filter_lead.range_rules) == 5
    
    def test_custom_rule_parser(self):
        filter_custom = MolecularFilter("custom")
        filter_custom.parse_rule_string("MW<500 AND LogP<5 AND HBD<=5")
        
        assert len(filter_custom.rules) == 3
        assert filter_custom.rules[0].descriptor == "MW"
        assert filter_custom.rules[0].operator == "<"
        assert filter_custom.rules[0].value == 500
    
    def test_filter_violations(self):
        large_mol = CriOSMolecule("C" * 50, "large")
        filter_obj = MolecularFilter("lipinski")
        
        violations = filter_obj.get_violations(large_mol)
        assert "MW" in violations


class TestMoleculeValidator:
    
    def test_smiles_validation(self):
        validator = MoleculeValidator(strict_mode=True)
        
        valid, error = validator.validate_smiles("CC(=O)Oc1ccccc1C(=O)O")
        assert valid == True
        assert error is None
        
        valid, error = validator.validate_smiles("INVALID")
        assert valid == False
        assert error is not None
    
    def test_standardization(self):
        validator = MoleculeValidator(standardize=True, remove_salts=True)
        
        canonical, error = validator.process_smiles("c1ccccc1")
        assert canonical == "c1ccccc1"
        assert error is None
        
        canonical, error = validator.process_smiles("C1=CC=CC=C1")
        assert canonical == "c1ccccc1"
    
    def test_batch_validation(self):
        validator = MoleculeValidator()
        
        smiles_list = [
            "CC(=O)Oc1ccccc1C(=O)O",
            "INVALID",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        ]
        
        results = validator.validate_batch(smiles_list)
        
        assert results[0][0] == True
        assert results[1][0] == False
        assert results[2][0] == True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])