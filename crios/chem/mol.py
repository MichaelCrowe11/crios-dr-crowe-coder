"""
CriOS Chemistry Substrate - Molecular Operations
RDKit wrappers for standardization, validation, and molecular processing
"""

import logging
from typing import List, Optional, Union, Dict, Any, Tuple
import hashlib
import warnings

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.SaltRemover import SaltRemover

from ..core.exceptions import CriosError, ChemistryError

# Suppress RDKit warnings for cleaner output
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

logger = logging.getLogger(__name__)


class MoleculeError(CriosError):
    """Molecule-specific errors"""
    pass


class StandardizationError(MoleculeError):
    """Standardization-specific errors"""
    pass


class Molecule:
    """
    Core molecule wrapper with RDKit integration
    Handles SMILES parsing, validation, standardization, and property calculation
    """
    
    def __init__(
        self,
        smiles: str,
        mol_id: Optional[str] = None,
        standardize: bool = True,
        validate: bool = True
    ):
        """
        Initialize molecule from SMILES
        
        Args:
            smiles: SMILES string
            mol_id: Unique molecule identifier  
            standardize: Apply standardization procedures
            validate: Validate molecule structure
        """
        self.original_smiles = smiles.strip()
        self.mol_id = mol_id or self._generate_id(self.original_smiles)
        self._mol = None
        self._standardized_smiles = None
        self._properties = {}
        self._fingerprints = {}
        
        # Parse and optionally standardize
        try:
            self._mol = Chem.MolFromSmiles(self.original_smiles)
            if self._mol is None:
                raise MoleculeError(f"Cannot parse SMILES: {self.original_smiles}")
            
            if standardize:
                self._mol = self.standardize_molecule(self._mol)
                self._standardized_smiles = Chem.MolToSmiles(self._mol)
            else:
                self._standardized_smiles = self.original_smiles
            
            if validate and not self.is_valid():
                raise MoleculeError(f"Invalid molecule structure: {self.original_smiles}")
                
        except Exception as e:
            logger.error(f"Failed to initialize molecule {self.mol_id}: {e}")
            raise MoleculeError(f"Molecule initialization failed: {e}")
    
    @staticmethod
    def _generate_id(smiles: str) -> str:
        """Generate unique ID from SMILES"""
        return hashlib.md5(smiles.encode()).hexdigest()[:12]
    
    @property
    def smiles(self) -> str:
        """Get standardized SMILES"""
        return self._standardized_smiles or self.original_smiles
    
    @property
    def mol(self) -> Optional[Chem.Mol]:
        """Get RDKit Mol object"""
        return self._mol
    
    def is_valid(self) -> bool:
        """Check if molecule is valid"""
        return self._mol is not None
    
    @staticmethod
    def standardize_molecule(mol: Chem.Mol) -> Chem.Mol:
        """
        Standardize molecule using RDKit
        - Neutralize charges
        - Sanitize 
        - Canonical tautomer
        - Strip salts
        """
        if mol is None:
            return None
            
        try:
            # Make a copy to avoid modifying original
            mol = Chem.Mol(mol)
            
            # Sanitize
            Chem.SanitizeMol(mol)
            
            # Remove salts
            salt_remover = SaltRemover()
            mol = salt_remover.StripMol(mol)
            
            # Neutralize charges
            neutralizer = rdMolStandardize.Uncharger()
            mol = neutralizer.uncharge(mol)
            
            # Canonical tautomer
            tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
            mol = tautomer_enumerator.Canonicalize(mol)
            
            return mol
            
        except Exception as e:
            logger.warning(f"Standardization failed: {e}")
            return mol  # Return original if standardization fails
    
    def get_property(self, name: str, calculate_if_missing: bool = True) -> Optional[float]:
        """Get molecular property, calculating if needed"""
        if name in self._properties:
            return self._properties[name]
        
        if not calculate_if_missing or not self.is_valid():
            return None
        
        try:
            value = None
            if name == "molecular_weight":
                value = Descriptors.MolWt(self._mol)
            elif name == "logp":
                value = Descriptors.MolLogP(self._mol)
            elif name == "hbd":
                value = Descriptors.NumHDonors(self._mol)
            elif name == "hba":
                value = Descriptors.NumHAcceptors(self._mol)
            elif name == "tpsa":
                value = Descriptors.TPSA(self._mol)
            elif name == "rotatable_bonds":
                value = Descriptors.NumRotatableBonds(self._mol)
            elif name == "aromatic_rings":
                value = Descriptors.NumAromaticRings(self._mol)
            elif name == "aliphatic_rings":
                value = Descriptors.NumAliphaticRings(self._mol)
            elif name == "heavy_atoms":
                value = self._mol.GetNumHeavyAtoms()
            
            if value is not None:
                self._properties[name] = value
                
            return value
            
        except Exception as e:
            logger.warning(f"Failed to calculate property {name}: {e}")
            return None
    
    def get_all_properties(self) -> Dict[str, float]:
        """Calculate and return all standard properties"""
        properties = [
            "molecular_weight", "logp", "hbd", "hba", "tpsa",
            "rotatable_bonds", "aromatic_rings", "aliphatic_rings", "heavy_atoms"
        ]
        
        result = {}
        for prop in properties:
            value = self.get_property(prop)
            if value is not None:
                result[prop] = value
        
        return result
    
    def passes_drug_like_filters(self, config: Optional[Dict[str, Any]] = None) -> bool:
        """
        Check if molecule passes drug-like filters
        Default: Lipinski Rule of Five + extended criteria
        """
        if not self.is_valid():
            return False
        
        # Default thresholds (Lipinski + extensions)
        defaults = {
            "molecular_weight": {"min": 150.0, "max": 500.0},
            "logp": {"min": -2.0, "max": 5.0},
            "hbd": {"max": 5},
            "hba": {"max": 10},
            "tpsa": {"min": 20.0, "max": 140.0},
            "rotatable_bonds": {"max": 10},
            "aromatic_rings": {"max": 4},
            "heavy_atoms": {"min": 10, "max": 70}
        }
        
        if config:
            defaults.update(config)
        
        # Check each filter
        for prop_name, thresholds in defaults.items():
            value = self.get_property(prop_name)
            if value is None:
                continue
                
            if "min" in thresholds and value < thresholds["min"]:
                return False
            if "max" in thresholds and value > thresholds["max"]:
                return False
        
        return True
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation"""
        data = {
            "mol_id": self.mol_id,
            "original_smiles": self.original_smiles,
            "smiles": self.smiles,
            "is_valid": self.is_valid(),
            "drug_like": self.passes_drug_like_filters() if self.is_valid() else False
        }
        
        # Add properties
        data.update(self.get_all_properties())
        
        return data
    
    def __str__(self) -> str:
        return f"Molecule(id={self.mol_id}, smiles={self.smiles[:50]}...)"
    
    def __repr__(self) -> str:
        return self.__str__()


def standardize_smiles_list(
    smiles_list: List[str],
    sanitize: bool = True,
    neutralize: bool = True,
    remove_salts: bool = True,
    canonical_tautomer: bool = True,
    error_handling: str = "skip"
) -> List[Optional[str]]:
    """
    Vectorized SMILES standardization
    
    Args:
        smiles_list: List of SMILES strings
        sanitize: Apply sanitization
        neutralize: Neutralize charges
        remove_salts: Remove salts
        canonical_tautomer: Use canonical tautomer
        error_handling: How to handle errors ("skip", "raise", "none")
    
    Returns:
        List of standardized SMILES (None for failed)
    """
    results = []
    salt_remover = SaltRemover() if remove_salts else None
    neutralizer = rdMolStandardize.Uncharger() if neutralize else None
    tautomer_enum = rdMolStandardize.TautomerEnumerator() if canonical_tautomer else None
    
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles.strip())
            if mol is None:
                if error_handling == "raise":
                    raise StandardizationError(f"Cannot parse SMILES: {smiles}")
                results.append(None)
                continue
            
            # Apply standardization steps
            if sanitize:
                Chem.SanitizeMol(mol)
            
            if salt_remover:
                mol = salt_remover.StripMol(mol)
            
            if neutralizer:
                mol = neutralizer.uncharge(mol)
            
            if tautomer_enum:
                mol = tautomer_enum.Canonicalize(mol)
            
            # Convert back to SMILES
            standardized = Chem.MolToSmiles(mol)
            results.append(standardized)
            
        except Exception as e:
            if error_handling == "raise":
                raise StandardizationError(f"Standardization failed for {smiles}: {e}")
            elif error_handling == "none":
                results.append(smiles)  # Return original
            else:  # skip
                results.append(None)
    
    return results


def validate_smiles_list(
    smiles_list: List[str],
    return_mols: bool = False
) -> Union[List[bool], Tuple[List[bool], List[Optional[Chem.Mol]]]]:
    """
    Vectorized SMILES validation
    
    Args:
        smiles_list: List of SMILES strings
        return_mols: Also return parsed Mol objects
    
    Returns:
        List of validation results, optionally with Mol objects
    """
    valid_flags = []
    mols = [] if return_mols else None
    
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles.strip())
            is_valid = mol is not None
            valid_flags.append(is_valid)
            
            if return_mols:
                mols.append(mol)
                
        except Exception:
            valid_flags.append(False)
            if return_mols:
                mols.append(None)
    
    if return_mols:
        return valid_flags, mols
    return valid_flags


def calculate_properties_batch(
    smiles_list: List[str],
    properties: List[str] = None,
    standardize: bool = True,
    error_handling: str = "skip"
) -> List[Dict[str, Optional[float]]]:
    """
    Vectorized property calculation
    
    Args:
        smiles_list: List of SMILES strings
        properties: Properties to calculate (None for all standard)
        standardize: Apply standardization
        error_handling: Error handling strategy
    
    Returns:
        List of property dictionaries
    """
    if properties is None:
        properties = [
            "molecular_weight", "logp", "hbd", "hba", "tpsa",
            "rotatable_bonds", "aromatic_rings", "aliphatic_rings", "heavy_atoms"
        ]
    
    results = []
    
    for smiles in smiles_list:
        result = {prop: None for prop in properties}
        
        try:
            mol = Molecule(smiles, standardize=standardize)
            if mol.is_valid():
                for prop in properties:
                    result[prop] = mol.get_property(prop)
            
        except Exception as e:
            if error_handling == "raise":
                raise MoleculeError(f"Property calculation failed for {smiles}: {e}")
            # Otherwise, leave as None values
        
        results.append(result)
    
    return results


def filter_drug_like_batch(
    smiles_list: List[str],
    config: Optional[Dict[str, Any]] = None,
    return_molecules: bool = False
) -> Union[List[bool], Tuple[List[bool], List[Optional[Molecule]]]]:
    """
    Vectorized drug-likeness filtering
    
    Args:
        smiles_list: List of SMILES strings
        config: Filter configuration
        return_molecules: Also return Molecule objects
    
    Returns:
        List of filter results, optionally with Molecule objects
    """
    results = []
    molecules = [] if return_molecules else None
    
    for smiles in smiles_list:
        try:
            mol = Molecule(smiles, standardize=True)
            passes = mol.passes_drug_like_filters(config)
            results.append(passes)
            
            if return_molecules:
                molecules.append(mol)
                
        except Exception:
            results.append(False)
            if return_molecules:
                molecules.append(None)
    
    if return_molecules:
        return results, molecules
    return results


# Utility functions
def smiles_to_inchi(smiles: str) -> Optional[str]:
    """Convert SMILES to InChI"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToInchi(mol)
    except:
        return None


def smiles_to_inchi_key(smiles: str) -> Optional[str]:
    """Convert SMILES to InChI Key"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToInchiKey(mol)
    except:
        return None


def canonical_smiles(smiles: str) -> Optional[str]:
    """Get canonical SMILES representation"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
    except:
        return None