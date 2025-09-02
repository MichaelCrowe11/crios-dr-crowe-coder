"""
CriOS Core Molecule Classes
RDKit-based molecular representation and processing
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union, Any, Iterator
from pathlib import Path
import pickle

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator

from ..config.schemas import SafetyCategory
from ..exceptions import MoleculeError, ValidationError

logger = logging.getLogger(__name__)


@dataclass
class MolecularProperties:
    """Comprehensive molecular property container"""
    
    # Basic properties
    molecular_weight: Optional[float] = None
    heavy_atom_count: Optional[int] = None
    molecular_formula: Optional[str] = None
    
    # Lipinski Rule of Five
    logp: Optional[float] = None
    hbd: Optional[int] = None  # Hydrogen bond donors
    hba: Optional[int] = None  # Hydrogen bond acceptors
    tpsa: Optional[float] = None  # Topological polar surface area
    
    # Additional drug-like properties
    rotatable_bonds: Optional[int] = None
    aromatic_rings: Optional[int] = None
    aliphatic_rings: Optional[int] = None
    
    # Complexity metrics
    bertz_complexity: Optional[float] = None
    mol_mr: Optional[float] = None  # Molar refractivity
    
    # Custom properties
    custom_properties: Dict[str, Any] = field(default_factory=dict)


class Molecule:
    """
    Core molecular representation with RDKit integration
    Handles SMILES, InChI, SDF parsing and molecular property calculation
    """
    
    def __init__(
        self,
        identifier: Union[str, Chem.Mol],
        mol_id: Optional[str] = None,
        name: Optional[str] = None,
        source: Optional[str] = None,
        validate: bool = True
    ):
        """
        Initialize molecule from SMILES, InChI string, or RDKit Mol object
        
        Args:
            identifier: SMILES string, InChI string, or RDKit Mol object
            mol_id: Unique identifier for the molecule
            name: Human-readable name
            source: Source database or origin
            validate: Whether to validate the molecule structure
        """
        self.mol_id = mol_id or f"mol_{id(self)}"
        self.name = name
        self.source = source
        self.properties = MolecularProperties()
        
        # Initialize RDKit molecule
        if isinstance(identifier, str):
            self._mol = self._parse_string(identifier)
            self._smiles = identifier if self._is_smiles(identifier) else None
        elif isinstance(identifier, Chem.Mol):
            self._mol = identifier
            self._smiles = None
        else:
            raise MoleculeError(f"Unsupported identifier type: {type(identifier)}")
        
        # Validation
        if validate and not self.is_valid():
            raise ValidationError(f"Invalid molecule: {identifier}")
        
        # Calculate basic properties on initialization
        self._calculate_basic_properties()
    
    def _parse_string(self, identifier: str) -> Chem.Mol:
        """Parse string identifier to RDKit Mol object"""
        try:
            if self._is_smiles(identifier):
                mol = Chem.MolFromSmiles(identifier)
            elif identifier.startswith("InChI"):
                mol = Chem.MolFromInchi(identifier)
            else:
                # Try SMILES first, then InChI
                mol = Chem.MolFromSmiles(identifier)
                if mol is None:
                    mol = Chem.MolFromInchi(identifier)
            
            if mol is None:
                raise MoleculeError(f"Could not parse identifier: {identifier}")
            
            return mol
            
        except Exception as e:
            raise MoleculeError(f"Failed to parse molecule: {e}")
    
    def _is_smiles(self, identifier: str) -> bool:
        """Heuristic to determine if string is SMILES format"""
        return not identifier.startswith("InChI")
    
    def _calculate_basic_properties(self) -> None:
        """Calculate basic molecular properties"""
        if not self.is_valid():
            return
        
        try:
            mol = self._mol
            
            # Basic properties
            self.properties.molecular_weight = Descriptors.MolWt(mol)
            self.properties.heavy_atom_count = mol.GetNumHeavyAtoms()
            self.properties.molecular_formula = CalcMolFormula(mol)
            
            # Lipinski properties
            self.properties.logp = Crippen.MolLogP(mol)
            self.properties.hbd = Lipinski.NumHDonors(mol)
            self.properties.hba = Lipinski.NumHAcceptors(mol)
            self.properties.tpsa = Descriptors.TPSA(mol)
            
            # Additional properties
            self.properties.rotatable_bonds = Descriptors.NumRotatableBonds(mol)
            self.properties.aromatic_rings = Descriptors.NumAromaticRings(mol)
            self.properties.aliphatic_rings = Descriptors.NumAliphaticRings(mol)
            
            # Complexity
            self.properties.bertz_complexity = Descriptors.BertzCT(mol)
            self.properties.mol_mr = Crippen.MolMR(mol)
            
        except Exception as e:
            logger.warning(f"Failed to calculate properties for {self.mol_id}: {e}")
    
    @property
    def rdkit_mol(self) -> Optional[Chem.Mol]:
        """Get RDKit Mol object"""
        return self._mol
    
    @property
    def smiles(self) -> Optional[str]:
        """Get canonical SMILES representation"""
        if self._smiles is None and self.is_valid():
            self._smiles = Chem.MolToSmiles(self._mol)
        return self._smiles
    
    @property
    def inchi(self) -> Optional[str]:
        """Get InChI representation"""
        if self.is_valid():
            return Chem.MolToInchi(self._mol)
        return None
    
    @property
    def inchi_key(self) -> Optional[str]:
        """Get InChI key"""
        if self.is_valid():
            return Chem.MolToInchiKey(self._mol)
        return None
    
    def is_valid(self) -> bool:
        """Check if molecule is valid"""
        return self._mol is not None
    
    def get_fingerprint(
        self,
        fp_type: str = "morgan",
        radius: int = 2,
        n_bits: int = 2048
    ) -> Optional[np.ndarray]:
        """
        Generate molecular fingerprint
        
        Args:
            fp_type: Fingerprint type ('morgan', 'rdkit', 'maccs')
            radius: Morgan fingerprint radius
            n_bits: Number of bits in fingerprint
        """
        if not self.is_valid():
            return None
        
        try:
            if fp_type.lower() == "morgan":
                fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
                    self._mol, radius, nBits=n_bits
                )
            elif fp_type.lower() == "rdkit":
                fp = Chem.RDKFingerprint(self._mol, fpSize=n_bits)
            elif fp_type.lower() == "maccs":
                fp = rdMolDescriptors.GetMACCSKeysFingerprint(self._mol)
            else:
                raise ValueError(f"Unsupported fingerprint type: {fp_type}")
            
            return np.array(fp)
            
        except Exception as e:
            logger.error(f"Failed to generate fingerprint for {self.mol_id}: {e}")
            return None
    
    def passes_lipinski(self) -> bool:
        """Check if molecule passes Lipinski's Rule of Five"""
        if not self.is_valid():
            return False
        
        props = self.properties
        return (
            props.molecular_weight <= 500 and
            props.logp <= 5 and
            props.hbd <= 5 and
            props.hba <= 10
        )
    
    def passes_veber(self) -> bool:
        """Check if molecule passes Veber's rules"""
        if not self.is_valid():
            return False
        
        props = self.properties
        return (
            props.rotatable_bonds <= 10 and
            props.tpsa <= 140
        )
    
    def is_drug_like(self) -> bool:
        """Check if molecule is drug-like (passes both Lipinski and Veber)"""
        return self.passes_lipinski() and self.passes_veber()
    
    def add_property(self, name: str, value: Any) -> None:
        """Add custom property to molecule"""
        self.properties.custom_properties[name] = value
    
    def get_property(self, name: str) -> Any:
        """Get custom property from molecule"""
        return self.properties.custom_properties.get(name)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert molecule to dictionary representation"""
        return {
            "mol_id": self.mol_id,
            "name": self.name,
            "source": self.source,
            "smiles": self.smiles,
            "inchi": self.inchi,
            "inchi_key": self.inchi_key,
            "molecular_weight": self.properties.molecular_weight,
            "heavy_atom_count": self.properties.heavy_atom_count,
            "molecular_formula": self.properties.molecular_formula,
            "logp": self.properties.logp,
            "hbd": self.properties.hbd,
            "hba": self.properties.hba,
            "tpsa": self.properties.tpsa,
            "rotatable_bonds": self.properties.rotatable_bonds,
            "aromatic_rings": self.properties.aromatic_rings,
            "passes_lipinski": self.passes_lipinski(),
            "passes_veber": self.passes_veber(),
            "is_drug_like": self.is_drug_like(),
            "custom_properties": self.properties.custom_properties
        }
    
    def __str__(self) -> str:
        return f"Molecule(id={self.mol_id}, name={self.name}, smiles={self.smiles})"
    
    def __repr__(self) -> str:
        return self.__str__()


class MoleculeCollection:
    """
    Collection of molecules with batch processing capabilities
    Optimized for high-throughput compound processing
    """
    
    def __init__(self, molecules: Optional[List[Molecule]] = None):
        """
        Initialize molecule collection
        
        Args:
            molecules: List of Molecule objects
        """
        self.molecules: Dict[str, Molecule] = {}
        self._index = 0
        
        if molecules:
            for mol in molecules:
                self.add_molecule(mol)
    
    def add_molecule(self, molecule: Molecule) -> None:
        """Add molecule to collection"""
        self.molecules[molecule.mol_id] = molecule
    
    def remove_molecule(self, mol_id: str) -> bool:
        """Remove molecule from collection"""
        if mol_id in self.molecules:
            del self.molecules[mol_id]
            return True
        return False
    
    def get_molecule(self, mol_id: str) -> Optional[Molecule]:
        """Get molecule by ID"""
        return self.molecules.get(mol_id)
    
    def filter_valid(self) -> 'MoleculeCollection':
        """Return collection containing only valid molecules"""
        valid_mols = [mol for mol in self.molecules.values() if mol.is_valid()]
        return MoleculeCollection(valid_mols)
    
    def filter_drug_like(self) -> 'MoleculeCollection':
        """Return collection containing only drug-like molecules"""
        drug_like_mols = [mol for mol in self.molecules.values() if mol.is_drug_like()]
        return MoleculeCollection(drug_like_mols)
    
    def get_fingerprint_matrix(
        self,
        fp_type: str = "morgan",
        radius: int = 2,
        n_bits: int = 2048
    ) -> Optional[np.ndarray]:
        """
        Get fingerprint matrix for all molecules
        
        Returns:
            Matrix where each row is a molecule fingerprint
        """
        valid_mols = [mol for mol in self.molecules.values() if mol.is_valid()]
        if not valid_mols:
            return None
        
        fps = []
        for mol in valid_mols:
            fp = mol.get_fingerprint(fp_type, radius, n_bits)
            if fp is not None:
                fps.append(fp)
        
        if not fps:
            return None
        
        return np.array(fps)
    
    def to_dataframe(self) -> pd.DataFrame:
        """Convert collection to pandas DataFrame"""
        data = []
        for mol in self.molecules.values():
            data.append(mol.to_dict())
        
        return pd.DataFrame(data)
    
    def to_sdf(self, filename: Union[str, Path]) -> None:
        """Export collection to SDF file"""
        filename = Path(filename)
        
        with Chem.SDWriter(str(filename)) as writer:
            for mol in self.molecules.values():
                if mol.is_valid():
                    # Add properties to molecule
                    rdkit_mol = mol.rdkit_mol
                    for prop_name, prop_value in mol.to_dict().items():
                        if prop_value is not None and prop_name != "rdkit_mol":
                            rdkit_mol.SetProp(prop_name, str(prop_value))
                    
                    writer.write(rdkit_mol)
    
    def save(self, filename: Union[str, Path]) -> None:
        """Save collection to pickle file"""
        filename = Path(filename)
        with open(filename, 'wb') as f:
            pickle.dump(self.molecules, f)
    
    @classmethod
    def load(cls, filename: Union[str, Path]) -> 'MoleculeCollection':
        """Load collection from pickle file"""
        filename = Path(filename)
        with open(filename, 'rb') as f:
            molecules_dict = pickle.load(f)
        
        collection = cls()
        collection.molecules = molecules_dict
        return collection
    
    @classmethod
    def from_smiles(
        cls,
        smiles_list: List[str],
        mol_ids: Optional[List[str]] = None,
        names: Optional[List[str]] = None
    ) -> 'MoleculeCollection':
        """Create collection from list of SMILES strings"""
        collection = cls()
        
        for i, smiles in enumerate(smiles_list):
            mol_id = mol_ids[i] if mol_ids else f"mol_{i}"
            name = names[i] if names else None
            
            try:
                mol = Molecule(smiles, mol_id=mol_id, name=name)
                collection.add_molecule(mol)
            except Exception as e:
                logger.warning(f"Failed to create molecule from SMILES {smiles}: {e}")
        
        return collection
    
    @classmethod
    def from_sdf(cls, filename: Union[str, Path]) -> 'MoleculeCollection':
        """Load collection from SDF file"""
        filename = Path(filename)
        collection = cls()
        
        supplier = Chem.SDMolSupplier(str(filename))
        for i, mol in enumerate(supplier):
            if mol is not None:
                mol_id = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{i}"
                molecule = Molecule(mol, mol_id=mol_id)
                
                # Extract properties from SDF
                for prop_name in mol.GetPropNames():
                    prop_value = mol.GetProp(prop_name)
                    molecule.add_property(prop_name, prop_value)
                
                collection.add_molecule(molecule)
        
        return collection
    
    def __len__(self) -> int:
        return len(self.molecules)
    
    def __iter__(self) -> Iterator[Molecule]:
        return iter(self.molecules.values())
    
    def __getitem__(self, mol_id: str) -> Molecule:
        return self.molecules[mol_id]
    
    def __contains__(self, mol_id: str) -> bool:
        return mol_id in self.molecules
    
    def __str__(self) -> str:
        return f"MoleculeCollection({len(self.molecules)} molecules)"
    
    def __repr__(self) -> str:
        return self.__str__()


# Utility functions
def validate_smiles(smiles: str) -> bool:
    """Validate SMILES string using RDKit"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False


def standardize_smiles(smiles: str) -> Optional[str]:
    """Standardize SMILES string to canonical form"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
    except:
        return None


def smiles_to_inchi(smiles: str) -> Optional[str]:
    """Convert SMILES to InChI"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToInchi(mol)
    except:
        return None