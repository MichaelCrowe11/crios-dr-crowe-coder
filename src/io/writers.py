from rdkit import Chem
from typing import List, Dict, Any, Union, Optional, TextIO
import pandas as pd
import csv
from pathlib import Path
import logging
from abc import ABC, abstractmethod

logger = logging.getLogger(__name__)


class MoleculeWriter(ABC):
    
    def __init__(self, file_path: Union[str, Path], overwrite: bool = True):
        self.file_path = Path(file_path)
        self.overwrite = overwrite
        self._molecule_count = 0
        
        if self.file_path.exists() and not self.overwrite:
            raise FileExistsError(f"File already exists: {self.file_path}")
        
        self.file_path.parent.mkdir(parents=True, exist_ok=True)
    
    @abstractmethod
    def write(self, molecules: List[Dict[str, Any]]) -> int:
        pass
    
    @abstractmethod
    def write_single(self, molecule: Dict[str, Any]) -> bool:
        pass
    
    def get_count(self) -> int:
        return self._molecule_count


class SMILESWriter(MoleculeWriter):
    
    def __init__(self, file_path: Union[str, Path], overwrite: bool = True,
                 delimiter: str = '\t', include_header: bool = False,
                 canonical: bool = True):
        super().__init__(file_path, overwrite)
        self.delimiter = delimiter
        self.include_header = include_header
        self.canonical = canonical
        self._file_handle: Optional[TextIO] = None
        self._header_written = False
    
    def _open_file(self):
        if self._file_handle is None:
            self._file_handle = open(self.file_path, 'w', encoding='utf-8')
            if self.include_header and not self._header_written:
                self._file_handle.write(f"SMILES{self.delimiter}ID\n")
                self._header_written = True
    
    def _close_file(self):
        if self._file_handle:
            self._file_handle.close()
            self._file_handle = None
    
    def write_single(self, molecule: Dict[str, Any]) -> bool:
        self._open_file()
        
        try:
            if "smiles" in molecule:
                smiles = molecule["smiles"]
                if self.canonical and "mol" in molecule and molecule["mol"]:
                    smiles = Chem.MolToSmiles(molecule["mol"], canonical=True)
            elif "mol" in molecule and molecule["mol"]:
                smiles = Chem.MolToSmiles(molecule["mol"], canonical=self.canonical)
            else:
                logger.warning(f"No SMILES or mol object for molecule {molecule.get('mol_id', 'unknown')}")
                return False
            
            mol_id = molecule.get("mol_id", f"MOL_{self._molecule_count:06d}")
            
            self._file_handle.write(f"{smiles}{self.delimiter}{mol_id}\n")
            self._molecule_count += 1
            return True
            
        except Exception as e:
            logger.error(f"Error writing molecule: {e}")
            return False
    
    def write(self, molecules: List[Dict[str, Any]]) -> int:
        self._open_file()
        count = 0
        
        for mol in molecules:
            if self.write_single(mol):
                count += 1
        
        self._close_file()
        return count
    
    def __del__(self):
        self._close_file()


class SDFWriter(MoleculeWriter):
    
    def __init__(self, file_path: Union[str, Path], overwrite: bool = True,
                 kekulize: bool = False):
        super().__init__(file_path, overwrite)
        self.kekulize = kekulize
        self._writer: Optional[Chem.SDWriter] = None
    
    def _open_writer(self):
        if self._writer is None:
            self._writer = Chem.SDWriter(str(self.file_path))
            self._writer.SetKekulize(self.kekulize)
    
    def _close_writer(self):
        if self._writer:
            self._writer.close()
            self._writer = None
    
    def write_single(self, molecule: Dict[str, Any]) -> bool:
        self._open_writer()
        
        try:
            if "mol" in molecule and molecule["mol"]:
                mol = molecule["mol"]
            elif "smiles" in molecule:
                mol = Chem.MolFromSmiles(molecule["smiles"])
                if mol is None:
                    logger.warning(f"Invalid SMILES: {molecule['smiles'][:50]}")
                    return False
            else:
                logger.warning("No mol object or SMILES in molecule data")
                return False
            
            if "mol_id" in molecule:
                mol.SetProp("_Name", str(molecule["mol_id"]))
            
            if "properties" in molecule:
                for prop_name, prop_value in molecule["properties"].items():
                    mol.SetProp(prop_name, str(prop_value))
            
            self._writer.write(mol)
            self._molecule_count += 1
            return True
            
        except Exception as e:
            logger.error(f"Error writing molecule: {e}")
            return False
    
    def write(self, molecules: List[Dict[str, Any]]) -> int:
        self._open_writer()
        count = 0
        
        for mol in molecules:
            if self.write_single(mol):
                count += 1
        
        self._close_writer()
        return count
    
    def __del__(self):
        self._close_writer()


class CSVWriter(MoleculeWriter):
    
    def __init__(self, file_path: Union[str, Path], overwrite: bool = True,
                 canonical_smiles: bool = True, include_descriptors: bool = True):
        super().__init__(file_path, overwrite)
        self.canonical_smiles = canonical_smiles
        self.include_descriptors = include_descriptors
        self._data_buffer = []
    
    def write_single(self, molecule: Dict[str, Any]) -> bool:
        try:
            row_data = {}
            
            row_data["mol_id"] = molecule.get("mol_id", f"MOL_{self._molecule_count:06d}")
            
            if "smiles" in molecule:
                smiles = molecule["smiles"]
                if self.canonical_smiles and "mol" in molecule and molecule["mol"]:
                    smiles = Chem.MolToSmiles(molecule["mol"], canonical=True)
            elif "mol" in molecule and molecule["mol"]:
                smiles = Chem.MolToSmiles(molecule["mol"], canonical=self.canonical_smiles)
            else:
                logger.warning(f"No SMILES or mol object for molecule {row_data['mol_id']}")
                return False
            
            row_data["smiles"] = smiles
            
            if self.include_descriptors and "descriptors" in molecule:
                row_data.update(molecule["descriptors"])
            
            if "properties" in molecule:
                row_data.update(molecule["properties"])
            
            self._data_buffer.append(row_data)
            self._molecule_count += 1
            return True
            
        except Exception as e:
            logger.error(f"Error processing molecule: {e}")
            return False
    
    def write(self, molecules: List[Dict[str, Any]]) -> int:
        self._data_buffer = []
        count = 0
        
        for mol in molecules:
            if self.write_single(mol):
                count += 1
        
        if self._data_buffer:
            df = pd.DataFrame(self._data_buffer)
            
            columns_order = ["mol_id", "smiles"]
            other_cols = [col for col in df.columns if col not in columns_order]
            df = df[columns_order + other_cols]
            
            df.to_csv(self.file_path, index=False)
        
        return count


class WriterFactory:
    
    @staticmethod
    def create_writer(file_path: Union[str, Path], format: Optional[str] = None,
                     **kwargs) -> MoleculeWriter:
        
        file_path = Path(file_path)
        
        if format is None:
            suffix = file_path.suffix.lower()
            format_map = {
                '.smi': 'smiles',
                '.smiles': 'smiles',
                '.sdf': 'sdf',
                '.mol': 'sdf',
                '.csv': 'csv',
                '.tsv': 'csv'
            }
            format = format_map.get(suffix)
            
            if format is None:
                raise ValueError(f"Cannot determine format from extension: {suffix}")
        
        format = format.lower()
        
        if format in ['smiles', 'smi']:
            return SMILESWriter(file_path, **kwargs)
        elif format in ['sdf', 'mol']:
            return SDFWriter(file_path, **kwargs)
        elif format in ['csv', 'tsv']:
            return CSVWriter(file_path, **kwargs)
        else:
            raise ValueError(f"Unsupported format: {format}")