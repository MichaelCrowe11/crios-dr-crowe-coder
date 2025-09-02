from rdkit import Chem
from typing import Iterator, Optional, List, Dict, Any, Union
import pandas as pd
import csv
from pathlib import Path
import logging
from abc import ABC, abstractmethod

logger = logging.getLogger(__name__)


class MoleculeReader(ABC):
    
    def __init__(self, file_path: Union[str, Path], validate: bool = True):
        self.file_path = Path(file_path)
        self.validate = validate
        self._molecule_count = 0
        
        if not self.file_path.exists():
            raise FileNotFoundError(f"File not found: {self.file_path}")
    
    @abstractmethod
    def read(self) -> Iterator[Dict[str, Any]]:
        pass
    
    def read_all(self) -> List[Dict[str, Any]]:
        return list(self.read())
    
    def get_count(self) -> int:
        return self._molecule_count


class SMILESReader(MoleculeReader):
    
    def __init__(self, file_path: Union[str, Path], validate: bool = True,
                 delimiter: str = '\t', has_header: bool = False,
                 smiles_column: int = 0, id_column: Optional[int] = 1):
        super().__init__(file_path, validate)
        self.delimiter = delimiter
        self.has_header = has_header
        self.smiles_column = smiles_column
        self.id_column = id_column
    
    def read(self) -> Iterator[Dict[str, Any]]:
        with open(self.file_path, 'r', encoding='utf-8') as f:
            if self.has_header:
                header = next(f).strip().split(self.delimiter)
            
            for line_num, line in enumerate(f, start=1 if not self.has_header else 2):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split(self.delimiter)
                
                if len(parts) <= self.smiles_column:
                    logger.warning(f"Line {line_num}: Not enough columns")
                    continue
                
                smiles = parts[self.smiles_column]
                mol_id = None
                
                if self.id_column is not None and len(parts) > self.id_column:
                    mol_id = parts[self.id_column]
                else:
                    mol_id = f"MOL_{line_num:06d}"
                
                if self.validate:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is None:
                        logger.warning(f"Line {line_num}: Invalid SMILES: {smiles[:50]}")
                        continue
                else:
                    mol = None
                
                self._molecule_count += 1
                
                yield {
                    "smiles": smiles,
                    "mol_id": mol_id,
                    "mol": mol,
                    "source_line": line_num,
                    "properties": {}
                }


class SDFReader(MoleculeReader):
    
    def __init__(self, file_path: Union[str, Path], validate: bool = True,
                 sanitize: bool = True, remove_hs: bool = False):
        super().__init__(file_path, validate)
        self.sanitize = sanitize
        self.remove_hs = remove_hs
    
    def read(self) -> Iterator[Dict[str, Any]]:
        supplier = Chem.SDMolSupplier(str(self.file_path), 
                                      sanitize=self.sanitize,
                                      removeHs=self.remove_hs)
        
        for idx, mol in enumerate(supplier):
            if mol is None:
                logger.warning(f"Molecule {idx}: Failed to parse")
                if self.validate:
                    continue
            
            mol_id = mol.GetProp("_Name") if mol and mol.HasProp("_Name") else f"MOL_{idx:06d}"
            
            properties = {}
            if mol:
                for prop_name in mol.GetPropNames():
                    if not prop_name.startswith("_"):
                        try:
                            properties[prop_name] = mol.GetProp(prop_name)
                        except:
                            pass
                
                smiles = Chem.MolToSmiles(mol)
            else:
                smiles = None
            
            self._molecule_count += 1
            
            yield {
                "smiles": smiles,
                "mol_id": mol_id,
                "mol": mol,
                "source_index": idx,
                "properties": properties
            }


class CSVReader(MoleculeReader):
    
    def __init__(self, file_path: Union[str, Path], validate: bool = True,
                 smiles_column: str = "smiles", id_column: Optional[str] = "mol_id",
                 property_columns: Optional[List[str]] = None):
        super().__init__(file_path, validate)
        self.smiles_column = smiles_column
        self.id_column = id_column
        self.property_columns = property_columns
    
    def read(self) -> Iterator[Dict[str, Any]]:
        df = pd.read_csv(self.file_path)
        
        if self.smiles_column not in df.columns:
            available_cols = ', '.join(df.columns[:5])
            raise ValueError(f"SMILES column '{self.smiles_column}' not found. Available: {available_cols}")
        
        for idx, row in df.iterrows():
            smiles = row[self.smiles_column]
            
            if pd.isna(smiles):
                logger.warning(f"Row {idx}: Empty SMILES")
                if self.validate:
                    continue
                smiles = ""
            
            if self.id_column and self.id_column in df.columns:
                mol_id = str(row[self.id_column])
            else:
                mol_id = f"MOL_{idx:06d}"
            
            if self.validate:
                mol = Chem.MolFromSmiles(str(smiles))
                if mol is None:
                    logger.warning(f"Row {idx}: Invalid SMILES: {str(smiles)[:50]}")
                    continue
            else:
                mol = None
            
            properties = {}
            if self.property_columns:
                for col in self.property_columns:
                    if col in df.columns:
                        properties[col] = row[col]
            else:
                for col in df.columns:
                    if col not in [self.smiles_column, self.id_column]:
                        properties[col] = row[col]
            
            self._molecule_count += 1
            
            yield {
                "smiles": str(smiles),
                "mol_id": mol_id,
                "mol": mol,
                "source_index": idx,
                "properties": properties
            }


class ReaderFactory:
    
    @staticmethod
    def create_reader(file_path: Union[str, Path], format: Optional[str] = None, 
                     **kwargs) -> MoleculeReader:
        
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
            return SMILESReader(file_path, **kwargs)
        elif format in ['sdf', 'mol']:
            return SDFReader(file_path, **kwargs)
        elif format in ['csv', 'tsv']:
            if file_path.suffix.lower() == '.tsv':
                kwargs.setdefault('delimiter', '\t')
            return CSVReader(file_path, **kwargs)
        else:
            raise ValueError(f"Unsupported format: {format}")