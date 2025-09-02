from typing import Union, Optional, Dict, Any
from pathlib import Path
import logging
from .readers import ReaderFactory
from .writers import WriterFactory

logger = logging.getLogger(__name__)


class FormatConverter:
    
    def __init__(self, validate: bool = True, canonical: bool = True,
                 remove_duplicates: bool = False):
        self.validate = validate
        self.canonical = canonical
        self.remove_duplicates = remove_duplicates
    
    def convert(self, input_path: Union[str, Path], output_path: Union[str, Path],
                input_format: Optional[str] = None, output_format: Optional[str] = None,
                batch_size: int = 1000, show_progress: bool = True) -> Dict[str, Any]:
        
        input_path = Path(input_path)
        output_path = Path(output_path)
        
        reader = ReaderFactory.create_reader(input_path, format=input_format, 
                                            validate=self.validate)
        writer = WriterFactory.create_writer(output_path, format=output_format,
                                            canonical_smiles=self.canonical)
        
        processed = 0
        written = 0
        errors = 0
        seen_smiles = set() if self.remove_duplicates else None
        
        batch = []
        
        molecules = reader.read()
        
        for mol_data in molecules:
            try:
                if self.remove_duplicates:
                    smiles = mol_data.get("smiles", "")
                    if smiles in seen_smiles:
                        continue
                    seen_smiles.add(smiles)
                
                batch.append(mol_data)
                processed += 1
                
                if len(batch) >= batch_size:
                    written += writer.write(batch)
                    batch = []
                    
            except Exception as e:
                logger.error(f"Error processing molecule: {e}")
                errors += 1
        
        if batch:
            written += writer.write(batch)
        
        return {
            "input_file": str(input_path),
            "output_file": str(output_path),
            "molecules_processed": processed,
            "molecules_written": written,
            "errors": errors,
            "duplicates_removed": processed - written if self.remove_duplicates else 0
        }
    
    @staticmethod
    def convert_file(input_path: Union[str, Path], output_path: Union[str, Path],
                    **kwargs) -> Dict[str, Any]:
        
        converter = FormatConverter(**kwargs)
        return converter.convert(input_path, output_path)