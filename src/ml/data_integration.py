"""
Data Integration Module for CriOS Nova Cheminformatics
Integrates ChEMBL, PubMed, and local drug discovery data
"""

import os
import json
import gzip
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
from datetime import datetime
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests
from Bio import Entrez
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import PandasTools
import sqlite3

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class IntegratedCompound:
    """Unified compound data structure"""
    smiles: str
    inchi: Optional[str]
    inchi_key: Optional[str]
    chembl_id: Optional[str]
    pubchem_cid: Optional[str]
    drugbank_id: Optional[str]
    name: Optional[str]
    synonyms: List[str]
    targets: List[Dict[str, Any]]
    bioactivities: List[Dict[str, Any]]
    publications: List[Dict[str, Any]]
    clinical_trials: List[Dict[str, Any]]
    properties: Dict[str, Any]
    source_databases: List[str]
    last_updated: datetime


class ChEMBLIntegrator:
    """Integrate ChEMBL data from local files and API"""

    def __init__(self, local_data_path: str = None):
        self.local_data_path = Path(local_data_path) if local_data_path else None
        self.molecule = new_client.molecule
        self.target = new_client.target
        self.activity = new_client.activity
        self.assay = new_client.assay
        self.document = new_client.document

    def load_local_targets(self, file_path: str) -> pd.DataFrame:
        """Load ChEMBL targets from local JSONL file"""
        targets = []
        file_path = Path(file_path)

        if file_path.suffix == '.gz':
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'

        with open_func(file_path, mode) as f:
            for line in f:
                try:
                    target_data = json.loads(line)
                    targets.append(self._parse_target(target_data))
                except Exception as e:
                    logger.error(f"Error parsing target: {e}")
                    continue

        return pd.DataFrame(targets)

    def _parse_target(self, target_data: Dict) -> Dict:
        """Parse ChEMBL target data"""
        return {
            'chembl_id': target_data.get('target_chembl_id'),
            'pref_name': target_data.get('pref_name'),
            'organism': target_data.get('organism'),
            'target_type': target_data.get('target_type'),
            'tax_id': target_data.get('tax_id'),
            'components': target_data.get('target_components', [])
        }

    def fetch_compound_activities(self, chembl_id: str) -> List[Dict]:
        """Fetch compound activities from ChEMBL API"""
        try:
            activities = self.activity.filter(
                molecule_chembl_id=chembl_id,
                assay_type='B'  # Binding assays
            ).only(['activity_id', 'assay_chembl_id', 'target_chembl_id',
                   'standard_type', 'standard_value', 'standard_units'])

            return list(activities)
        except Exception as e:
            logger.error(f"Error fetching activities for {chembl_id}: {e}")
            return []

    def fetch_molecule_properties(self, chembl_id: str) -> Dict:
        """Fetch molecular properties from ChEMBL"""
        try:
            mol = self.molecule.get(chembl_id)
            if mol:
                return {
                    'smiles': mol.get('molecule_structures', {}).get('canonical_smiles'),
                    'inchi': mol.get('molecule_structures', {}).get('standard_inchi'),
                    'inchi_key': mol.get('molecule_structures', {}).get('standard_inchi_key'),
                    'molecular_weight': mol.get('molecule_properties', {}).get('mw_freebase'),
                    'logp': mol.get('molecule_properties', {}).get('alogp'),
                    'psa': mol.get('molecule_properties', {}).get('psa'),
                    'num_ro5_violations': mol.get('molecule_properties', {}).get('num_ro5_violations'),
                    'max_phase': mol.get('max_phase'),
                    'therapeutic_flags': mol.get('therapeutic_flags', []),
                    'molecule_type': mol.get('molecule_type')
                }
        except Exception as e:
            logger.error(f"Error fetching molecule {chembl_id}: {e}")
        return {}


class PubMedIntegrator:
    """Integrate PubMed literature data"""

    def __init__(self, email: str = "crios.nova@research.com"):
        Entrez.email = email
        self.search_cache = {}

    def search_compound_literature(self, compound_name: str,
                                  max_results: int = 100) -> List[Dict]:
        """Search PubMed for compound-related literature"""
        if compound_name in self.search_cache:
            return self.search_cache[compound_name]

        try:
            # Search PubMed
            search_query = f"{compound_name}[Title/Abstract] AND drug[Title/Abstract]"
            search_handle = Entrez.esearch(
                db="pubmed",
                term=search_query,
                retmax=max_results,
                sort="relevance"
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()

            # Fetch article details
            id_list = search_results["IdList"]
            if not id_list:
                return []

            fetch_handle = Entrez.efetch(
                db="pubmed",
                id=",".join(id_list),
                rettype="medline",
                retmode="xml"
            )
            articles = Entrez.read(fetch_handle)
            fetch_handle.close()

            # Parse articles
            parsed_articles = []
            for article in articles['PubmedArticle']:
                parsed_articles.append(self._parse_article(article))

            self.search_cache[compound_name] = parsed_articles
            return parsed_articles

        except Exception as e:
            logger.error(f"Error searching PubMed for {compound_name}: {e}")
            return []

    def _parse_article(self, article: Dict) -> Dict:
        """Parse PubMed article data"""
        medline = article.get('MedlineCitation', {})
        article_data = medline.get('Article', {})

        return {
            'pmid': str(medline.get('PMID', '')),
            'title': article_data.get('ArticleTitle', ''),
            'abstract': self._get_abstract(article_data),
            'authors': self._get_authors(article_data),
            'journal': article_data.get('Journal', {}).get('Title', ''),
            'pub_date': self._get_pub_date(article_data),
            'keywords': self._get_keywords(medline),
            'mesh_terms': self._get_mesh_terms(medline)
        }

    def _get_abstract(self, article_data: Dict) -> str:
        """Extract abstract text"""
        abstract = article_data.get('Abstract', {})
        if isinstance(abstract, dict):
            text_parts = abstract.get('AbstractText', [])
            if isinstance(text_parts, list):
                return ' '.join([str(part) for part in text_parts])
            return str(text_parts)
        return ''

    def _get_authors(self, article_data: Dict) -> List[str]:
        """Extract author names"""
        author_list = article_data.get('AuthorList', [])
        authors = []
        for author in author_list:
            if isinstance(author, dict):
                last_name = author.get('LastName', '')
                initials = author.get('Initials', '')
                if last_name:
                    authors.append(f"{last_name} {initials}".strip())
        return authors

    def _get_pub_date(self, article_data: Dict) -> str:
        """Extract publication date"""
        pub_date = article_data.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
        year = pub_date.get('Year', '')
        month = pub_date.get('Month', '')
        day = pub_date.get('Day', '')
        return f"{year}-{month}-{day}".strip('-')

    def _get_keywords(self, medline: Dict) -> List[str]:
        """Extract keywords"""
        keyword_list = medline.get('KeywordList', [])
        keywords = []
        for kw_set in keyword_list:
            if isinstance(kw_set, list):
                keywords.extend([str(kw) for kw in kw_set])
        return keywords

    def _get_mesh_terms(self, medline: Dict) -> List[str]:
        """Extract MeSH terms"""
        mesh_list = medline.get('MeshHeadingList', [])
        mesh_terms = []
        for mesh in mesh_list:
            if isinstance(mesh, dict):
                descriptor = mesh.get('DescriptorName', {})
                if isinstance(descriptor, dict):
                    mesh_terms.append(str(descriptor))
        return mesh_terms


class LocalDataIntegrator:
    """Integrate local drug discovery data files"""

    def __init__(self, base_path: str = "C:/Users/micha"):
        self.base_path = Path(base_path)
        self.data_cache = {}

    def load_synapse_examples(self) -> List[Dict]:
        """Load drug discovery examples from Synapse language files"""
        synapse_files = [
            "synapse-lang/examples/drug_discovery.syn",
            "synapse-lang/examples/drug_discovery_comprehensive.syn",
            "synapse-lang/examples/quantum_drug_discovery.syn",
            "synapse-lang/examples/crios_nova_drug_discovery.syn"
        ]

        examples = []
        for file_path in synapse_files:
            full_path = self.base_path / file_path
            if full_path.exists():
                try:
                    with open(full_path, 'r') as f:
                        content = f.read()
                        examples.append({
                            'file': file_path,
                            'content': content,
                            'type': 'synapse_example'
                        })
                except Exception as e:
                    logger.error(f"Error loading {file_path}: {e}")

        return examples

    def load_ai_workflow(self) -> Dict:
        """Load AI drug discovery workflow"""
        workflow_file = self.base_path / "Downloads/ai_drug_discovery_workflow.py"
        if workflow_file.exists():
            try:
                with open(workflow_file, 'r') as f:
                    return {
                        'content': f.read(),
                        'type': 'python_workflow',
                        'targets': ['CB1 (cryptic site)', 'Pks13 (covalent SuFEx target)']
                    }
            except Exception as e:
                logger.error(f"Error loading workflow: {e}")
        return {}

    def load_chembl_fasta(self) -> Optional[str]:
        """Load ChEMBL FASTA data"""
        fasta_file = self.base_path / "Downloads/chembl_35.fa.gz"
        if fasta_file.exists():
            try:
                with gzip.open(fasta_file, 'rt') as f:
                    return f.read()
            except Exception as e:
                logger.error(f"Error loading FASTA file: {e}")
        return None


class DataWarehouse:
    """Central data warehouse for integrated drug discovery data"""

    def __init__(self, db_path: str = "./crios_nova_warehouse.db"):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self._initialize_schema()

    def _initialize_schema(self):
        """Create database schema"""
        cursor = self.conn.cursor()

        # Compounds table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS compounds (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                smiles TEXT UNIQUE NOT NULL,
                inchi TEXT,
                inchi_key TEXT,
                chembl_id TEXT UNIQUE,
                pubchem_cid TEXT,
                drugbank_id TEXT,
                name TEXT,
                molecular_weight REAL,
                logp REAL,
                psa REAL,
                num_ro5_violations INTEGER,
                max_phase INTEGER,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)

        # Targets table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS targets (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                chembl_id TEXT UNIQUE NOT NULL,
                name TEXT,
                organism TEXT,
                target_type TEXT,
                tax_id INTEGER,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)

        # Activities table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS activities (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                compound_id INTEGER,
                target_id INTEGER,
                activity_type TEXT,
                activity_value REAL,
                activity_units TEXT,
                assay_chembl_id TEXT,
                document_id TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (compound_id) REFERENCES compounds(id),
                FOREIGN KEY (target_id) REFERENCES targets(id)
            )
        """)

        # Publications table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS publications (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pmid TEXT UNIQUE NOT NULL,
                title TEXT,
                abstract TEXT,
                authors TEXT,
                journal TEXT,
                pub_date TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)

        # Compound-Publication junction table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS compound_publications (
                compound_id INTEGER,
                publication_id INTEGER,
                relevance_score REAL,
                PRIMARY KEY (compound_id, publication_id),
                FOREIGN KEY (compound_id) REFERENCES compounds(id),
                FOREIGN KEY (publication_id) REFERENCES publications(id)
            )
        """)

        self.conn.commit()

    def insert_compound(self, compound_data: Dict) -> int:
        """Insert compound into database"""
        cursor = self.conn.cursor()

        try:
            cursor.execute("""
                INSERT OR REPLACE INTO compounds
                (smiles, inchi, inchi_key, chembl_id, pubchem_cid, drugbank_id,
                 name, molecular_weight, logp, psa, num_ro5_violations, max_phase)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                compound_data.get('smiles'),
                compound_data.get('inchi'),
                compound_data.get('inchi_key'),
                compound_data.get('chembl_id'),
                compound_data.get('pubchem_cid'),
                compound_data.get('drugbank_id'),
                compound_data.get('name'),
                compound_data.get('molecular_weight'),
                compound_data.get('logp'),
                compound_data.get('psa'),
                compound_data.get('num_ro5_violations'),
                compound_data.get('max_phase')
            ))
            self.conn.commit()
            return cursor.lastrowid
        except Exception as e:
            logger.error(f"Error inserting compound: {e}")
            return -1

    def insert_target(self, target_data: Dict) -> int:
        """Insert target into database"""
        cursor = self.conn.cursor()

        try:
            cursor.execute("""
                INSERT OR REPLACE INTO targets
                (chembl_id, name, organism, target_type, tax_id)
                VALUES (?, ?, ?, ?, ?)
            """, (
                target_data.get('chembl_id'),
                target_data.get('pref_name'),
                target_data.get('organism'),
                target_data.get('target_type'),
                target_data.get('tax_id')
            ))
            self.conn.commit()
            return cursor.lastrowid
        except Exception as e:
            logger.error(f"Error inserting target: {e}")
            return -1

    def query_compounds_by_target(self, target_chembl_id: str) -> pd.DataFrame:
        """Query compounds active against a specific target"""
        query = """
            SELECT c.*, a.activity_type, a.activity_value, a.activity_units
            FROM compounds c
            JOIN activities a ON c.id = a.compound_id
            JOIN targets t ON a.target_id = t.id
            WHERE t.chembl_id = ?
            ORDER BY a.activity_value ASC
        """
        return pd.read_sql_query(query, self.conn, params=(target_chembl_id,))

    def get_compound_profile(self, chembl_id: str) -> Dict:
        """Get complete compound profile"""
        cursor = self.conn.cursor()

        # Get compound data
        cursor.execute("SELECT * FROM compounds WHERE chembl_id = ?", (chembl_id,))
        compound = cursor.fetchone()

        if not compound:
            return {}

        # Get activities
        cursor.execute("""
            SELECT t.name, a.activity_type, a.activity_value, a.activity_units
            FROM activities a
            JOIN targets t ON a.target_id = t.id
            WHERE a.compound_id = ?
        """, (compound[0],))
        activities = cursor.fetchall()

        # Get publications
        cursor.execute("""
            SELECT p.pmid, p.title, p.journal, cp.relevance_score
            FROM publications p
            JOIN compound_publications cp ON p.id = cp.publication_id
            WHERE cp.compound_id = ?
            ORDER BY cp.relevance_score DESC
        """, (compound[0],))
        publications = cursor.fetchall()

        return {
            'compound': compound,
            'activities': activities,
            'publications': publications
        }

    def close(self):
        """Close database connection"""
        self.conn.close()


class IntegratedPipeline:
    """Main pipeline for integrating all data sources"""

    def __init__(self):
        self.chembl = ChEMBLIntegrator()
        self.pubmed = PubMedIntegrator()
        self.local = LocalDataIntegrator()
        self.warehouse = DataWarehouse()

    def integrate_all_sources(self):
        """Integrate data from all available sources"""
        logger.info("Starting data integration from all sources...")

        # Load local ChEMBL targets
        chembl_file = "C:/Users/micha/Downloads/chembl_targets_all.jsonl"
        if os.path.exists(chembl_file):
            targets_df = self.chembl.load_local_targets(chembl_file)
            logger.info(f"Loaded {len(targets_df)} targets from local ChEMBL file")

            # Insert targets into warehouse
            for _, target in targets_df.iterrows():
                self.warehouse.insert_target(target.to_dict())

        # Load local drug discovery examples
        synapse_examples = self.local.load_synapse_examples()
        logger.info(f"Loaded {len(synapse_examples)} Synapse language examples")

        # Load AI workflow
        ai_workflow = self.local.load_ai_workflow()
        if ai_workflow:
            logger.info("Loaded AI drug discovery workflow")

        logger.info("Data integration complete!")

    def search_compound_data(self, compound_name: str) -> IntegratedCompound:
        """Search and integrate data for a specific compound"""
        # Search ChEMBL
        chembl_data = {}
        try:
            molecules = self.chembl.molecule.search(compound_name)
            if molecules:
                chembl_data = self.chembl.fetch_molecule_properties(
                    molecules[0]['molecule_chembl_id']
                )
        except:
            pass

        # Search PubMed
        publications = self.pubmed.search_compound_literature(compound_name)

        # Create integrated compound object
        return IntegratedCompound(
            smiles=chembl_data.get('smiles', ''),
            inchi=chembl_data.get('inchi'),
            inchi_key=chembl_data.get('inchi_key'),
            chembl_id=chembl_data.get('chembl_id'),
            pubchem_cid=None,
            drugbank_id=None,
            name=compound_name,
            synonyms=[],
            targets=[],
            bioactivities=[],
            publications=publications,
            clinical_trials=[],
            properties=chembl_data,
            source_databases=['ChEMBL', 'PubMed'],
            last_updated=datetime.now()
        )


def main():
    """Main execution"""
    pipeline = IntegratedPipeline()
    pipeline.integrate_all_sources()

    # Example: Search for a specific compound
    compound_data = pipeline.search_compound_data("Imatinib")
    logger.info(f"Found data for compound: {compound_data.name}")
    logger.info(f"Publications found: {len(compound_data.publications)}")

    return pipeline


if __name__ == "__main__":
    pipeline = main()