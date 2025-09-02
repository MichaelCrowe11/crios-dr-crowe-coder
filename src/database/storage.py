import sqlite3
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import json
import pickle
import logging
from datetime import datetime
from contextlib import contextmanager

logger = logging.getLogger(__name__)


class CompoundStorage:
    """SQLite storage for compound library"""
    
    def __init__(self, db_path: str = "compound_library.db"):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_database()
    
    @contextmanager
    def get_connection(self):
        """Context manager for database connections"""
        conn = sqlite3.connect(str(self.db_path))
        conn.row_factory = sqlite3.Row
        try:
            yield conn
        finally:
            conn.close()
    
    def _init_database(self):
        """Initialize database schema"""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            # Main compounds table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS compounds (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    smiles TEXT UNIQUE NOT NULL,
                    canonical_smiles TEXT,
                    mol_id TEXT,
                    mw REAL,
                    logp REAL,
                    tpsa REAL,
                    hba INTEGER,
                    hbd INTEGER,
                    rotatable_bonds INTEGER,
                    fingerprint BLOB,
                    source TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
            
            # Descriptors table for flexible storage
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS descriptors (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    compound_id INTEGER NOT NULL,
                    descriptor_name TEXT NOT NULL,
                    descriptor_value REAL,
                    FOREIGN KEY (compound_id) REFERENCES compounds(id),
                    UNIQUE(compound_id, descriptor_name)
                )
            """)
            
            # Properties table for additional metadata
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS properties (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    compound_id INTEGER NOT NULL,
                    property_name TEXT NOT NULL,
                    property_value TEXT,
                    FOREIGN KEY (compound_id) REFERENCES compounds(id),
                    UNIQUE(compound_id, property_name)
                )
            """)
            
            # Activities table for bioassay data
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS activities (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    compound_id INTEGER NOT NULL,
                    target TEXT,
                    activity_type TEXT,
                    activity_value REAL,
                    activity_units TEXT,
                    assay_id TEXT,
                    source TEXT,
                    FOREIGN KEY (compound_id) REFERENCES compounds(id)
                )
            """)
            
            # Clusters table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS clusters (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    cluster_name TEXT UNIQUE,
                    cluster_method TEXT,
                    parameters TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
            
            # Cluster members table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS cluster_members (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    cluster_id INTEGER NOT NULL,
                    compound_id INTEGER NOT NULL,
                    is_representative BOOLEAN DEFAULT 0,
                    FOREIGN KEY (cluster_id) REFERENCES clusters(id),
                    FOREIGN KEY (compound_id) REFERENCES compounds(id),
                    UNIQUE(cluster_id, compound_id)
                )
            """)
            
            # Create indices for better performance
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_smiles ON compounds(smiles)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_mw ON compounds(mw)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_logp ON compounds(logp)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_source ON compounds(source)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_descriptors ON descriptors(compound_id, descriptor_name)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_activities ON activities(compound_id, target)")
            
            conn.commit()
    
    def add_compound(self, smiles: str, mol_id: Optional[str] = None,
                    descriptors: Optional[Dict[str, float]] = None,
                    properties: Optional[Dict[str, Any]] = None,
                    fingerprint: Optional[Any] = None,
                    source: Optional[str] = None) -> int:
        """Add a compound to the database"""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            # Insert main compound data
            try:
                cursor.execute("""
                    INSERT INTO compounds (smiles, mol_id, mw, logp, tpsa, hba, hbd, 
                                         rotatable_bonds, fingerprint, source)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    smiles,
                    mol_id,
                    descriptors.get('MW') if descriptors else None,
                    descriptors.get('LogP') if descriptors else None,
                    descriptors.get('TPSA') if descriptors else None,
                    descriptors.get('HBA') if descriptors else None,
                    descriptors.get('HBD') if descriptors else None,
                    descriptors.get('RotatableBonds') if descriptors else None,
                    pickle.dumps(fingerprint) if fingerprint else None,
                    source
                ))
                
                compound_id = cursor.lastrowid
                
                # Add additional descriptors
                if descriptors:
                    for name, value in descriptors.items():
                        if name not in ['MW', 'LogP', 'TPSA', 'HBA', 'HBD', 'RotatableBonds']:
                            cursor.execute("""
                                INSERT OR REPLACE INTO descriptors (compound_id, descriptor_name, descriptor_value)
                                VALUES (?, ?, ?)
                            """, (compound_id, name, value))
                
                # Add properties
                if properties:
                    for name, value in properties.items():
                        cursor.execute("""
                            INSERT OR REPLACE INTO properties (compound_id, property_name, property_value)
                            VALUES (?, ?, ?)
                        """, (compound_id, name, json.dumps(value) if not isinstance(value, str) else value))
                
                conn.commit()
                return compound_id
                
            except sqlite3.IntegrityError:
                # Compound already exists, get its ID
                cursor.execute("SELECT id FROM compounds WHERE smiles = ?", (smiles,))
                row = cursor.fetchone()
                return row['id'] if row else -1
    
    def get_compound(self, compound_id: Optional[int] = None, 
                    smiles: Optional[str] = None) -> Optional[Dict[str, Any]]:
        """Get compound by ID or SMILES"""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            if compound_id:
                cursor.execute("SELECT * FROM compounds WHERE id = ?", (compound_id,))
            elif smiles:
                cursor.execute("SELECT * FROM compounds WHERE smiles = ?", (smiles,))
            else:
                return None
            
            row = cursor.fetchone()
            if not row:
                return None
            
            compound = dict(row)
            
            # Get additional descriptors
            cursor.execute("""
                SELECT descriptor_name, descriptor_value 
                FROM descriptors WHERE compound_id = ?
            """, (compound['id'],))
            
            descriptors = {row['descriptor_name']: row['descriptor_value'] 
                          for row in cursor.fetchall()}
            compound['descriptors'] = descriptors
            
            # Get properties
            cursor.execute("""
                SELECT property_name, property_value 
                FROM properties WHERE compound_id = ?
            """, (compound['id'],))
            
            properties = {}
            for row in cursor.fetchall():
                try:
                    properties[row['property_name']] = json.loads(row['property_value'])
                except:
                    properties[row['property_name']] = row['property_value']
            compound['properties'] = properties
            
            # Deserialize fingerprint if present
            if compound.get('fingerprint'):
                compound['fingerprint'] = pickle.loads(compound['fingerprint'])
            
            return compound
    
    def search_compounds(self, filters: Dict[str, Any], limit: int = 100) -> List[Dict[str, Any]]:
        """Search compounds with filters"""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            where_clauses = []
            params = []
            
            # Build WHERE clause
            for field, value in filters.items():
                if field in ['mw', 'logp', 'tpsa', 'hba', 'hbd', 'rotatable_bonds']:
                    if isinstance(value, tuple) and len(value) == 2:
                        # Range query
                        if value[0] is not None:
                            where_clauses.append(f"{field} >= ?")
                            params.append(value[0])
                        if value[1] is not None:
                            where_clauses.append(f"{field} <= ?")
                            params.append(value[1])
                    else:
                        where_clauses.append(f"{field} = ?")
                        params.append(value)
                elif field == 'source':
                    where_clauses.append("source = ?")
                    params.append(value)
            
            query = "SELECT * FROM compounds"
            if where_clauses:
                query += " WHERE " + " AND ".join(where_clauses)
            query += f" LIMIT {limit}"
            
            cursor.execute(query, params)
            
            compounds = []
            for row in cursor.fetchall():
                compound = dict(row)
                compounds.append(compound)
            
            return compounds
    
    def add_activity(self, compound_id: int, target: str, activity_type: str,
                    activity_value: float, activity_units: str = None,
                    assay_id: str = None, source: str = None):
        """Add bioactivity data for a compound"""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO activities (compound_id, target, activity_type, 
                                      activity_value, activity_units, assay_id, source)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (compound_id, target, activity_type, activity_value, 
                 activity_units, assay_id, source))
            conn.commit()
    
    def save_clustering(self, cluster_name: str, method: str, clusters: List[List[int]],
                       parameters: Dict[str, Any] = None) -> int:
        """Save clustering results"""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            # Save cluster metadata
            cursor.execute("""
                INSERT INTO clusters (cluster_name, cluster_method, parameters)
                VALUES (?, ?, ?)
            """, (cluster_name, method, json.dumps(parameters) if parameters else None))
            
            cluster_id = cursor.lastrowid
            
            # Save cluster members
            for cluster_idx, member_ids in enumerate(clusters):
                for idx, compound_id in enumerate(member_ids):
                    is_representative = (idx == 0)  # First member as representative
                    cursor.execute("""
                        INSERT INTO cluster_members (cluster_id, compound_id, is_representative)
                        VALUES (?, ?, ?)
                    """, (cluster_id, compound_id, is_representative))
            
            conn.commit()
            return cluster_id
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get database statistics"""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            stats = {}
            
            # Total compounds
            cursor.execute("SELECT COUNT(*) as count FROM compounds")
            stats['total_compounds'] = cursor.fetchone()['count']
            
            # Compounds by source
            cursor.execute("""
                SELECT source, COUNT(*) as count 
                FROM compounds 
                GROUP BY source
            """)
            stats['by_source'] = {row['source']: row['count'] for row in cursor.fetchall()}
            
            # MW distribution
            cursor.execute("""
                SELECT MIN(mw) as min_mw, MAX(mw) as max_mw, AVG(mw) as avg_mw
                FROM compounds WHERE mw IS NOT NULL
            """)
            row = cursor.fetchone()
            if row:
                stats['mw_stats'] = dict(row)
            
            # LogP distribution
            cursor.execute("""
                SELECT MIN(logp) as min_logp, MAX(logp) as max_logp, AVG(logp) as avg_logp
                FROM compounds WHERE logp IS NOT NULL
            """)
            row = cursor.fetchone()
            if row:
                stats['logp_stats'] = dict(row)
            
            # Activity data
            cursor.execute("SELECT COUNT(DISTINCT compound_id) as count FROM activities")
            stats['compounds_with_activities'] = cursor.fetchone()['count']
            
            cursor.execute("SELECT COUNT(DISTINCT target) as count FROM activities")
            stats['unique_targets'] = cursor.fetchone()['count']
            
            return stats