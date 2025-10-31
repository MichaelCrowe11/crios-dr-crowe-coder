"""
CriOS Nova Cheminformatics - Machine Learning Pipeline for Compound Discovery
Integrates ChEMBL, PubMed, and drug discovery data for advanced molecular modeling
"""

import os
import json
import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass, field
from pathlib import Path
import joblib
from datetime import datetime

# Scientific computing
from scipy import stats
from scipy.spatial import distance

# RDKit for molecular operations
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, Crippen
from rdkit.Chem import rdMolDescriptors, DataStructs, rdFingerprintGenerator
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

# Machine Learning
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.ensemble import (
    RandomForestClassifier, RandomForestRegressor,
    GradientBoostingClassifier, GradientBoostingRegressor,
    ExtraTreesClassifier, ExtraTreesRegressor
)
from sklearn.svm import SVC, SVR
from sklearn.neural_network import MLPClassifier, MLPRegressor
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, mean_squared_error, r2_score, mean_absolute_error
)
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering

# Deep Learning
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from torch.optim import Adam, AdamW
from torch.optim.lr_scheduler import CosineAnnealingLR, ReduceLROnPlateau

# XGBoost and LightGBM for advanced gradient boosting
try:
    import xgboost as xgb
    XGBOOST_AVAILABLE = True
except ImportError:
    XGBOOST_AVAILABLE = False

try:
    import lightgbm as lgb
    LIGHTGBM_AVAILABLE = True
except ImportError:
    LIGHTGBM_AVAILABLE = False

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class CompoundData:
    """Container for compound information from various sources"""
    smiles: str
    chembl_id: Optional[str] = None
    pubmed_ids: List[str] = field(default_factory=list)
    target_info: Dict[str, Any] = field(default_factory=dict)
    activity_data: Dict[str, float] = field(default_factory=dict)
    descriptors: Optional[np.ndarray] = None
    fingerprint: Optional[np.ndarray] = None
    predicted_properties: Dict[str, float] = field(default_factory=dict)
    source: str = "unknown"


class MolecularFeaturizer:
    """Advanced molecular featurization for ML models"""

    def __init__(self, fingerprint_type: str = "morgan", fp_size: int = 2048):
        self.fingerprint_type = fingerprint_type
        self.fp_size = fp_size
        self.descriptor_calculator = self._init_descriptors()
        self.fingerprint_generator = self._init_fingerprint_generator()

    def _init_descriptors(self):
        """Initialize RDKit descriptor calculator"""
        descriptor_names = [desc[0] for desc in Descriptors.descList]
        return MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)

    def _init_fingerprint_generator(self):
        """Initialize fingerprint generator based on type"""
        generators = {
            'morgan': rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=self.fp_size),
            'topological': rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=self.fp_size),
            'avalon': rdFingerprintGenerator.GetAvalonGenerator(fpSize=self.fp_size),
            'rdkit': rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=self.fp_size),
            'atom_pair': rdFingerprintGenerator.GetAtomPairGenerator(fpSize=self.fp_size)
        }
        return generators.get(self.fingerprint_type, generators['morgan'])

    def featurize_molecule(self, smiles: str) -> Optional[Dict[str, np.ndarray]]:
        """Generate comprehensive molecular features"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            # Calculate descriptors
            descriptors = np.array(self.descriptor_calculator.CalcDescriptors(mol))

            # Generate fingerprint
            fp = self.fingerprint_generator.GetFingerprint(mol)
            fp_array = np.zeros(self.fp_size)
            DataStructs.ConvertToNumpyArray(fp, fp_array)

            # Calculate additional features
            features = {
                'descriptors': descriptors,
                'fingerprint': fp_array,
                'combined': np.concatenate([descriptors, fp_array]),
                'molecular_weight': Descriptors.MolWt(mol),
                'logp': Crippen.MolLogP(mol),
                'tpsa': rdMolDescriptors.CalcTPSA(mol),
                'num_rotatable_bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
                'num_h_donors': rdMolDescriptors.CalcNumHBD(mol),
                'num_h_acceptors': rdMolDescriptors.CalcNumHBA(mol),
                'num_aromatic_rings': rdMolDescriptors.CalcNumAromaticRings(mol),
                'qed': self._calculate_qed(mol),
                'sa_score': self._calculate_sa_score(mol)
            }

            return features

        except Exception as e:
            logger.error(f"Error featurizing molecule {smiles}: {e}")
            return None

    def _calculate_qed(self, mol):
        """Calculate Quantitative Estimate of Drug-likeness"""
        try:
            from rdkit.Chem import QED
            return QED.qed(mol)
        except:
            return 0.5

    def _calculate_sa_score(self, mol):
        """Calculate Synthetic Accessibility Score"""
        # Simplified SA score calculation
        complexity = (
            rdMolDescriptors.CalcNumBridgeheadAtoms(mol) * 2 +
            rdMolDescriptors.CalcNumSpiroAtoms(mol) * 2 +
            rdMolDescriptors.CalcNumHeteroatoms(mol) * 0.5 +
            rdMolDescriptors.CalcNumRings(mol) * 0.5
        )
        return min(10.0, max(1.0, complexity))


class ChEMBLDataLoader:
    """Load and process ChEMBL target and compound data"""

    def __init__(self, data_path: str):
        self.data_path = Path(data_path)
        self.targets_data = []
        self.compounds = []

    def load_targets(self, file_path: str) -> pd.DataFrame:
        """Load ChEMBL targets from JSONL file"""
        targets = []
        with open(file_path, 'r') as f:
            for line in f:
                try:
                    target = json.loads(line)
                    targets.append({
                        'chembl_id': target.get('target_chembl_id'),
                        'pref_name': target.get('pref_name'),
                        'organism': target.get('organism'),
                        'target_type': target.get('target_type'),
                        'tax_id': target.get('tax_id')
                    })
                except:
                    continue

        return pd.DataFrame(targets)

    def load_compound_activities(self, compounds_file: str) -> pd.DataFrame:
        """Load compound activity data"""
        # This would load actual compound-target activity data
        # For now, creating a placeholder structure
        return pd.DataFrame({
            'compound_chembl_id': [],
            'target_chembl_id': [],
            'standard_type': [],
            'standard_value': [],
            'standard_units': [],
            'activity_comment': []
        })


class GraphNeuralNetwork(nn.Module):
    """Graph Neural Network for molecular property prediction"""

    def __init__(self, input_dim: int, hidden_dim: int = 256, output_dim: int = 1,
                 num_layers: int = 3, dropout: float = 0.2):
        super(GraphNeuralNetwork, self).__init__()
        self.num_layers = num_layers
        self.dropout = dropout

        # Input layer
        self.input_layer = nn.Linear(input_dim, hidden_dim)

        # Hidden layers
        self.hidden_layers = nn.ModuleList([
            nn.Linear(hidden_dim, hidden_dim) for _ in range(num_layers - 1)
        ])

        # Batch normalization layers
        self.batch_norms = nn.ModuleList([
            nn.BatchNorm1d(hidden_dim) for _ in range(num_layers)
        ])

        # Output layer
        self.output_layer = nn.Linear(hidden_dim, output_dim)

        # Dropout
        self.dropout_layer = nn.Dropout(dropout)

    def forward(self, x):
        # Input transformation
        x = F.relu(self.batch_norms[0](self.input_layer(x)))
        x = self.dropout_layer(x)

        # Hidden layers with residual connections
        for i, (layer, bn) in enumerate(zip(self.hidden_layers, self.batch_norms[1:])):
            residual = x
            x = F.relu(bn(layer(x)))
            x = self.dropout_layer(x)
            if i % 2 == 0:  # Add residual connection every other layer
                x = x + residual

        # Output
        x = self.output_layer(x)
        return x


class CompoundDiscoveryPipeline:
    """Main pipeline for drug discovery ML workflows"""

    def __init__(self, config: Dict[str, Any] = None):
        self.config = config or self._default_config()
        self.featurizer = MolecularFeaturizer(
            fingerprint_type=self.config.get('fingerprint_type', 'morgan')
        )
        self.models = {}
        self.scalers = {}
        self.results = {}
        self.chembl_loader = None

    def _default_config(self) -> Dict[str, Any]:
        """Default configuration for the pipeline"""
        return {
            'fingerprint_type': 'morgan',
            'fp_size': 2048,
            'test_size': 0.2,
            'random_state': 42,
            'cv_folds': 5,
            'models': ['rf', 'xgb', 'nn'],
            'tasks': ['classification', 'regression'],
            'device': 'cuda' if torch.cuda.is_available() else 'cpu'
        }

    def load_chembl_data(self, chembl_file: str) -> pd.DataFrame:
        """Load ChEMBL data from file"""
        self.chembl_loader = ChEMBLDataLoader(os.path.dirname(chembl_file))
        targets_df = self.chembl_loader.load_targets(chembl_file)
        logger.info(f"Loaded {len(targets_df)} ChEMBL targets")
        return targets_df

    def prepare_training_data(self, compounds: List[CompoundData]) -> Tuple[np.ndarray, np.ndarray]:
        """Prepare training data from compound list"""
        X = []
        y = []

        for compound in compounds:
            features = self.featurizer.featurize_molecule(compound.smiles)
            if features is not None:
                X.append(features['combined'])
                # Extract target value (e.g., IC50, Ki, etc.)
                target_value = compound.activity_data.get('pIC50', 0)
                y.append(target_value)

        return np.array(X), np.array(y)

    def train_ensemble_models(self, X: np.ndarray, y: np.ndarray, task: str = 'regression'):
        """Train ensemble of ML models"""
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=self.config['test_size'],
            random_state=self.config['random_state']
        )

        # Scale features
        scaler = RobustScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        self.scalers[task] = scaler

        # Train models
        models_config = {
            'rf': RandomForestRegressor(n_estimators=200, max_depth=15, random_state=42)
                  if task == 'regression' else RandomForestClassifier(n_estimators=200, max_depth=15),
            'gb': GradientBoostingRegressor(n_estimators=150, learning_rate=0.1, max_depth=7)
                  if task == 'regression' else GradientBoostingClassifier(n_estimators=150, learning_rate=0.1),
            'et': ExtraTreesRegressor(n_estimators=200, max_depth=15, random_state=42)
                  if task == 'regression' else ExtraTreesClassifier(n_estimators=200, max_depth=15),
        }

        # Add XGBoost if available
        if XGBOOST_AVAILABLE:
            models_config['xgb'] = xgb.XGBRegressor(
                n_estimators=200, learning_rate=0.1, max_depth=7
            ) if task == 'regression' else xgb.XGBClassifier(
                n_estimators=200, learning_rate=0.1, max_depth=7
            )

        # Add LightGBM if available
        if LIGHTGBM_AVAILABLE:
            models_config['lgb'] = lgb.LGBMRegressor(
                n_estimators=200, learning_rate=0.1, max_depth=7
            ) if task == 'regression' else lgb.LGBMClassifier(
                n_estimators=200, learning_rate=0.1, max_depth=7
            )

        results = {}
        for name, model in models_config.items():
            logger.info(f"Training {name} model for {task}...")

            # Train model
            model.fit(X_train_scaled, y_train)

            # Evaluate
            y_pred = model.predict(X_test_scaled)

            if task == 'regression':
                metrics = {
                    'mse': mean_squared_error(y_test, y_pred),
                    'rmse': np.sqrt(mean_squared_error(y_test, y_pred)),
                    'mae': mean_absolute_error(y_test, y_pred),
                    'r2': r2_score(y_test, y_pred)
                }
            else:
                metrics = {
                    'accuracy': accuracy_score(y_test, y_pred),
                    'precision': precision_score(y_test, y_pred, average='weighted'),
                    'recall': recall_score(y_test, y_pred, average='weighted'),
                    'f1': f1_score(y_test, y_pred, average='weighted')
                }
                if len(np.unique(y)) == 2:  # Binary classification
                    y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
                    metrics['auc'] = roc_auc_score(y_test, y_pred_proba)

            results[name] = {
                'model': model,
                'metrics': metrics,
                'feature_importance': self._get_feature_importance(model, name)
            }

            logger.info(f"{name} metrics: {metrics}")

        self.models[task] = results
        return results

    def _get_feature_importance(self, model, model_name: str) -> Optional[np.ndarray]:
        """Extract feature importance from trained model"""
        try:
            if hasattr(model, 'feature_importances_'):
                return model.feature_importances_
            elif model_name == 'xgb' and XGBOOST_AVAILABLE:
                return model.feature_importances_
            elif model_name == 'lgb' and LIGHTGBM_AVAILABLE:
                return model.feature_importances_
        except:
            pass
        return None

    def train_deep_learning_model(self, X: np.ndarray, y: np.ndarray,
                                  epochs: int = 100, batch_size: int = 32):
        """Train deep learning model for molecular property prediction"""
        device = torch.device(self.config['device'])

        # Prepare data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=self.config['test_size'],
            random_state=self.config['random_state']
        )

        # Scale features
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)

        # Convert to tensors
        X_train_tensor = torch.FloatTensor(X_train_scaled).to(device)
        y_train_tensor = torch.FloatTensor(y_train.reshape(-1, 1)).to(device)
        X_test_tensor = torch.FloatTensor(X_test_scaled).to(device)
        y_test_tensor = torch.FloatTensor(y_test.reshape(-1, 1)).to(device)

        # Create model
        input_dim = X_train.shape[1]
        model = GraphNeuralNetwork(
            input_dim=input_dim,
            hidden_dim=256,
            output_dim=1,
            num_layers=4,
            dropout=0.2
        ).to(device)

        # Training setup
        criterion = nn.MSELoss()
        optimizer = AdamW(model.parameters(), lr=0.001, weight_decay=1e-5)
        scheduler = CosineAnnealingLR(optimizer, T_max=epochs, eta_min=1e-6)

        # Training loop
        model.train()
        train_losses = []

        for epoch in range(epochs):
            # Forward pass
            outputs = model(X_train_tensor)
            loss = criterion(outputs, y_train_tensor)

            # Backward pass
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            scheduler.step()

            train_losses.append(loss.item())

            if (epoch + 1) % 10 == 0:
                logger.info(f"Epoch [{epoch+1}/{epochs}], Loss: {loss.item():.4f}")

        # Evaluation
        model.eval()
        with torch.no_grad():
            test_outputs = model(X_test_tensor)
            test_loss = criterion(test_outputs, y_test_tensor)

            # Convert back to numpy for metrics
            y_pred = test_outputs.cpu().numpy()
            y_true = y_test_tensor.cpu().numpy()

            metrics = {
                'mse': mean_squared_error(y_true, y_pred),
                'rmse': np.sqrt(mean_squared_error(y_true, y_pred)),
                'mae': mean_absolute_error(y_true, y_pred),
                'r2': r2_score(y_true, y_pred)
            }

        logger.info(f"Deep Learning Model - Test Metrics: {metrics}")

        self.models['deep_learning'] = {
            'model': model,
            'scaler': scaler,
            'metrics': metrics,
            'training_history': train_losses
        }

        return model, metrics

    def virtual_screening(self, compound_library: List[str],
                         target_model: str = 'rf',
                         top_k: int = 100) -> List[Tuple[str, float]]:
        """Screen compound library using trained models"""
        if target_model not in self.models.get('regression', {}):
            raise ValueError(f"Model {target_model} not found. Train models first.")

        model = self.models['regression'][target_model]['model']
        scaler = self.scalers.get('regression')

        predictions = []
        for smiles in compound_library:
            features = self.featurizer.featurize_molecule(smiles)
            if features is not None:
                X = features['combined'].reshape(1, -1)
                X_scaled = scaler.transform(X)
                pred = model.predict(X_scaled)[0]
                predictions.append((smiles, pred))

        # Sort by predicted activity (higher is better for pIC50)
        predictions.sort(key=lambda x: x[1], reverse=True)

        return predictions[:top_k]

    def scaffold_analysis(self, compounds: List[str]) -> Dict[str, List[str]]:
        """Analyze molecular scaffolds in compound set"""
        scaffolds = {}

        for smiles in compounds:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_smiles = Chem.MolToSmiles(scaffold)

                if scaffold_smiles not in scaffolds:
                    scaffolds[scaffold_smiles] = []
                scaffolds[scaffold_smiles].append(smiles)

        # Sort by frequency
        sorted_scaffolds = dict(sorted(scaffolds.items(),
                                     key=lambda x: len(x[1]),
                                     reverse=True))

        return sorted_scaffolds

    def cluster_compounds(self, compounds: List[str], n_clusters: int = 10) -> Dict[int, List[str]]:
        """Cluster compounds based on molecular similarity"""
        # Generate features for all compounds
        features = []
        valid_compounds = []

        for smiles in compounds:
            feat = self.featurizer.featurize_molecule(smiles)
            if feat is not None:
                features.append(feat['fingerprint'])
                valid_compounds.append(smiles)

        if not features:
            return {}

        X = np.array(features)

        # Perform clustering
        kmeans = KMeans(n_clusters=min(n_clusters, len(valid_compounds)),
                        random_state=42)
        cluster_labels = kmeans.fit_predict(X)

        # Organize results
        clusters = {}
        for compound, label in zip(valid_compounds, cluster_labels):
            if label not in clusters:
                clusters[label] = []
            clusters[label].append(compound)

        return clusters

    def save_pipeline(self, output_dir: str):
        """Save trained pipeline to disk"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Save models
        for task, models in self.models.items():
            for model_name, model_data in models.items():
                if model_name != 'deep_learning':
                    model_file = output_path / f"{task}_{model_name}_model.pkl"
                    joblib.dump(model_data['model'], model_file)
                else:
                    # Save PyTorch model
                    model_file = output_path / f"{task}_deep_learning_model.pt"
                    torch.save(model_data['model'].state_dict(), model_file)

        # Save scalers
        for task, scaler in self.scalers.items():
            scaler_file = output_path / f"{task}_scaler.pkl"
            joblib.dump(scaler, scaler_file)

        # Save configuration
        config_file = output_path / "pipeline_config.json"
        with open(config_file, 'w') as f:
            json.dump(self.config, f, indent=2)

        logger.info(f"Pipeline saved to {output_path}")

    def load_pipeline(self, input_dir: str):
        """Load saved pipeline from disk"""
        input_path = Path(input_dir)

        # Load configuration
        config_file = input_path / "pipeline_config.json"
        if config_file.exists():
            with open(config_file, 'r') as f:
                self.config = json.load(f)

        # Load models
        for model_file in input_path.glob("*_model.pkl"):
            parts = model_file.stem.split('_')
            task = parts[0]
            model_name = '_'.join(parts[1:-1])

            if task not in self.models:
                self.models[task] = {}

            self.models[task][model_name] = {
                'model': joblib.load(model_file),
                'metrics': {},
                'feature_importance': None
            }

        # Load scalers
        for scaler_file in input_path.glob("*_scaler.pkl"):
            task = scaler_file.stem.replace('_scaler', '')
            self.scalers[task] = joblib.load(scaler_file)

        logger.info(f"Pipeline loaded from {input_path}")


def main():
    """Main execution function"""
    # Initialize pipeline
    pipeline = CompoundDiscoveryPipeline({
        'fingerprint_type': 'morgan',
        'fp_size': 2048,
        'models': ['rf', 'gb', 'xgb', 'lgb'],
        'device': 'cuda' if torch.cuda.is_available() else 'cpu'
    })

    # Load ChEMBL data
    chembl_file = "C:/Users/micha/Downloads/chembl_targets_all.jsonl"
    if os.path.exists(chembl_file):
        targets_df = pipeline.load_chembl_data(chembl_file)
        logger.info(f"Loaded {len(targets_df)} targets from ChEMBL")

    # Example: Create synthetic training data for demonstration
    # In production, this would come from actual ChEMBL activity data
    example_compounds = [
        CompoundData(
            smiles="CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
            chembl_id="CHEMBL25",
            activity_data={'pIC50': 5.2}
        ),
        CompoundData(
            smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
            chembl_id="CHEMBL113",
            activity_data={'pIC50': 4.8}
        ),
        CompoundData(
            smiles="CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
            chembl_id="CHEMBL521",
            activity_data={'pIC50': 5.5}
        )
    ]

    # Prepare training data
    if example_compounds:
        X, y = pipeline.prepare_training_data(example_compounds)

        if len(X) > 0:
            # Train ensemble models
            results = pipeline.train_ensemble_models(X, y, task='regression')

            # Train deep learning model
            if len(X) >= 10:  # Need enough samples for DL
                model, metrics = pipeline.train_deep_learning_model(X, y, epochs=50)

            # Save pipeline
            pipeline.save_pipeline("./crios_nova_ml_models")

    logger.info("CriOS Nova ML Pipeline initialized successfully!")

    return pipeline


if __name__ == "__main__":
    pipeline = main()