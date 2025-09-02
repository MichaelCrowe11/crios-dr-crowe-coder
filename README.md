# CriOS - Crowe Discovery Framework

## 🌿 Regenerative Discovery System for Molecular Intelligence

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![CriOS Nova](https://img.shields.io/badge/CriOS%20Nova-AI%20Agent-purple.svg)](docs/crowe-discovery-framework.md)

**CriOS (Crowe Research Intelligence OS)** is the computational foundation of the **Crowe Discovery Framework** - a regenerative discovery system that combines AI-guided hypothesis generation with live experimental validation to accelerate breakthrough compounds across neurodegeneration, environmental health, and sustainable materials.

> *Science before status. Discovery before profit. Earth and people above extraction.*

## 🚀 Quick Start

```bash
# Clone repository
git clone https://github.com/michaelcrowe11/crios.git
cd crios

# Install dependencies
pip install -e .

# Run demo
python scripts/demo.py

# Start CLI
crios --help

# Start web interface
crios web serve
```

## 🧠 CriOS Discovery Engine Architecture

### Universal Compound Discovery Platform

**CriOS** is a production-grade compound discovery platform covering:

- **🌿 Natural Products**: Mining fungal metabolites, plant compounds, and microbial libraries for bioactive scaffolds
- **⚗️ Synthetic Chemistry**: AI-guided design of novel small molecules with optimized drug-like properties  
- **🔬 Hybrid Approaches**: Nature-inspired synthetic analogs combining natural scaffolds with synthetic optimization
- **🎯 Multi-Target Discovery**: Compounds addressing neurodegeneration, cancer, aging, and metabolic disorders

### Core Capabilities

#### Molecular Processing Engine
- **Structure Validation**: RDKit-based SMILES parsing and molecular validation
- **Descriptor Calculation**: 50+ molecular descriptors including custom Crowe metrics
- **Fingerprint Generation**: Morgan, RDKit, and MACCS molecular fingerprints
- **Similarity Search**: High-performance Tanimoto similarity with biological scaling

#### Discovery Algorithms  
- **Crowe Scoring**: Multi-objective compound prioritization (novelty, drug-likeness, safety, ethics)
- **Clustering Analysis**: Butina and hierarchical clustering for chemical space exploration
- **Diversity Selection**: MaxMin algorithm for diverse compound subset selection
- **Multi-Objective Optimization**: Pareto ranking for balanced compound optimization

#### Ethical AI Framework
- **Mission-Locked IP**: Ethical licensing ensuring science serves healing
- **Safety Screening**: Automated toxicity and dual-use assessment
- **Bias Mitigation**: Fairness metrics and diverse training data requirements
- **Transparency**: Explainable AI with confidence quantification

## 📖 System Architecture

```
CriOS/
├── crios/                   # Core Python package
│   ├── core/               # Molecular processing engine
│   │   ├── molecule.py     # RDKit molecular representation
│   │   ├── compound.py     # Extended compound with bio data
│   │   ├── similarity.py   # Similarity search engine
│   │   ├── descriptors.py  # Descriptor calculation
│   │   └── clustering.py   # Clustering algorithms
│   ├── scoring/            # Multi-objective scoring
│   │   ├── crowe_score.py  # Crowe Discovery methodology
│   │   ├── drug_likeness.py# ADMET and drug-like filters
│   │   ├── novelty.py      # Structural novelty assessment
│   │   └── ethics.py       # Ethical compliance scoring
│   ├── cli/               # Command line interface
│   │   ├── main.py        # Typer-based CLI
│   │   └── commands/      # Individual command modules
│   ├── web/               # FastAPI web interface
│   │   ├── main.py        # REST API endpoints
│   │   └── models.py      # Pydantic data models
│   └── config/            # Configuration management
│       ├── config.yaml    # System configuration
│       ├── ethics.yaml    # Ethical policies
│       └── schemas.py     # Pydantic config models
├── platform/              # Web platform (Next.js)
│   ├── frontend/          # React components
│   │   └── components/    # Immersive IDE components
│   └── backend/           # Additional FastAPI services
├── scripts/               # Utility scripts
│   └── demo.py           # Complete system demonstration
├── docs/                  # Documentation
│   ├── crowe-discovery-framework.md
│   └── implementation-roadmap.md
└── tests/                # Comprehensive test suite
```

## 🔧 Command Line Interface

### Core Discovery Commands

```bash
# Validate molecular structures
crios validate -i compounds.csv --standardize

# Calculate molecular descriptors
crios process -i molecules.sdf --descriptors --crowe-score

# Perform similarity search
crios search -q "CCO" -d database.sdf -t 0.7

# Cluster compounds by structure
crios cluster -i compounds.sdf -a butina -t 0.6

# Generate comprehensive analysis
crios analyze -i library.sdf --type detailed

# Database operations
crios db import -i compounds.sdf
crios db export -o results.sdf

# Web interface
crios web serve --port 8000
```

### Natural Product Discovery Workflow
```bash
# Extract bioactive scaffolds from natural sources
crios process --source fungi --bioactivity neuroprotective
crios search --scaffold-similarity --extract-pharmacophores
crios analyze --natural-product-score --optimize-synthesis
```

### Synthetic Chemistry Design Workflow  
```bash
# AI-guided small molecule design with drug-like optimization
crios generate --target GPCR --properties "BBB,oral-bioavailability"
crios process --lipinski-compliant --synthetic-accessibility
crios score --crowe-methodology --multi-objective
```

## 🌐 Web Interface & API

### REST API Endpoints

```python
# Molecule validation
POST /validate
{
  "smiles": "CCO",
  "standardize": true
}

# Similarity search
POST /similarity
{
  "query_smiles": "CCO",
  "database_smiles": ["CC(C)O", "CCC"],
  "threshold": 0.7
}

# Crowe scoring
POST /crowe-score
{
  "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
  "compound_origin": "synthetic",
  "therapeutic_areas": ["neurodegeneration"]
}

# Ethics validation
POST /ethics-check
{
  "smiles": "CCO"
}
```

### Web Interface
```bash
# Launch full web platform
crios web serve

# Access at http://localhost:8000
# - Interactive compound analysis
# - Real-time similarity search  
# - Crowe scoring dashboard
# - Ethical compliance monitoring
```

## 🧪 Python SDK Usage

```python
from crios import Molecule, CompoundLibrary, CroweScorer, SimilaritySearch

# Create and validate molecules
molecule = Molecule("CCO")  # Ethanol
print(f"Valid: {molecule.is_valid()}")
print(f"Drug-like: {molecule.is_drug_like()}")
print(f"MW: {molecule.properties.molecular_weight}")

# Build compound library
from crios.core.compound import Compound, CompoundOrigin
compound = Compound("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", origin=CompoundOrigin.SYNTHETIC)
compound.metadata.therapeutic_areas = ["neurodegeneration"]

library = CompoundLibrary([compound], name="Demo Library")

# Similarity search
search = SimilaritySearch(fingerprint_type="morgan")
results = search.search_single(molecule, library, threshold=0.7)

# Crowe scoring
scorer = CroweScorer()
score, components = scorer.score_single(compound, return_components=True)
print(f"Crowe Score: {score:.3f}")
print(f"Novelty: {components.novelty:.3f}")
print(f"Drug-likeness: {components.drug_likeness:.3f}")
```

## 📊 Performance Benchmarks

| Operation | Dataset Size | Processing Time | Memory Usage |
|-----------|--------------|-----------------|---------------|
| SMILES Validation | 10,000 compounds | 8 seconds | 120 MB |
| Similarity Search | 50,000 compounds | 12 seconds | 240 MB |
| Crowe Scoring | 5,000 compounds | 25 seconds | 180 MB |
| Clustering (Butina) | 5,000 compounds | 15 seconds | 200 MB |
| Descriptor Calculation | 10,000 compounds | 35 seconds | 280 MB |

*Benchmarked on Intel i7-8700K, 16GB RAM*

## ⚖️ Ethical AI Framework

CriOS implements comprehensive ethical constraints:

### Core Principles
- **Science before status**
- **Discovery before profit** 
- **Earth and people above extraction**
- **Healing over harm**

### Safety & Compliance
- Automated dual-use research screening
- Toxicity prediction and safety categorization
- Patent landscape and freedom-to-operate analysis
- Environmental impact assessment

### Mission-Locked IP
All discoveries flow into an ethical IP structure ensuring:
- Fair licensing terms with profit caps
- Public benefit clauses
- Research exemptions
- Open access preferences

## 🤝 Contributing

We welcome contributions aligned with our ethical framework:

1. **Fork the repository**
2. **Create feature branch** (`git checkout -b feature/amazing-discovery`)
3. **Commit changes** (`git commit -m 'Add amazing discovery'`)
4. **Push to branch** (`git push origin feature/amazing-discovery`)
5. **Open Pull Request**

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

## 📄 License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## 🌟 Citation

If you use CriOS in your research, please cite:

```bibtex
@software{crios2025,
  title={CriOS: Crowe Discovery Framework for Universal Compound Discovery},
  author={Crowe, Michael B.},
  year={2025},
  url={https://github.com/michaelcrowe11/crios},
  note={Universal compound discovery platform with ethical AI constraints}
}
```

## 📬 Contact

- **Website**: [www.crios.ai](https://www.crios.ai)
- **Email**: michael@crowelogic.com
- **Issues**: [GitHub Issues](https://github.com/michaelcrowe11/crios/issues)

## 🙏 Acknowledgments

Built with:
- [RDKit](https://www.rdkit.org/) - Cheminformatics toolkit
- [FastAPI](https://fastapi.tiangolo.com/) - Modern web framework
- [Typer](https://typer.tiangolo.com/) - CLI framework
- [Rich](https://github.com/Textualize/rich) - Beautiful terminal output
- [Pydantic](https://pydantic-docs.helpmanual.io/) - Data validation

---

*"We are not building a lab. We are seeding a forest of discovery."*

**Crowe Discovery Framework** - Where pain becomes pattern, and the wisdom of collapse births the architectures of renewal.

© 2025 Crowe BioSystems. Mission-locked IP structure ensures science is not siloed, stolen, or scaled without soul.