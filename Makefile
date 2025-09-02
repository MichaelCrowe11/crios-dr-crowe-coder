# CriOS Discovery Engine - Development Makefile
# Automation for development, testing, and deployment

.PHONY: help install dev test demo docs clean lint format type-check security
.DEFAULT_GOAL := help

# Variables
PYTHON = python
PIP = pip
PYTEST = pytest
BLACK = black
ISORT = isort
FLAKE8 = flake8
MYPY = mypy
BANDIT = bandit

# Source directories
SRC_DIRS = crios tests scripts

help: ## Show this help message
	@echo "🧬 CriOS Discovery Engine - Development Commands"
	@echo ""
	@echo "Usage: make <target>"
	@echo ""
	@echo "Targets:"
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

# Installation targets
install: ## Install package in development mode
	$(PIP) install -e .

install-dev: ## Install development dependencies
	$(PIP) install -e ".[dev]"

install-all: ## Install all optional dependencies
	$(PIP) install -e ".[all]"

# Development targets
dev: install-dev ## Set up development environment
	pre-commit install
	@echo "✅ Development environment ready!"
	@echo "💡 Try: make demo"

# Testing targets
test: ## Run all tests
	$(PYTEST) tests/

test-unit: ## Run unit tests only
	$(PYTEST) tests/ -m "unit"

test-integration: ## Run integration tests only
	$(PYTEST) tests/ -m "integration"

test-fast: ## Run fast tests (exclude slow tests)
	$(PYTEST) tests/ -m "not slow"

test-cov: ## Run tests with coverage report
	$(PYTEST) tests/ --cov=crios --cov-report=html --cov-report=term

# Code quality targets
lint: ## Run all linting checks
	$(FLAKE8) $(SRC_DIRS)
	$(BLACK) --check $(SRC_DIRS)
	$(ISORT) --check-only $(SRC_DIRS)

format: ## Format code with black and isort
	$(BLACK) $(SRC_DIRS)
	$(ISORT) $(SRC_DIRS)
	@echo "✅ Code formatted!"

type-check: ## Run type checking with mypy
	$(MYPY) crios/

security: ## Run security checks with bandit
	$(BANDIT) -r crios/ -f json -o security-report.json || true
	$(BANDIT) -r crios/ || true

# Quality assurance (run all checks)
qa: lint type-check security test-fast ## Run all quality assurance checks
	@echo "✅ Quality assurance complete!"

# Demo and examples
demo: ## Run end-to-end demonstration
	@echo "🧬 Running CriOS Discovery Engine Demo..."
	@echo ""
	
	# Initialize project
	@echo "📁 Initializing demo project..."
	@mkdir -p demo_run
	@cd demo_run && $(PYTHON) -m crios.cli.main init --overwrite
	
	# Run featurization
	@echo ""
	@echo "🔬 Running featurization demo..."
	@cd demo_run && $(PYTHON) -m crios.cli.main featurize --in data/examples/molecules.smi --out data/outputs/features.csv
	
	# Run scoring
	@echo ""
	@echo "📊 Running scoring demo..."
	@cd demo_run && $(PYTHON) -m crios.cli.main score --in data/examples/molecules.smi --out data/outputs/scored.csv --engine synthetic --target-class GPCR
	
	# Run ethics check
	@echo ""
	@echo "⚖️ Running ethics compliance demo..."
	@cd demo_run && $(PYTHON) -m crios.cli.main explain --in data/examples/molecules.smi --policy configs/ethics.yaml --out data/outputs/ethics_report.json
	
	@echo ""
	@echo "🎉 Demo complete! Check demo_run/data/outputs/ for results."

demo-clean: ## Clean up demo artifacts
	rm -rf demo_run/
	@echo "🧹 Demo artifacts cleaned!"

# Python demo script
demo-python: ## Run Python SDK demonstration
	$(PYTHON) scripts/demo.py

# Documentation targets
docs: ## Build documentation
	@echo "📚 Building documentation..."
	@echo "💡 Documentation build not yet implemented"
	@echo "📋 TODO: Set up mkdocs or sphinx"

docs-serve: ## Serve documentation locally
	@echo "🌐 Serving documentation..."
	@echo "💡 Documentation serving not yet implemented"

# Web interface targets
web: ## Start web interface
	@echo "🌐 Starting CriOS web interface..."
	$(PYTHON) -m uvicorn crios.api.main:app --reload --port 8000

web-prod: ## Start production web server
	@echo "🚀 Starting production web server..."
	$(PYTHON) -m uvicorn crios.api.main:app --host 0.0.0.0 --port 8000

# Build and distribution targets
build: clean ## Build distribution packages
	$(PYTHON) -m build

dist: build ## Create distribution (alias for build)

upload: build ## Upload to PyPI (requires twine)
	twine upload dist/*

upload-test: build ## Upload to Test PyPI
	twine upload --repository testpypi dist/*

# Maintenance targets
clean: ## Clean build artifacts and cache files
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .coverage
	rm -rf htmlcov/
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf security-report.json
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	@echo "🧹 Cleaned build artifacts!"

clean-all: clean demo-clean ## Clean everything including demo artifacts

# Development workflow targets
check: qa ## Run all checks (alias for qa)

fix: format ## Fix code issues (alias for format)

# CI/CD simulation
ci: ## Simulate CI/CD pipeline
	@echo "🔄 Simulating CI/CD pipeline..."
	make clean
	make install-dev
	make lint
	make type-check
	make security
	make test-fast
	make build
	@echo "✅ CI/CD simulation complete!"

# Performance and profiling
profile: ## Run performance profiling
	@echo "📈 Performance profiling not yet implemented"
	@echo "💡 TODO: Add profiling with py-spy or similar"

benchmark: ## Run benchmarks
	@echo "⚡ Running performance benchmarks..."
	$(PYTHON) -c "
import time
from crios.scoring.design import crowe_score
test_smiles = ['CCO', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C']
start = time.time()
for smiles in test_smiles * 100:
    result = crowe_score(smiles)
elapsed = time.time() - start
print(f'Scored {len(test_smiles) * 100} molecules in {elapsed:.2f}s')
print(f'Average: {(elapsed / (len(test_smiles) * 100)) * 1000:.1f}ms per molecule')
"

# Docker targets
docker-build: ## Build Docker image
	docker build -t crios:latest -f docker/Dockerfile .

docker-run: ## Run Docker container
	docker run -it --rm -p 8000:8000 crios:latest

docker-dev: ## Run Docker container in development mode
	docker-compose -f docker/docker-compose.yml up --build

# Release targets
release-patch: ## Bump patch version and tag
	@echo "🏷️ Creating patch release..."
	@echo "💡 Version bumping not yet implemented"

release-minor: ## Bump minor version and tag
	@echo "🏷️ Creating minor release..."
	@echo "💡 Version bumping not yet implemented"

release-major: ## Bump major version and tag
	@echo "🏷️ Creating major release..."
	@echo "💡 Version bumping not yet implemented"

# Environment setup
setup: ## Initial setup for new contributors
	@echo "🚀 Setting up CriOS development environment..."
	@echo ""
	@echo "📋 Prerequisites:"
	@echo "  - Python 3.11+"
	@echo "  - Git"
	@echo "  - Make"
	@echo ""
	make install-dev
	@echo ""
	@echo "✅ Setup complete! Try 'make demo' to test the installation."

# Utility targets
info: ## Show project information
	@echo "🧬 CriOS Discovery Engine"
	@echo "========================"
	@echo "Version: $$($(PYTHON) -c 'import crios; print(crios.__version__)' 2>/dev/null || echo 'Not installed')"
	@echo "Python: $$($(PYTHON) --version)"
	@echo "Platform: $$(uname -s)"
	@echo "Working Directory: $$(pwd)"
	@echo "Git Branch: $$(git branch --show-current 2>/dev/null || echo 'Not a git repository')"
	@echo "Git Commit: $$(git rev-parse --short HEAD 2>/dev/null || echo 'Not a git repository')"

# Help with common tasks
recipes: ## Show common development recipes
	@echo "🧬 CriOS Development Recipes"
	@echo "============================"
	@echo ""
	@echo "🚀 First time setup:"
	@echo "  make setup"
	@echo ""
	@echo "🔧 Development workflow:"
	@echo "  make dev          # Set up development environment"
	@echo "  make demo         # Run demonstration"
	@echo "  make test-fast    # Run fast tests during development"
	@echo "  make format       # Format code"
	@echo "  make qa           # Run all quality checks"
	@echo ""
	@echo "🏗️ Before committing:"
	@echo "  make qa           # Run all checks"
	@echo "  make test         # Run full test suite"
	@echo ""
	@echo "🌐 Web development:"
	@echo "  make web          # Start development server"
	@echo "  make web-prod     # Start production server"
	@echo ""
	@echo "📦 Building and distribution:"
	@echo "  make build        # Build distribution packages"
	@echo "  make ci           # Simulate CI/CD pipeline"