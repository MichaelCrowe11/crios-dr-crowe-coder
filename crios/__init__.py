"""
CriOS Discovery Engine
Crowe Discovery Framework - Universal Compound Discovery Platform

"Science before status. Discovery before profit. Earth and people above extraction."
"""

__version__ = "0.1.0"
__author__ = "Dr. Michael Crowe"
__email__ = "michael@crowelogic.com"

from .core import *
from .config import ConfigManager

__all__ = [
    "ConfigManager",
    "Molecule", 
    "CompoundLibrary",
    "SimilaritySearch",
    "ClusterAnalyzer"
]