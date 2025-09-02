"""
CriOS Discovery Engine Exceptions
Custom exception classes for error handling
"""


class CriosError(Exception):
    """Base exception for all CriOS errors"""
    pass


class ValidationError(CriosError):
    """Raised when data validation fails"""
    pass


class MoleculeError(CriosError):
    """Raised when molecule processing fails"""
    pass


class CompoundError(CriosError):
    """Raised when compound operations fail"""
    pass


class SimilarityError(CriosError):
    """Raised when similarity calculations fail"""
    pass


class ClusteringError(CriosError):
    """Raised when clustering operations fail"""
    pass


class DatabaseError(CriosError):
    """Raised when database operations fail"""
    pass


class APIError(CriosError):
    """Raised when external API calls fail"""
    pass


class EthicsViolationError(CriosError):
    """Raised when ethical policies are violated"""
    pass


class ConfigurationError(CriosError):
    """Raised when configuration is invalid"""
    pass


class CalculationError(CriosError):
    """Raised when mathematical calculations fail"""
    pass


class FilterError(CriosError):
    """Raised when filtering operations fail"""
    pass


class GeneratorError(CriosError):
    """Raised when compound generation fails"""
    pass