"""
CriOS Ethics and Safety Filters
Structural alerts, controlled substance detection, and safety screening
"""

import logging
import hashlib
from typing import List, Dict, Optional, Any, Tuple
import yaml
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from ..core.exceptions import CriosError
from ..chem.mol import Molecule

logger = logging.getLogger(__name__)


class EthicsError(CriosError):
    """Ethics and safety related errors"""
    pass


class SafetyViolation:
    """Represents a safety policy violation"""
    
    def __init__(
        self,
        rule_id: str,
        rule_name: str,
        severity: str,
        description: str,
        triggered_by: str,
        confidence: float = 1.0
    ):
        self.rule_id = rule_id
        self.rule_name = rule_name
        self.severity = severity
        self.description = description
        self.triggered_by = triggered_by
        self.confidence = confidence
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "rule_id": self.rule_id,
            "rule_name": self.rule_name,
            "severity": self.severity,
            "description": self.description,
            "triggered_by": self.triggered_by,
            "confidence": self.confidence
        }


class StructuralAlertFilter:
    """Filter for structural alerts and toxicophores"""
    
    def __init__(self, alerts: Optional[Dict[str, Any]] = None):
        """
        Initialize structural alert filter
        
        Args:
            alerts: Dictionary of alert definitions
        """
        self.alerts = alerts or self._get_default_alerts()
    
    def _get_default_alerts(self) -> Dict[str, Any]:
        """Get default structural alerts (examples only)"""
        return {
            "nitro_aromatic": {
                "smarts": "[$(c1ccccc1)][N+](=O)[O-]",
                "severity": "medium",
                "description": "Aromatic nitro groups - potential mutagenicity concern"
            },
            "aromatic_amine": {
                "smarts": "c[NH2]",
                "severity": "low",
                "description": "Aromatic amine - monitor for potential carcinogenicity"
            },
            "epoxide": {
                "smarts": "C1OC1",
                "severity": "high",
                "description": "Epoxide group - potential DNA reactivity"
            },
            "aldehyde": {
                "smarts": "[CH]=O",
                "severity": "low",
                "description": "Aldehyde group - potential protein reactivity"
            },
            "acyl_halide": {
                "smarts": "C(=O)[F,Cl,Br,I]",
                "severity": "high",
                "description": "Acyl halide - highly reactive"
            },
            "isocyanate": {
                "smarts": "[N-]=[C+]=[O]",
                "severity": "high",
                "description": "Isocyanate - respiratory sensitizer"
            }
        }
    
    def check_molecule(self, smiles: str) -> List[SafetyViolation]:
        """
        Check molecule against structural alerts
        
        Args:
            smiles: SMILES string
        
        Returns:
            List of violations
        """
        violations = []
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return violations
            
            for alert_id, alert_info in self.alerts.items():
                smarts = alert_info.get("smarts")
                if not smarts:
                    continue
                
                try:
                    pattern = Chem.MolFromSmarts(smarts)
                    if pattern is None:
                        logger.warning(f"Invalid SMARTS pattern: {smarts}")
                        continue
                    
                    matches = mol.GetSubstructMatches(pattern)
                    if matches:
                        violation = SafetyViolation(
                            rule_id=f"structural_alert_{alert_id}",
                            rule_name=alert_info.get("name", alert_id),
                            severity=alert_info.get("severity", "medium"),
                            description=alert_info.get("description", f"Structural alert: {alert_id}"),
                            triggered_by=f"Substructure matches: {len(matches)}"
                        )
                        violations.append(violation)
                        
                except Exception as e:
                    logger.warning(f"Error checking alert {alert_id}: {e}")
                    continue
            
            return violations
            
        except Exception as e:
            logger.error(f"Structural alert check failed for {smiles}: {e}")
            return violations


class ControlledSubstanceFilter:
    """Filter for controlled substances and illicit compounds"""
    
    def __init__(self, blocklists: Optional[Dict[str, Any]] = None):
        """
        Initialize controlled substance filter
        
        Args:
            blocklists: Dictionary of blocklist definitions
        """
        self.blocked_hashes = blocklists.get("blocked_hashes", []) if blocklists else []
        self.blocked_patterns = blocklists.get("blocked_patterns", []) if blocklists else []
        self.similarity_blocklist = blocklists.get("similarity_blocklist", {}) if blocklists else {}
    
    def check_molecule(self, smiles: str) -> List[SafetyViolation]:
        """
        Check molecule against controlled substance filters
        
        Args:
            smiles: SMILES string
        
        Returns:
            List of violations
        """
        violations = []
        
        try:
            # Hash-based checking
            canonical_smiles = self._canonicalize_smiles(smiles)
            if canonical_smiles:
                smiles_hash = hashlib.sha256(canonical_smiles.encode()).hexdigest()
                
                if smiles_hash in self.blocked_hashes:
                    violation = SafetyViolation(
                        rule_id="controlled_substance_hash",
                        rule_name="Controlled Substance (Hash Match)",
                        severity="critical",
                        description="Molecule matches known controlled substance",
                        triggered_by=f"Hash: {smiles_hash[:16]}..."
                    )
                    violations.append(violation)
            
            # Pattern-based checking
            violations.extend(self._check_blocked_patterns(smiles))
            
            # Similarity-based checking
            violations.extend(self._check_similarity_blocklist(smiles))
            
            return violations
            
        except Exception as e:
            logger.error(f"Controlled substance check failed for {smiles}: {e}")
            return violations
    
    def _canonicalize_smiles(self, smiles: str) -> Optional[str]:
        """Canonicalize SMILES string"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            return Chem.MolToSmiles(mol)
        except:
            return None
    
    def _check_blocked_patterns(self, smiles: str) -> List[SafetyViolation]:
        """Check against blocked SMARTS patterns"""
        violations = []
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return violations
            
            for pattern in self.blocked_patterns:
                try:
                    pattern_mol = Chem.MolFromSmarts(pattern)
                    if pattern_mol is None:
                        continue
                    
                    if mol.HasSubstructMatch(pattern_mol):
                        violation = SafetyViolation(
                            rule_id="controlled_substance_pattern",
                            rule_name="Controlled Substance (Pattern Match)",
                            severity="critical",
                            description="Molecule contains blocked substructure",
                            triggered_by=f"Pattern: {pattern}"
                        )
                        violations.append(violation)
                        
                except Exception as e:
                    logger.warning(f"Error checking pattern {pattern}: {e}")
                    continue
            
            return violations
            
        except Exception as e:
            logger.warning(f"Pattern checking failed: {e}")
            return violations
    
    def _check_similarity_blocklist(self, smiles: str) -> List[SafetyViolation]:
        """Check against similarity blocklist"""
        violations = []
        
        # This would require implementing similarity comparison
        # against reference structures. For now, return empty list
        # as no reference structures are provided by default.
        
        return violations


class ContextualGuardrailFilter:
    """Filter for contextual safety guardrails"""
    
    def __init__(self, guardrails: Optional[Dict[str, Any]] = None):
        """
        Initialize contextual guardrail filter
        
        Args:
            guardrails: Dictionary of guardrail definitions
        """
        self.guardrails = guardrails or self._get_default_guardrails()
    
    def _get_default_guardrails(self) -> Dict[str, Any]:
        """Get default guardrails"""
        return {
            "blocked_keywords": [
                "synthesis of controlled substances",
                "drug manufacturing",
                "illegal drug production",
                "weapon development",
                "bioweapon",
                "chemical weapon"
            ],
            "forbidden_classes": [
                "explosive precursors",
                "nerve agents",
                "biological toxins"
            ]
        }
    
    def check_context(self, context: str) -> List[SafetyViolation]:
        """
        Check context (e.g., user prompt) against guardrails
        
        Args:
            context: Context string to check
        
        Returns:
            List of violations
        """
        violations = []
        
        context_lower = context.lower()
        
        # Check blocked keywords
        for keyword in self.guardrails.get("blocked_keywords", []):
            if keyword.lower() in context_lower:
                violation = SafetyViolation(
                    rule_id="blocked_keyword",
                    rule_name="Blocked Keyword",
                    severity="critical",
                    description=f"Context contains blocked keyword: {keyword}",
                    triggered_by=keyword
                )
                violations.append(violation)
        
        # Check forbidden classes
        for forbidden_class in self.guardrails.get("forbidden_classes", []):
            if forbidden_class.lower() in context_lower:
                violation = SafetyViolation(
                    rule_id="forbidden_class",
                    rule_name="Forbidden Class",
                    severity="critical",
                    description=f"Context mentions forbidden class: {forbidden_class}",
                    triggered_by=forbidden_class
                )
                violations.append(violation)
        
        return violations


class EthicsFilter:
    """Comprehensive ethics and safety filter"""
    
    def __init__(self, policy_file: Optional[Path] = None):
        """
        Initialize ethics filter
        
        Args:
            policy_file: Path to ethics policy YAML file
        """
        self.policy = self._load_policy(policy_file)
        
        # Initialize sub-filters
        self.structural_filter = StructuralAlertFilter(
            self.policy.get("structural_alerts", {})
        )
        self.controlled_filter = ControlledSubstanceFilter(
            self.policy.get("controlled_substances", {})
        )
        self.guardrail_filter = ContextualGuardrailFilter(
            self.policy.get("guardrails", {})
        )
    
    def _load_policy(self, policy_file: Optional[Path]) -> Dict[str, Any]:
        """Load ethics policy from file"""
        if policy_file and policy_file.exists():
            try:
                with open(policy_file, 'r') as f:
                    return yaml.safe_load(f)
            except Exception as e:
                logger.warning(f"Failed to load policy file {policy_file}: {e}")
        
        # Return minimal default policy
        return {
            "metadata": {
                "version": "1.0.0",
                "policy_name": "Default CriOS Safety Policy"
            },
            "enforcement": {
                "mode": "strict",
                "auto_reject_violations": True
            }
        }
    
    def check_molecule(self, smiles: str) -> Tuple[bool, List[SafetyViolation]]:
        """
        Comprehensive ethics check for a molecule
        
        Args:
            smiles: SMILES string
        
        Returns:
            Tuple of (passed, violations)
        """
        all_violations = []
        
        # Structural alerts
        structural_violations = self.structural_filter.check_molecule(smiles)
        all_violations.extend(structural_violations)
        
        # Controlled substances
        controlled_violations = self.controlled_filter.check_molecule(smiles)
        all_violations.extend(controlled_violations)
        
        # Determine if molecule passes
        enforcement_mode = self.policy.get("enforcement", {}).get("mode", "strict")
        
        if enforcement_mode == "strict":
            # Any critical violation fails
            critical_violations = [v for v in all_violations if v.severity == "critical"]
            passed = len(critical_violations) == 0
        elif enforcement_mode == "warning":
            # Only block if explicitly configured
            passed = True
        else:  # permissive
            passed = True
        
        return passed, all_violations
    
    def check_context(self, context: str) -> Tuple[bool, List[SafetyViolation]]:
        """
        Check context against guardrails
        
        Args:
            context: Context string
        
        Returns:
            Tuple of (passed, violations)
        """
        violations = self.guardrail_filter.check_context(context)
        
        # Context violations are always blocking in strict mode
        enforcement_mode = self.policy.get("enforcement", {}).get("mode", "strict")
        passed = len(violations) == 0 or enforcement_mode == "permissive"
        
        return passed, violations
    
    def explain_violations(self, violations: List[SafetyViolation]) -> Dict[str, Any]:
        """
        Generate explanation for violations
        
        Args:
            violations: List of violations
        
        Returns:
            Explanation dictionary
        """
        if not violations:
            return {
                "total_violations": 0,
                "severity_breakdown": {},
                "explanations": []
            }
        
        # Severity breakdown
        severity_counts = {}
        for violation in violations:
            severity_counts[violation.severity] = severity_counts.get(violation.severity, 0) + 1
        
        # Detailed explanations
        explanations = [violation.to_dict() for violation in violations]
        
        return {
            "total_violations": len(violations),
            "severity_breakdown": severity_counts,
            "explanations": explanations,
            "recommendation": self._get_recommendation(violations)
        }
    
    def _get_recommendation(self, violations: List[SafetyViolation]) -> str:
        """Get recommendation based on violations"""
        if not violations:
            return "No safety concerns identified."
        
        critical_count = len([v for v in violations if v.severity == "critical"])
        high_count = len([v for v in violations if v.severity == "high"])
        
        if critical_count > 0:
            return "REJECT: Critical safety violations detected. This compound should not be pursued."
        elif high_count > 0:
            return "CAUTION: High-severity safety concerns. Extensive safety evaluation required."
        else:
            return "PROCEED WITH CAUTION: Minor safety alerts detected. Consider structural modifications."


def check_ethics_compliance(
    smiles: str,
    policy_file: Optional[Path] = None,
    context: Optional[str] = None
) -> Dict[str, Any]:
    """
    Convenience function for ethics compliance checking
    
    Args:
        smiles: SMILES string
        policy_file: Path to ethics policy file
        context: Optional context string
    
    Returns:
        Ethics compliance result
    """
    ethics_filter = EthicsFilter(policy_file)
    
    # Check molecule
    mol_passed, mol_violations = ethics_filter.check_molecule(smiles)
    
    # Check context if provided
    context_passed = True
    context_violations = []
    if context:
        context_passed, context_violations = ethics_filter.check_context(context)
    
    # Combined result
    all_violations = mol_violations + context_violations
    overall_passed = mol_passed and context_passed
    
    return {
        "smiles": smiles,
        "passed": overall_passed,
        "molecule_passed": mol_passed,
        "context_passed": context_passed,
        "violations": [v.to_dict() for v in all_violations],
        "explanation": ethics_filter.explain_violations(all_violations)
    }