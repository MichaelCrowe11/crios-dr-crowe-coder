from __future__ import annotations
from typing import Dict, Tuple
import math
import re

# Built-in rule bounds: (min, max) inclusive where not None
_BUILTIN = {
    "lipinski": {"MW": (None, 500), "LogP": (None, 5), "HBD": (None, 5), "HBA": (None, 10)},
    "lead_like": {"MW": (200, 350), "LogP": (-1, 3), "RotatableBonds": (None, 7)},
}

class MolecularFilter:
    """Lipinski/Lead-like + simple rule parser: e.g., 'MW<500 AND LogP<5 AND HBD<=5'."""
    def __init__(self, name_or_expr: str):
        self.name_or_expr = name_or_expr.strip()

    def _check_bounds(self, desc: Dict[str, float], bounds: Dict[str, Tuple[float|None,float|None]]) -> Tuple[bool, Dict[str,str]]:
        violations: Dict[str,str] = {}
        for k, (lo, hi) in bounds.items():
            v = desc.get(k, float("nan"))
            if math.isnan(v):
                violations[k] = "NA"
                continue
            if lo is not None and v < lo: violations[k] = f"{v:.2f} < {lo}"
            if hi is not None and v > hi: violations[k] = f"{v:.2f} > {hi}"
        return (len(violations) == 0), violations

    def evaluate(self, mol) -> bool:
        # Ensure required descriptors exist
        from .descriptors import DescriptorCalculator
        needed = set(["MW","LogP","HBD","HBA","RotatableBonds"])
        desc = mol.calculate_descriptors(list(needed))

        if self.name_or_expr.lower() in _BUILTIN:
            ok, _ = self._check_bounds(desc, _BUILTIN[self.name_or_expr.lower()])
            return ok

        # Parse simple boolean expression (AND/OR) of comparisons
        expr = self.name_or_expr.upper()
        tokens = re.split(r"\s+(AND|OR)\s+", expr)
        def eval_comp(comp: str) -> bool:
            m = re.fullmatch(r"\s*([A-Z0-9_]+)\s*(<=|>=|<|>|=|==)\s*([-\d\.]+)\s*", comp, re.I)
            if not m: return False
            key, op, rhs = m.group(1), m.group(2), float(m.group(3))
            val = desc.get(key)
            if val is None or math.isnan(val): return False
            if op in ("=", "=="): return float(val) == rhs
            if op == "<": return val < rhs
            if op == "<=": return val <= rhs
            if op == ">": return val > rhs
            if op == ">=": return val >= rhs
            return False
        # Left-to-right with AND/OR (no parentheses)
        result = eval_comp(tokens[0])
        i = 1
        while i < len(tokens):
            op, comp = tokens[i], tokens[i+1]
            val = eval_comp(comp)
            result = (result and val) if op == "AND" else (result or val)
            i += 2
        return bool(result)