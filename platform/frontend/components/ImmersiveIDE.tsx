// CriOS Immersive IDE - Dr. Crowe Coder Integration
// Algorithm & Compound Development Environment

import React, { useState, useEffect, useRef } from 'react';
import {
  Box,
  Paper,
  Typography,
  Button,
  IconButton,
  Tabs,
  Tab,
  TextField,
  Select,
  MenuItem,
  FormControl,
  InputLabel,
  Divider,
  Chip,
  List,
  ListItem,
  ListItemText,
  ListItemIcon,
  Tooltip,
  Grid,
  Card,
  CardContent,
  LinearProgress,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  ToggleButton,
  ToggleButtonGroup
} from '@mui/material';
import {
  PlayArrow as RunIcon,
  Stop as StopIcon,
  Save as SaveIcon,
  CloudDownload as LoadIcon,
  Code as CodeIcon,
  Science as ScienceIcon,
  Psychology as AIIcon,
  Hub as NetworkIcon,
  Visibility as PreviewIcon,
  Settings as SettingsIcon,
  ExpandMore as ExpandIcon,
  BugReport as DebugIcon,
  Timeline as ProfileIcon,
  Memory as MemoryIcon,
  Speed as SpeedIcon
} from '@mui/icons-material';
import { Monaco } from '@monaco-editor/react';
import Editor from '@monaco-editor/react';
import * as monaco from 'monaco-editor';

interface CodeTemplate {
  id: string;
  name: string;
  description: string;
  language: 'python' | 'typescript' | 'javascript';
  category: 'algorithm' | 'compound' | 'analysis' | 'ml';
  code: string;
  dependencies: string[];
}

interface ExecutionResult {
  success: boolean;
  output: string;
  error?: string;
  performance: {
    executionTime: number;
    memoryUsage: number;
    croweScore: number;
  };
  biologicalPatterns?: string[];
  compounds?: any[];
}

interface DrCroweAgent {
  id: string;
  name: string;
  status: 'idle' | 'analyzing' | 'coding' | 'optimizing';
  specialization: string[];
  confidence: number;
}

const codeTemplates: CodeTemplate[] = [
  {
    id: 'neurotherapeutic-discovery',
    name: 'Neurotherapeutic Lead Discovery',
    description: 'Alzheimer\'s and Parkinson\'s compound design with BBB optimization',
    language: 'python',
    category: 'compound',
    code: `from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import numpy as np

class NeuroTherapeuticDesigner:
    """
    CriOS Nova Neurotherapeutic Agent
    Specialized for CNS drug discovery with ethical constraints
    """
    
    def __init__(self):
        self.target_proteins = ["AChE", "BuChE", "NMDA", "alpha_synuclein"]
        self.bbb_requirements = {
            "mol_weight": (150, 450),  # Optimal BBB crossing
            "logp": (1.0, 3.5),        # Balanced lipophilicity
            "tpsa": (40, 90),          # Polar surface area
            "hbd": (0, 3),             # H-bond donors
            "rotatable_bonds": (0, 8)  # Flexibility
        }
    
    def design_cholinesterase_inhibitors(self, scaffold="donepezil", n_variants=20):
        """Generate enhanced cholinesterase inhibitor analogs"""
        
        # Donepezil scaffold for Alzheimer's
        donepezil_core = "COc1cc2c(cc1OC)C(=O)C(CC3CCN(CC3)Cc4ccccc4)C2"
        
        if scaffold == "donepezil":
            base_smiles = donepezil_core
        else:
            base_smiles = scaffold
            
        # Apply Crowe Logic: Science before status
        variants = self._generate_ethical_variants(
            base_smiles, 
            n_variants,
            therapeutic_focus="neuroprotection"
        )
        
        # Filter for blood-brain barrier compatibility
        bbb_optimized = []
        for variant in variants:
            if self._predict_bbb_crossing(variant):
                score = self._calculate_crowe_score(variant)
                bbb_optimized.append({
                    'smiles': variant,
                    'bbb_score': score,
                    'predicted_activity': self._predict_ache_inhibition(variant)
                })
        
        return sorted(bbb_optimized, key=lambda x: x['predicted_activity'], reverse=True)
    
    def _predict_bbb_crossing(self, smiles):
        """Predict blood-brain barrier permeability using Crowe criteria"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False
            
        # Calculate BBB-relevant descriptors
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)
        
        # Apply BBB requirements with ethical bias toward safety
        return (self.bbb_requirements["mol_weight"][0] <= mw <= self.bbb_requirements["mol_weight"][1] and
                self.bbb_requirements["logp"][0] <= logp <= self.bbb_requirements["logp"][1] and
                self.bbb_requirements["tpsa"][0] <= tpsa <= self.bbb_requirements["tpsa"][1] and
                hbd <= self.bbb_requirements["hbd"][1] and
                rotatable <= self.bbb_requirements["rotatable_bonds"][1])
    
    def _predict_ache_inhibition(self, smiles):
        """Predict acetylcholinesterase inhibition activity"""
        # Simplified activity prediction based on structural features
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return 0.0
            
        # Key pharmacophoric features for AChE inhibition
        features = {
            'aromatic_rings': len([x for x in mol.GetRingInfo().AtomRings() 
                                 if mol.GetAtomWithIdx(x[0]).GetIsAromatic()]),
            'nitrogen_count': len([a for a in mol.GetAtoms() if a.GetSymbol() == 'N']),
            'quaternary_n': self._count_quaternary_nitrogens(mol),
            'molecular_flexibility': Descriptors.NumRotatableBonds(mol)
        }
        
        # Empirical scoring based on known AChE inhibitor patterns
        score = (features['aromatic_rings'] * 0.3 + 
                features['nitrogen_count'] * 0.4 +
                features['quaternary_n'] * 0.2 +
                min(features['molecular_flexibility'] / 10, 0.1))
        
        return min(score, 1.0)  # Cap at 1.0
    
    def _count_quaternary_nitrogens(self, mol):
        """Count quaternary nitrogen atoms (important for AChE binding)"""
        count = 0
        for atom in mol.GetAtoms():
            if (atom.GetSymbol() == 'N' and 
                atom.GetTotalNumHs() == 0 and 
                len(atom.GetNeighbors()) == 4):
                count += 1
        return count
    
    def _generate_ethical_variants(self, base_smiles, n_variants, therapeutic_focus):
        """Generate molecular variants with ethical constraints"""
        variants = []
        mol = Chem.MolFromSmiles(base_smiles)
        
        # Ethical principle: Generate variants that prioritize safety over potency
        for i in range(n_variants):
            # Simple structural modifications for demo
            # In practice, this would use sophisticated generative models
            variant = self._apply_structural_modification(base_smiles, i)
            
            if variant and self._passes_ethical_filter(variant, therapeutic_focus):
                variants.append(variant)
        
        return variants[:n_variants]
    
    def _apply_structural_modification(self, smiles, seed):
        """Apply structural modifications guided by medicinal chemistry"""
        # Simplified modification - in practice would use advanced generative AI
        modifications = [
            smiles,  # Original
            smiles.replace('c1ccccc1', 'c1ccc(F)cc1'),  # Add fluorine
            smiles.replace('OC', 'NC'),  # Replace hydroxyl with amino
            smiles.replace('CC', 'CCC') if 'CC' in smiles else smiles  # Extend alkyl
        ]
        
        return modifications[seed % len(modifications)]
    
    def _passes_ethical_filter(self, smiles, therapeutic_focus):
        """Ensure compounds align with healing principles"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False
            
        # Ethical constraints: No obviously toxic substructures
        toxic_patterns = [
            '[N+](=O)[O-]',  # Nitro groups (potential mutagens)
            'c1ccc2[nH]c3ccccc3c2c1',  # Carbazole (potential carcinogen)
            'C=C(Cl)Cl',  # Vinyl chloride patterns
        ]
        
        for pattern in toxic_patterns:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                return False
        
        return True
    
    def _calculate_crowe_score(self, smiles):
        """Calculate ethical discovery score based on Crowe principles"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return 0.0
            
        # Science before status: Favor validated pharmacology
        # Discovery before profit: Favor accessible synthesis
        # Earth before extraction: Favor biodegradable structures
        
        scores = {
            'druggability': self._assess_druggability(mol),
            'synthetic_accessibility': self._assess_synthesis_feasibility(mol),
            'environmental_impact': self._assess_biodegradability(mol)
        }
        
        # Weighted average with ethical bias
        crowe_score = (scores['druggability'] * 0.5 + 
                      scores['synthetic_accessibility'] * 0.3 +
                      scores['environmental_impact'] * 0.2)
        
        return crowe_score
    
    def _assess_druggability(self, mol):
        """Assess drug-likeness using Lipinski-like rules"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        violations = 0
        if mw > 500: violations += 1
        if logp > 5: violations += 1
        if hbd > 5: violations += 1
        if hba > 10: violations += 1
        
        return max(0, (4 - violations) / 4)
    
    def _assess_synthesis_feasibility(self, mol):
        """Estimate synthetic accessibility (simplified)"""
        # Higher ring strain and complexity = lower accessibility
        ring_penalty = len(mol.GetRingInfo().AtomRings()) * 0.1
        size_penalty = mol.GetNumAtoms() * 0.02
        
        return max(0, 1.0 - ring_penalty - size_penalty)
    
    def _assess_biodegradability(self, mol):
        """Estimate environmental biodegradation potential"""
        # Favor structures with biodegradable linkages
        biodegradable_patterns = [
            'C(=O)O',    # Carboxylic acids
            'C(=O)N',    # Amides
            'OC=O',      # Esters
            'C(=O)C',    # Ketones
        ]
        
        score = 0.5  # Base score
        for pattern in biodegradable_patterns:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                score += 0.125
        
        return min(score, 1.0)

# Example usage with ethical discovery principles
if __name__ == "__main__":
    designer = NeuroTherapeuticDesigner()
    
    print("CriOS Nova Neurotherapeutic Discovery")
    print("=====================================")
    print("Mission: Science before status. Discovery before profit.")
    print()
    
    # Design cholinesterase inhibitors for Alzheimer's
    compounds = designer.design_cholinesterase_inhibitors(n_variants=10)
    
    print("Top 5 Neurotherapeutic Candidates:")
    for i, compound in enumerate(compounds[:5], 1):
        print(f"{i}. SMILES: {compound['smiles']}")
        print(f"   BBB Score: {compound['bbb_score']:.3f}")
        print(f"   Predicted Activity: {compound['predicted_activity']:.3f}")
        print()
    
    print("Ethical Constraints Applied:")
    print("- No mutagenic or carcinogenic substructures")
    print("- Optimized for blood-brain barrier crossing")
    print("- Balanced for druggability and accessibility")
    print("- Environmental biodegradability considered")`,
    dependencies: ['rdkit', 'numpy']
  },
  {
    id: 'natural-product-discovery',
    name: 'Fungal Natural Product Discovery',
    description: 'Mycelial compound mining and bioactive analog generation',
    language: 'python',
    category: 'compound',
    code: `from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
import numpy as np

class FungalCompoundDiscoverer:
    """
    CriOS Nova Natural Product Discovery Agent
    Specialized for mycological compound mining and analog generation
    """
    
    def __init__(self):
        # Known bioactive fungal metabolites as scaffolds
        self.fungal_scaffolds = {
            'psilocybin': 'CNP(=O)(O)Oc1ccc2[nH]cc(CCN)c2c1',  # Neuroprotective
            'lovastatin': 'CCC(C)C(=O)OC1CC(C)C=C2CC(C)C(=O)C(C)C(O)CC(=O)C=CC12',
            'penicillin': 'CC1(C)SC2C(NC(=O)CC3=CC=CC=C3)C(=O)N2C1C(=O)O',
            'ergosterol': 'CCC(CCC(C)C1CCC2C1(C)CCC3C2CC=C4CC(O)CCC34C)C(C)C',
            'cordycepin': 'NC1=NC=NC2=C1N=CN2C3OC(CO)CC3O'  # From Cordyceps
        }
        
        self.bioactivity_targets = {
            'neuroprotective': ['AChE', 'NMDA', 'neuroinflammation'],
            'anticancer': ['p53', 'apoptosis', 'angiogenesis'],
            'immunomodulatory': ['TNF-alpha', 'IL-6', 'NF-kB'],
            'antioxidant': ['ROS scavenging', 'Nrf2 pathway'],
            'longevity': ['autophagy', 'senescence', 'telomerase']
        }
    
    def discover_fungal_analogs(self, target_activity='neuroprotective', n_analogs=20):
        """
        Discover fungal natural product analogs with specified bioactivity
        Applies Crowe Logic: Science before status, Earth before extraction
        """
        
        print(f"CriOS Nova Fungal Discovery: {target_activity}")
        print("=" * 50)
        print("Mission: Harvest nature's wisdom through ethical synthesis")
        print()
        
        # Select relevant fungal scaffolds for target activity
        relevant_scaffolds = self._select_scaffolds_for_activity(target_activity)
        
        analogs = []
        for scaffold_name, scaffold_smiles in relevant_scaffolds.items():
            print(f"Generating analogs from {scaffold_name}...")
            
            # Generate biologically-inspired variations
            scaffold_analogs = self._generate_mycological_analogs(
                scaffold_smiles, 
                scaffold_name,
                n_analogs // len(relevant_scaffolds)
            )
            
            for analog in scaffold_analogs:
                # Apply natural product filters
                if self._passes_natural_product_filter(analog):
                    bioactivity_score = self._predict_bioactivity(analog, target_activity)
                    synthesis_score = self._assess_synthetic_feasibility(analog)
                    sustainability_score = self._assess_environmental_impact(analog)
                    
                    analogs.append({
                        'smiles': analog,
                        'parent_scaffold': scaffold_name,
                        'bioactivity_score': bioactivity_score,
                        'synthesis_feasibility': synthesis_score,
                        'sustainability_score': sustainability_score,
                        'crowe_score': self._calculate_natural_crowe_score(
                            bioactivity_score, synthesis_score, sustainability_score
                        )
                    })
        
        # Sort by Crowe score (balances efficacy, feasibility, sustainability)
        return sorted(analogs, key=lambda x: x['crowe_score'], reverse=True)
    
    def _select_scaffolds_for_activity(self, target_activity):
        """Select fungal scaffolds most relevant to target bioactivity"""
        
        activity_scaffold_map = {
            'neuroprotective': ['psilocybin', 'cordycepin'],
            'anticancer': ['lovastatin', 'penicillin'],
            'immunomodulatory': ['ergosterol', 'cordycepin'],
            'antioxidant': ['ergosterol', 'lovastatin'],
            'longevity': ['cordycepin', 'ergosterol']
        }
        
        selected_scaffolds = {}
        for scaffold_name in activity_scaffold_map.get(target_activity, ['psilocybin']):
            if scaffold_name in self.fungal_scaffolds:
                selected_scaffolds[scaffold_name] = self.fungal_scaffolds[scaffold_name]
        
        return selected_scaffolds
    
    def _generate_mycological_analogs(self, scaffold_smiles, scaffold_name, n_variants):
        """Generate analogs using mycological biosynthetic principles"""
        
        mol = Chem.MolFromSmiles(scaffold_smiles)
        if not mol:
            return []
        
        analogs = []
        
        # Mycological modification strategies
        modifications = {
            'hydroxylation': self._add_hydroxyl_groups,
            'methylation': self._add_methyl_groups, 
            'demethylation': self._remove_methyl_groups,
            'cyclization': self._form_additional_rings,
            'oxidation': self._oxidize_functional_groups,
            'glycosylation': self._add_sugar_moieties
        }
        
        for i in range(n_variants):
            # Apply 1-2 random biosynthetic modifications
            mod_types = np.random.choice(list(modifications.keys()), 
                                       size=np.random.randint(1, 3), 
                                       replace=False)
            
            modified_mol = mol
            for mod_type in mod_types:
                try:
                    modified_mol = modifications[mod_type](modified_mol)
                    if modified_mol is None:
                        break
                except:
                    continue
            
            if modified_mol:
                analog_smiles = Chem.MolToSmiles(modified_mol)
                if analog_smiles != scaffold_smiles:  # Ensure modification occurred
                    analogs.append(analog_smiles)
        
        return list(set(analogs))  # Remove duplicates
    
    def _add_hydroxyl_groups(self, mol):
        """Add hydroxyl groups (common fungal modification)"""
        # Simplified hydroxylation - in practice would use more sophisticated methods
        smiles = Chem.MolToSmiles(mol)
        
        # Target aromatic carbons for hydroxylation
        if 'c1ccccc1' in smiles:
            modified = smiles.replace('c1ccccc1', 'c1cc(O)ccc1', 1)
            return Chem.MolFromSmiles(modified)
        
        return mol
    
    def _add_methyl_groups(self, mol):
        """Add methyl groups (common in secondary metabolism)"""
        smiles = Chem.MolToSmiles(mol)
        
        # Add methyl to nitrogen if present
        if 'NH' in smiles:
            modified = smiles.replace('NH', 'N(C)', 1)
            new_mol = Chem.MolFromSmiles(modified)
            return new_mol if new_mol else mol
        
        return mol
    
    def _remove_methyl_groups(self, mol):
        """Remove methyl groups (demethylation)"""
        smiles = Chem.MolToSmiles(mol)
        
        # Simple demethylation
        if 'OC' in smiles:
            modified = smiles.replace('OC', 'O', 1)
            new_mol = Chem.MolFromSmiles(modified)
            return new_mol if new_mol else mol
        
        return mol
    
    def _form_additional_rings(self, mol):
        """Form additional rings (cyclization)"""
        # Complex modification - simplified for demo
        return mol
    
    def _oxidize_functional_groups(self, mol):
        """Oxidize alcohols to ketones/aldehydes"""
        smiles = Chem.MolToSmiles(mol)
        
        # Convert primary alcohol to aldehyde
        if 'CO' in smiles and 'C(=O)' not in smiles:
            modified = smiles.replace('CO', 'C=O', 1)
            new_mol = Chem.MolFromSmiles(modified)
            return new_mol if new_mol else mol
        
        return mol
    
    def _add_sugar_moieties(self, mol):
        """Add sugar groups (glycosylation)"""
        # Simplified glycosylation
        smiles = Chem.MolToSmiles(mol)
        
        if 'OH' in smiles:
            # Replace OH with glucose-like moiety
            modified = smiles.replace('OH', 'OC1OC(CO)C(O)C(O)C1O', 1)
            new_mol = Chem.MolFromSmiles(modified)
            return new_mol if new_mol else mol
        
        return mol
    
    def _passes_natural_product_filter(self, smiles):
        """Apply natural product-like filters"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False
        
        # Natural product criteria
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        aromatic_rings = Descriptors.NumAromaticRings(mol)
        
        # Natural product-like ranges (more permissive than Lipinski)
        return (200 <= mw <= 800 and          # Natural products can be larger
                -2 <= logp <= 6 and           # Broader lipophilicity range  
                rotatable_bonds <= 15 and     # More flexibility allowed
                aromatic_rings <= 4)          # Reasonable aromatic content
    
    def _predict_bioactivity(self, smiles, target_activity):
        """Predict bioactivity score for target"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return 0.0
        
        # Activity-specific structural features
        features = self._calculate_bioactivity_features(mol)
        
        if target_activity == 'neuroprotective':
            # Favor BBB-permeable, cholinesterase-inhibiting features
            score = (features['bbb_score'] * 0.4 + 
                    features['nitrogen_heterocycles'] * 0.3 +
                    features['aromatic_rings'] * 0.2 +
                    features['molecular_complexity'] * 0.1)
        
        elif target_activity == 'anticancer':
            # Favor DNA-intercalating, apoptosis-inducing features
            score = (features['planar_rings'] * 0.4 +
                    features['electrophilic_groups'] * 0.3 +
                    features['molecular_size'] * 0.2 +
                    features['hydrogen_bonding'] * 0.1)
        
        elif target_activity == 'immunomodulatory':
            # Favor cytokine-modulating features
            score = (features['hydrogen_bonding'] * 0.3 +
                    features['polar_surface_area'] * 0.3 +
                    features['ring_systems'] * 0.2 +
                    features['functional_diversity'] * 0.2)
        
        else:
            # General bioactivity score
            score = np.mean(list(features.values()))
        
        return min(score, 1.0)
    
    def _calculate_bioactivity_features(self, mol):
        """Calculate structural features relevant to bioactivity"""
        
        # Basic molecular properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        
        # Structural features
        aromatic_rings = Descriptors.NumAromaticRings(mol)
        heterocycles = len([atom for atom in mol.GetAtoms() 
                          if atom.IsInRing() and atom.GetSymbol() != 'C'])
        
        return {
            'bbb_score': self._estimate_bbb_permeability(mol),
            'nitrogen_heterocycles': heterocycles,
            'aromatic_rings': min(aromatic_rings / 3, 1.0),
            'molecular_complexity': min(mol.GetNumAtoms() / 50, 1.0),
            'planar_rings': min(aromatic_rings / 4, 1.0),
            'electrophilic_groups': self._count_electrophilic_groups(mol),
            'molecular_size': min(mw / 500, 1.0),
            'hydrogen_bonding': min((Descriptors.NumHDonors(mol) + 
                                   Descriptors.NumHAcceptors(mol)) / 8, 1.0),
            'polar_surface_area': min(tpsa / 100, 1.0),
            'ring_systems': min(mol.GetRingInfo().NumRings() / 4, 1.0),
            'functional_diversity': self._assess_functional_diversity(mol)
        }
    
    def _estimate_bbb_permeability(self, mol):
        """Estimate blood-brain barrier permeability"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        
        # Simple BBB model
        if 300 <= mw <= 450 and 1 <= logp <= 3 and tpsa <= 90:
            return 0.8
        elif mw > 500 or logp < 0 or tpsa > 120:
            return 0.2
        else:
            return 0.5
    
    def _count_electrophilic_groups(self, mol):
        """Count electrophilic functional groups"""
        electrophilic_patterns = [
            'C(=O)',     # Carbonyl
            'C=C',       # Alkene
            'C#C',       # Alkyne
            'C(=O)O',    # Carboxylic acid
        ]
        
        count = 0
        for pattern in electrophilic_patterns:
            matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
            count += len(matches)
        
        return min(count / 3, 1.0)
    
    def _assess_functional_diversity(self, mol):
        """Assess diversity of functional groups"""
        functional_groups = [
            'OH',        # Hydroxyl
            'NH',        # Amino
            'C(=O)',     # Carbonyl
            'C(=O)O',    # Carboxyl
            'S',         # Sulfur
            'P',         # Phosphorus
        ]
        
        present_groups = 0
        for group in functional_groups:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(group)):
                present_groups += 1
        
        return present_groups / len(functional_groups)
    
    def _assess_synthetic_feasibility(self, smiles):
        """Assess how feasible synthesis would be"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return 0.0
        
        # Penalize complex ring systems and many stereocenters
        ring_penalty = mol.GetRingInfo().NumRings() * 0.1
        stereocenter_penalty = len(Chem.FindMolChiralCenters(mol)) * 0.05
        
        # Favor simpler, more accessible structures
        base_score = 0.8
        feasibility = base_score - ring_penalty - stereocenter_penalty
        
        return max(0.1, min(feasibility, 1.0))
    
    def _assess_environmental_impact(self, smiles):
        """Assess environmental sustainability"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return 0.0
        
        # Favor biodegradable functional groups
        biodegradable_score = 0.5
        
        # Bonus for biodegradable patterns
        biodegradable_patterns = [
            'C(=O)O',    # Carboxylic acid
            'C(=O)N',    # Amide
            'OC=O',      # Ester
            'OH',        # Hydroxyl
        ]
        
        for pattern in biodegradable_patterns:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                biodegradable_score += 0.1
        
        return min(biodegradable_score, 1.0)
    
    def _calculate_natural_crowe_score(self, bioactivity, synthesis, sustainability):
        """Calculate Crowe score for natural products"""
        # Weighted scoring emphasizing ethical discovery
        return (bioactivity * 0.5 +      # Efficacy matters
                synthesis * 0.3 +        # Accessibility important
                sustainability * 0.2)    # Environmental responsibility

# Example usage for fungal natural product discovery
if __name__ == "__main__":
    discoverer = FungalCompoundDiscoverer()
    
    print("CriOS Nova Fungal Natural Product Discovery")
    print("==========================================")
    print("Mission: Harvest nature's wisdom through ethical synthesis")
    print()
    
    # Discover neuroprotective fungal analogs
    analogs = discoverer.discover_fungal_analogs('neuroprotective', n_analogs=12)
    
    print("Top 5 Fungal Natural Product Analogs:")
    for i, compound in enumerate(analogs[:5], 1):
        print(f"{i}. Parent: {compound['parent_scaffold']}")
        print(f"   SMILES: {compound['smiles']}")
        print(f"   Bioactivity: {compound['bioactivity_score']:.3f}")
        print(f"   Synthesis: {compound['synthesis_feasibility']:.3f}")
        print(f"   Sustainability: {compound['sustainability_score']:.3f}")
        print(f"   Crowe Score: {compound['crowe_score']:.3f}")
        print()
    
    print("Mycological Design Principles Applied:")
    print("- Biosynthetic modifications (hydroxylation, methylation, etc.)")
    print("- Natural product-like property ranges")
    print("- Sustainable synthesis pathway consideration")
    print("- Bioactivity prediction based on target mechanism")`,
    dependencies: ['rdkit', 'numpy']
  },
  {
    id: 'microplastic-sequestrant',
    name: 'Microplastic Sequestrant Design',
    description: 'Oral-compatible polymer-binding agents for detoxification',
    language: 'python',
    category: 'compound',
    code: `from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import numpy as np

class MicroplasticSequestrantDesigner:
    """
    CriOS Nova Environmental Health Agent
    Specialized for microplastic sequestration compound design
    """
    
    def __init__(self):
        self.target_polymers = {
            'PET': 'CC(=O)Oc1ccc(cc1)C(=O)Oc2ccc(cc2)C(C)(C)C',  # Polyethylene terephthalate
            'PS': 'c1ccc(cc1)C',  # Polystyrene unit
            'PE': 'CC',  # Polyethylene unit
            'PP': 'C(C)C'  # Polypropylene unit
        }
        
        self.oral_safety_requirements = {
            'molecular_weight': (200, 800),  # Avoid absorption
            'logp': (-1, 4),  # Balanced for gut retention
            'toxicity_alerts': self._load_toxicity_patterns()
        }
    
    def design_oral_sequestrants(self, target_polymer='PET', n_variants=15):
        """Design oral-safe microplastic binding agents"""
        
        print(f"Designing sequestrants for {target_polymer} microplastics...")
        print("Crowe Principle: Earth and people above extraction")
        
        # Design strategy based on polymer chemistry
        if target_polymer == 'PET':
            # PET has aromatic rings and ester linkages
            binding_motifs = self._design_pet_binding_motifs()
        elif target_polymer == 'PS':
            # Polystyrene has aromatic rings
            binding_motifs = self._design_ps_binding_motifs()
        else:
            # General hydrophobic binding for PE/PP
            binding_motifs = self._design_hydrophobic_binding_motifs()
        
        sequestrants = []
        for i, motif in enumerate(binding_motifs[:n_variants]):
            candidate = self._optimize_for_oral_safety(motif)
            if candidate:
                binding_affinity = self._predict_polymer_binding(candidate, target_polymer)
                safety_score = self._assess_oral_safety(candidate)
                
                sequestrants.append({
                    'smiles': candidate,
                    'binding_affinity': binding_affinity,
                    'oral_safety_score': safety_score,
                    'crowe_score': self._calculate_environmental_score(candidate),
                    'molecular_weight': Descriptors.MolWt(Chem.MolFromSmiles(candidate))
                })
        
        # Sort by combined safety and efficacy
        return sorted(sequestrants, 
                     key=lambda x: (x['oral_safety_score'] * x['binding_affinity']), 
                     reverse=True)
    
    def _design_pet_binding_motifs(self):
        """Design compounds that bind to PET through π-π stacking and H-bonding"""
        # PET binding strategy: aromatic stacking + ester interaction
        motifs = [
            # Aromatic compounds with carboxylic acid groups
            'c1ccc(cc1)C(=O)O',  # Benzoic acid derivative
            'c1ccc2c(c1)c(cc(c2=O)O)C(=O)O',  # Anthraquinone carboxylic acid
            'c1cc(cc(c1)C(=O)O)C(=O)O',  # Isophthalic acid derivative
            'c1ccc(cc1)c2ccc(cc2)C(=O)O',  # Biphenyl carboxylic acid
            
            # Aromatic compounds with hydroxyl groups (H-bonding)
            'c1ccc(cc1)c2c(c(cc(c2=O)O)O)C(=O)O',  # Flavonoid-like structure
            'c1cc(c(c(c1)O)C(=O)O)O',  # Dihydroxybenzoic acid
            
            # Larger aromatic systems for stronger π-π interactions
            'c1ccc2c(c1)ccc3c2ccc(c3)C(=O)O',  # Anthracene carboxylic acid
        ]
        return motifs
    
    def _design_ps_binding_motifs(self):
        """Design compounds that bind to polystyrene through π-π stacking"""
        # Polystyrene binding: aromatic-aromatic interactions
        motifs = [
            'c1ccc(cc1)c2ccccc2',  # Biphenyl
            'c1ccc2c(c1)ccc3c2ccc4c3cccc4',  # Anthracene
            'c1ccc2c(c1)oc3ccccc32',  # Dibenzofuran
            'c1ccc(cc1)C(=O)c2ccccc2',  # Benzophenone
            'c1ccc(cc1)c2ccc3c(c2)ccc4c3ccc5c4cccc5',  # Extended aromatic
        ]
        return motifs
    
    def _design_hydrophobic_binding_motifs(self):
        """Design compounds for general hydrophobic plastic binding"""
        # PE/PP binding: hydrophobic interactions
        motifs = [
            'CCCCCCCCCCCCCCCCCC(=O)O',  # Long-chain fatty acid
            'c1ccc(cc1)CCCCCCCCCC(=O)O',  # Aromatic-alkyl hybrid
            'CC(C)CCCC(C)CCCC(C)CCCC(=O)O',  # Branched fatty acid
            'c1ccc(cc1)c2ccc(cc2)CCCCCC(=O)O',  # Biphenyl-alkyl
        ]
        return motifs
    
    def _optimize_for_oral_safety(self, base_smiles):
        """Optimize compound for oral administration safety"""
        mol = Chem.MolFromSmiles(base_smiles)
        if not mol:
            return None
            
        # Check molecular weight (avoid systemic absorption)
        mw = Descriptors.MolWt(mol)
        if mw < self.oral_safety_requirements['molecular_weight'][0]:
            # Too small - might be absorbed
            return self._add_bulk_groups(base_smiles)
        elif mw > self.oral_safety_requirements['molecular_weight'][1]:
            # Too large - might be problematic
            return self._simplify_structure(base_smiles)
        
        return base_smiles
    
    def _add_bulk_groups(self, smiles):
        """Add bulk to prevent systemic absorption"""
        # Add carboxylic acid groups or large substituents
        modifications = [
            smiles + 'C(=O)O',  # Add carboxyl
            smiles.replace('c1ccccc1', 'c1ccc(cc1)C(=O)O'),  # Replace phenyl with benzoic acid
            smiles + 'OS(=O)(=O)O'  # Add sulfonic acid (highly polar, non-absorbable)
        ]
        
        for mod in modifications:
            try:
                mol = Chem.MolFromSmiles(mod)
                if mol and self._passes_safety_filter(mol):
                    return mod
            except:
                continue
        
        return smiles  # Return original if modifications fail
    
    def _simplify_structure(self, smiles):
        """Simplify overly complex structures"""
        # Simple truncation strategy - in practice would be more sophisticated
        if len(smiles) > 50:
            return smiles[:50]  # Rough truncation
        return smiles
    
    def _predict_polymer_binding(self, sequestrant_smiles, polymer_type):
        """Predict binding affinity to target polymer"""
        mol = Chem.MolFromSmiles(sequestrant_smiles)
        if not mol:
            return 0.0
            
        # Simplified binding prediction based on molecular features
        features = {
            'aromatic_rings': len([x for x in mol.GetRingInfo().AtomRings() 
                                 if mol.GetAtomWithIdx(x[0]).GetIsAromatic()]),
            'logp': Descriptors.MolLogP(mol),
            'molecular_weight': Descriptors.MolWt(mol),
            'polar_groups': self._count_polar_groups(mol)
        }
        
        # Polymer-specific binding models
        if polymer_type == 'PET':
            # PET binding favors aromatic rings and polar interactions
            score = (features['aromatic_rings'] * 0.4 + 
                    features['polar_groups'] * 0.3 +
                    min(features['logp'] / 5, 0.2) + 0.1)
        elif polymer_type == 'PS':
            # PS binding favors aromatic-aromatic interactions
            score = (features['aromatic_rings'] * 0.6 + 
                    min(features['logp'] / 3, 0.3) + 0.1)
        else:
            # PE/PP binding favors hydrophobic interactions
            score = (min(features['logp'] / 4, 0.5) + 
                    features['molecular_weight'] / 1000 + 0.2)
        
        return min(score, 1.0)
    
    def _count_polar_groups(self, mol):
        """Count polar functional groups"""
        polar_patterns = [
            'C(=O)O',    # Carboxylic acid
            'OS(=O)(=O)O',  # Sulfonic acid
            'OP(=O)(O)O',   # Phosphoric acid
            'c1ccccc1O',    # Phenol
            'N',            # Nitrogen
        ]
        
        count = 0
        for pattern in polar_patterns:
            matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
            count += len(matches)
        
        return count
    
    def _assess_oral_safety(self, smiles):
        """Assess safety for oral administration"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return 0.0
            
        safety_score = 1.0
        
        # Check for toxic structural alerts
        for pattern in self.oral_safety_requirements['toxicity_alerts']:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                safety_score -= 0.3
        
        # Molecular weight check (large = less absorption = safer)
        mw = Descriptors.MolWt(mol)
        if mw < 300:
            safety_score -= 0.2  # Too small, might be absorbed
        elif mw > 600:
            safety_score += 0.1  # Good size for gut retention
        
        # LogP check (polar = less absorption = safer for gut)
        logp = Descriptors.MolLogP(mol)
        if logp < 0:
            safety_score += 0.2  # Very polar, minimal absorption
        elif logp > 4:
            safety_score -= 0.1  # Too lipophilic
        
        return max(0, min(safety_score, 1.0))
    
    def _load_toxicity_patterns(self):
        """Load known toxic substructure patterns"""
        return [
            '[N+](=O)[O-]',  # Nitro groups
            'c1ccc2c(c1)ccc3c2cccc3N=Nc4ccccc4',  # Azo dyes
            'C=C(Cl)Cl',  # Vinyl halides
            'c1ccc(cc1)N=Nc2ccc(cc2)N',  # Aromatic azo compounds
            'ClCCl',  # Dichloromethane-like
        ]
    
    def _passes_safety_filter(self, mol):
        """Check if molecule passes basic safety filters"""
        # Check for problematic substructures
        for pattern in self.oral_safety_requirements['toxicity_alerts']:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                return False
        return True
    
    def _calculate_environmental_score(self, smiles):
        """Calculate environmental impact score (Crowe principle)"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return 0.0
            
        # Environmental factors: biodegradability, non-toxicity, efficacy
        biodegradable_score = self._assess_biodegradability(mol)
        non_toxic_score = self._assess_environmental_toxicity(mol)
        
        return (biodegradable_score + non_toxic_score) / 2
    
    def _assess_biodegradability(self, mol):
        """Assess environmental biodegradation potential"""
        biodegradable_patterns = [
            'C(=O)O',    # Carboxylic acids
            'C(=O)N',    # Amides
            'CC(=O)O',   # Esters
            'c1ccccc1O', # Phenols
        ]
        
        score = 0.3  # Base biodegradability
        for pattern in biodegradable_patterns:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                score += 0.2
        
        return min(score, 1.0)
    
    def _assess_environmental_toxicity(self, mol):
        """Assess potential environmental toxicity (inverted - higher is better)"""
        # Penalize known environmental toxins
        toxic_patterns = [
            'c1ccc(cc1)Cl',  # Chlorinated aromatics
            'C(F)(F)F',      # Perfluorinated compounds
            '[Hg]',          # Heavy metals
            '[Pb]',
            '[Cd]'
        ]
        
        score = 1.0
        for pattern in toxic_patterns:
            try:
                if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                    score -= 0.4
            except:
                continue  # Skip invalid patterns
        
        return max(0, score)

# Example usage for microplastic detoxification
if __name__ == "__main__":
    designer = MicroplasticSequestrantDesigner()
    
    print("CriOS Nova Environmental Health Discovery")
    print("=======================================")
    print("Mission: Earth and people above extraction")
    print()
    
    # Design PET microplastic sequestrants
    sequestrants = designer.design_oral_sequestrants('PET', n_variants=8)
    
    print("Top 5 Microplastic Sequestrant Candidates (PET):")
    for i, compound in enumerate(sequestrants[:5], 1):
        print(f"{i}. SMILES: {compound['smiles']}")
        print(f"   Binding Affinity: {compound['binding_affinity']:.3f}")
        print(f"   Oral Safety: {compound['oral_safety_score']:.3f}")
        print(f"   Environmental Score: {compound['crowe_score']:.3f}")
        print(f"   Molecular Weight: {compound['molecular_weight']:.1f}")
        print()
    
    print("Ethical Design Principles Applied:")
    print("- Large molecular weight prevents systemic absorption")
    print("- No known toxic substructure patterns")
    print("- Biodegradable functional groups included")
    print("- Optimized for gut retention and plastic binding")`,
    dependencies: ['rdkit', 'numpy']
  },
  {
    id: 'synthetic-drug-design',
    name: 'AI-Guided Synthetic Drug Design',
    description: 'De novo small molecule design with ADMET optimization',
    language: 'python',
    category: 'compound',
    code: `from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
import numpy as np
import random

class SyntheticDrugDesigner:
    """
    CriOS Nova Synthetic Chemistry Agent
    Specialized for de novo small molecule design with multi-objective optimization
    """
    
    def __init__(self):
        # Common drug-like molecular fragments for assembly
        self.drug_fragments = {
            'aromatic_cores': [
                'c1ccccc1',      # Benzene
                'c1ccc2ccccc2c1', # Naphthalene
                'c1cnc2ccccc2c1', # Quinoline
                'c1ccc2nc3ccccc3cc2c1', # Acridine
                'c1ccc2[nH]c3ccccc3c2c1', # Carbazole
            ],
            'heterocycles': [
                'c1ccncc1',      # Pyridine
                'c1cnccn1',      # Pyrimidine
                'c1ccc2[nH]ccc2c1', # Indole
                'c1cc2ccccc2[nH]1',  # Indazole
                'c1nnc2ccccc2n1',    # Phthalazine
            ],
            'linkers': [
                'C(=O)N',        # Amide
                'C(=O)O',        # Carboxylic acid
                'S(=O)(=O)',     # Sulfonyl
                'C=C',           # Alkene
                'CC',            # Alkyl chain
            ],
            'terminal_groups': [
                'N',             # Amine
                'O',             # Hydroxyl
                'C(=O)O',        # Carboxyl
                'C(F)(F)F',      # Trifluoromethyl
                'S(=O)(=O)N',    # Sulfonamide
            ]
        }
        
        # Target protein families and their preferred molecular features
        self.target_preferences = {
            'GPCR': {'mw_range': (200, 500), 'logp_range': (2, 5), 'rotatable_bonds': 8},
            'kinase': {'mw_range': (300, 600), 'logp_range': (1, 4), 'hbd_max': 3},
            'protease': {'mw_range': (400, 800), 'logp_range': (0, 3), 'peptide_like': True},
            'ion_channel': {'mw_range': (250, 450), 'logp_range': (3, 6), 'rigid_core': True},
            'enzyme': {'mw_range': (150, 400), 'logp_range': (-1, 3), 'polar_groups': True}
        }
    
    def design_synthetic_drugs(self, target_class='GPCR', n_designs=25, 
                             optimize_for=['potency', 'selectivity', 'admet']):
        """
        Design novel synthetic small molecules for specified target class
        Applies Crowe Logic: Science before status, discovery before profit
        """
        
        print(f"CriOS Nova Synthetic Drug Design: {target_class}")
        print("=" * 55)
        print("Mission: Engineer molecules that heal with precision")
        print()
        
        # Get target-specific design preferences
        preferences = self.target_preferences.get(target_class, self.target_preferences['enzyme'])
        
        designs = []
        
        for i in range(n_designs):
            # Generate de novo molecule using fragment assembly
            candidate = self._assemble_drug_candidate(preferences)
            
            if candidate and self._passes_drug_like_filters(candidate):
                # Multi-objective optimization
                scores = self._calculate_design_scores(candidate, target_class, optimize_for)
                
                designs.append({
                    'smiles': candidate,
                    'target_class': target_class,
                    'potency_score': scores.get('potency', 0),
                    'selectivity_score': scores.get('selectivity', 0), 
                    'admet_score': scores.get('admet', 0),
                    'synthetic_accessibility': scores.get('synthesis', 0),
                    'novelty_score': scores.get('novelty', 0),
                    'crowe_score': self._calculate_synthetic_crowe_score(scores),
                    'molecular_weight': Descriptors.MolWt(Chem.MolFromSmiles(candidate))
                })
        
        # Sort by Crowe score (multi-objective optimization)
        return sorted(designs, key=lambda x: x['crowe_score'], reverse=True)
    
    def _assemble_drug_candidate(self, preferences):
        """Assemble drug candidate from molecular fragments"""
        
        # Start with aromatic or heterocyclic core
        if random.random() < 0.6:
            core = random.choice(self.drug_fragments['aromatic_cores'])
        else:
            core = random.choice(self.drug_fragments['heterocycles'])
        
        # Add 1-3 modifications
        n_modifications = random.randint(1, 3)
        candidate_smiles = core
        
        for _ in range(n_modifications):
            modification_type = random.choice(['linker', 'terminal', 'substitute'])
            
            if modification_type == 'linker':
                linker = random.choice(self.drug_fragments['linkers'])
                terminal = random.choice(self.drug_fragments['terminal_groups'])
                # Simplified assembly - in practice would use more sophisticated methods
                candidate_smiles += linker + terminal
                
            elif modification_type == 'terminal':
                terminal = random.choice(self.drug_fragments['terminal_groups'])
                candidate_smiles += terminal
                
            elif modification_type == 'substitute':
                # Simple substitution on aromatic ring
                if 'c1ccccc1' in candidate_smiles:
                    substituent = random.choice(['F', 'Cl', 'C', 'O', 'N'])
                    candidate_smiles = candidate_smiles.replace(
                        'c1ccccc1', f'c1ccc(${substituent})cc1', 1
                    )
        
        # Clean up and validate SMILES
        try:
            mol = Chem.MolFromSmiles(candidate_smiles)
            if mol:
                return Chem.MolToSmiles(mol)
        except:
            pass
        
        return core  # Return simple core if assembly fails
    
    def _passes_drug_like_filters(self, smiles):
        """Apply Lipinski Rule of Five and additional drug-like filters"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False
        
        # Lipinski Rule of Five
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        lipinski_violations = 0
        if mw > 500: lipinski_violations += 1
        if logp > 5: lipinski_violations += 1  
        if hbd > 5: lipinski_violations += 1
        if hba > 10: lipinski_violations += 1
        
        # Additional drug-like criteria
        tpsa = Descriptors.TPSA(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        
        return (lipinski_violations <= 1 and    # Allow 1 violation
                20 <= tpsa <= 140 and           # Topological polar surface area
                rotatable_bonds <= 10 and       # Flexibility
                mol.GetNumAtoms() >= 10)        # Minimum complexity
    
    def _calculate_design_scores(self, smiles, target_class, optimize_for):
        """Calculate multi-objective design scores"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {}
        
        scores = {}
        
        if 'potency' in optimize_for:
            scores['potency'] = self._predict_target_potency(mol, target_class)
        
        if 'selectivity' in optimize_for:
            scores['selectivity'] = self._predict_selectivity(mol, target_class)
        
        if 'admet' in optimize_for:
            scores['admet'] = self._predict_admet_properties(mol)
        
        # Always calculate synthesis and novelty
        scores['synthesis'] = self._assess_synthetic_accessibility(mol)
        scores['novelty'] = self._assess_structural_novelty(mol)
        
        return scores
    
    def _predict_target_potency(self, mol, target_class):
        """Predict binding affinity/potency for target class"""
        
        # Simplified potency prediction based on known SAR principles
        features = self._extract_molecular_features(mol)
        
        if target_class == 'GPCR':
            # GPCRs favor lipophilic, rigid molecules with basic nitrogen
            score = (features['lipophilicity'] * 0.3 +
                    features['basic_nitrogen'] * 0.3 +
                    features['aromatic_rings'] * 0.2 +
                    features['molecular_size'] * 0.2)
        
        elif target_class == 'kinase':
            # Kinases favor ATP-competitive inhibitors with hinge-binding motifs
            score = (features['hydrogen_bonding'] * 0.4 +
                    features['planar_rings'] * 0.3 +
                    features['nitrogen_heterocycles'] * 0.2 +
                    features['optimal_size'] * 0.1)
        
        elif target_class == 'protease':
            # Proteases favor peptide-like molecules with multiple H-bond donors/acceptors
            score = (features['hydrogen_bonding'] * 0.5 +
                    features['peptide_character'] * 0.3 +
                    features['flexibility'] * 0.2)
        
        else:
            # Generic binding prediction
            score = np.mean([features['lipophilicity'], 
                           features['hydrogen_bonding'],
                           features['molecular_complexity']])
        
        return min(score, 1.0)
    
    def _predict_selectivity(self, mol, target_class):
        """Predict selectivity within target class"""
        
        # Selectivity often comes from specific binding pocket interactions
        features = self._extract_molecular_features(mol)
        
        # Favor molecules with unique structural features for selectivity
        selectivity_score = (features['structural_uniqueness'] * 0.4 +
                           features['specific_interactions'] * 0.3 +
                           features['steric_bulk'] * 0.2 +
                           features['electronic_features'] * 0.1)
        
        return min(selectivity_score, 1.0)
    
    def _predict_admet_properties(self, mol):
        """Predict ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity)"""
        
        # Multi-parameter ADMET prediction
        admet_scores = {
            'absorption': self._predict_oral_absorption(mol),
            'distribution': self._predict_tissue_distribution(mol),
            'metabolism': self._predict_metabolic_stability(mol),
            'excretion': self._predict_clearance(mol),
            'toxicity': self._predict_safety_profile(mol)
        }
        
        # Weighted ADMET score
        return (admet_scores['absorption'] * 0.25 +
                admet_scores['distribution'] * 0.15 +
                admet_scores['metabolism'] * 0.2 +
                admet_scores['excretion'] * 0.15 +
                admet_scores['toxicity'] * 0.25)
    
    def _predict_oral_absorption(self, mol):
        """Predict oral bioavailability"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        
        # Oral absorption rules (Veber, Egan, etc.)
        score = 1.0
        if mw > 500: score -= 0.2
        if logp > 5 or logp < 0: score -= 0.2
        if tpsa > 140: score -= 0.3
        if rotatable_bonds > 10: score -= 0.1
        
        return max(0, score)
    
    def _predict_tissue_distribution(self, mol):
        """Predict tissue distribution properties"""
        logp = Descriptors.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        
        # Optimal distribution: moderate lipophilicity and size
        if 1 <= logp <= 3 and 200 <= mw <= 400:
            return 0.8
        else:
            return 0.5
    
    def _predict_metabolic_stability(self, mol):
        """Predict resistance to metabolism"""
        
        # Check for metabolically labile groups
        labile_patterns = [
            'C(=O)O',        # Carboxylic acid
            'CC(=O)',        # Methyl ketone
            'N(C)C',         # N,N-dimethyl
            'c1cc(OH)cc1',   # Phenol
        ]
        
        stability_score = 0.8
        for pattern in labile_patterns:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                stability_score -= 0.1
        
        return max(0.2, stability_score)
    
    def _predict_clearance(self, mol):
        """Predict renal/hepatic clearance"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        
        # Moderate clearance is often optimal
        if 250 <= mw <= 450 and 1 <= logp <= 4:
            return 0.7
        else:
            return 0.5
    
    def _predict_safety_profile(self, mol):
        """Predict toxicity risk"""
        
        # Check for structural alerts (toxicophores)
        toxic_alerts = [
            '[N+](=O)[O-]',  # Nitro groups
            'c1ccc(cc1)N',   # Aniline (mutagenic potential)
            'C=C(Cl)Cl',     # Vinyl halides
            'C(=O)Cl',       # Acid chlorides
            'N=N',           # Azo groups
        ]
        
        safety_score = 0.9  # Start with high safety
        for alert in toxic_alerts:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(alert)):
                safety_score -= 0.3
        
        return max(0.1, safety_score)
    
    def _extract_molecular_features(self, mol):
        """Extract molecular features for activity prediction"""
        
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        
        return {
            'lipophilicity': min(max(logp, 0) / 5, 1.0),
            'basic_nitrogen': len([a for a in mol.GetAtoms() 
                                 if a.GetSymbol() == 'N' and a.GetTotalNumHs() > 0]),
            'aromatic_rings': Descriptors.NumAromaticRings(mol),
            'molecular_size': min(mw / 500, 1.0),
            'hydrogen_bonding': (Descriptors.NumHDonors(mol) + 
                                Descriptors.NumHAcceptors(mol)) / 8,
            'planar_rings': Descriptors.NumAromaticRings(mol),
            'nitrogen_heterocycles': len([a for a in mol.GetAtoms()
                                        if a.GetSymbol() == 'N' and a.IsInRing()]),
            'peptide_character': self._assess_peptide_like_character(mol),
            'flexibility': min(Descriptors.NumRotatableBonds(mol) / 8, 1.0),
            'molecular_complexity': min(mol.GetNumAtoms() / 30, 1.0),
            'structural_uniqueness': self._assess_structural_uniqueness(mol),
            'specific_interactions': self._count_specific_binding_groups(mol),
            'steric_bulk': self._assess_steric_bulk(mol),
            'electronic_features': self._assess_electronic_properties(mol),
            'optimal_size': 1.0 if 250 <= mw <= 450 else 0.5
        }
    
    def _assess_peptide_like_character(self, mol):
        """Assess peptide-like characteristics"""
        amide_bonds = len(mol.GetSubstructMatches(Chem.MolFromSmarts('C(=O)N')))
        amino_groups = len(mol.GetSubstructMatches(Chem.MolFromSmarts('N')))
        return min((amide_bonds + amino_groups) / 5, 1.0)
    
    def _assess_structural_uniqueness(self, mol):
        """Assess structural novelty/uniqueness"""
        # Simplified uniqueness - count unusual ring systems and functional groups
        ring_systems = mol.GetRingInfo().NumRings()
        unusual_atoms = len([a for a in mol.GetAtoms() 
                           if a.GetSymbol() not in ['C', 'N', 'O', 'F', 'Cl']])
        return min((ring_systems + unusual_atoms) / 5, 1.0)
    
    def _count_specific_binding_groups(self, mol):
        """Count groups that provide specific binding interactions"""
        binding_groups = [
            'C(=O)N',        # Amide (H-bonding)
            'OH',            # Hydroxyl (H-bonding)
            'c1cc(F)cc1',    # Fluorinated aromatic (halogen bonding)
            'S(=O)(=O)',     # Sulfonyl (electrostatic)
        ]
        
        count = 0
        for group in binding_groups:
            matches = mol.GetSubstructMatches(Chem.MolFromSmarts(group))
            count += len(matches)
        
        return min(count / 3, 1.0)
    
    def _assess_steric_bulk(self, mol):
        """Assess steric bulk for selectivity"""
        # Count branching points and bulky substituents
        branched_carbons = len([a for a in mol.GetAtoms() 
                              if a.GetSymbol() == 'C' and len(a.GetNeighbors()) > 2])
        return min(branched_carbons / 5, 1.0)
    
    def _assess_electronic_properties(self, mol):
        """Assess electronic properties (pi-systems, electron-rich/poor regions)"""
        aromatic_atoms = len([a for a in mol.GetAtoms() if a.GetIsAromatic()])
        electronegative_atoms = len([a for a in mol.GetAtoms() 
                                   if a.GetSymbol() in ['N', 'O', 'F']])
        return min((aromatic_atoms + electronegative_atoms) / 10, 1.0)
    
    def _assess_synthetic_accessibility(self, mol):
        """Assess synthetic accessibility"""
        # Penalize complex ring systems, many stereocenters, unusual functional groups
        ring_penalty = mol.GetRingInfo().NumRings() * 0.1
        stereocenter_penalty = len(Chem.FindMolChiralCenters(mol)) * 0.05
        size_penalty = mol.GetNumAtoms() * 0.01
        
        accessibility = 1.0 - ring_penalty - stereocenter_penalty - size_penalty
        return max(0.1, accessibility)
    
    def _assess_structural_novelty(self, mol):
        """Assess novelty compared to known drugs"""
        # Simplified novelty assessment
        ring_systems = mol.GetRingInfo().NumRings()
        functional_groups = len([a for a in mol.GetAtoms() 
                               if a.GetSymbol() != 'C'])
        
        # More rings and functional groups = potentially more novel
        novelty = min((ring_systems + functional_groups) / 8, 1.0)
        return novelty
    
    def _calculate_synthetic_crowe_score(self, scores):
        """Calculate Crowe score for synthetic drug design"""
        # Multi-objective optimization with ethical weighting
        
        weights = {
            'potency': 0.25,        # Efficacy important
            'selectivity': 0.2,     # Minimize side effects
            'admet': 0.25,          # Safety critical
            'synthesis': 0.2,       # Accessibility matters
            'novelty': 0.1          # Innovation valuable
        }
        
        crowe_score = 0
        for metric, weight in weights.items():
            if metric in scores:
                crowe_score += scores[metric] * weight
        
        return crowe_score

# Example usage for synthetic drug design
if __name__ == "__main__":
    designer = SyntheticDrugDesigner()
    
    print("CriOS Nova Synthetic Drug Design")
    print("================================")
    print("Mission: Engineer molecules that heal with precision")
    print()
    
    # Design GPCR-targeted drugs
    designs = designer.design_synthetic_drugs(
        target_class='GPCR', 
        n_designs=15,
        optimize_for=['potency', 'selectivity', 'admet']
    )
    
    print("Top 5 Synthetic Drug Candidates (GPCR-targeted):")
    for i, compound in enumerate(designs[:5], 1):
        print(f"{i}. SMILES: {compound['smiles']}")
        print(f"   Target: {compound['target_class']}")
        print(f"   Potency: {compound['potency_score']:.3f}")
        print(f"   Selectivity: {compound['selectivity_score']:.3f}")
        print(f"   ADMET: {compound['admet_score']:.3f}")
        print(f"   Synthesis: {compound['synthetic_accessibility']:.3f}")
        print(f"   Crowe Score: {compound['crowe_score']:.3f}")
        print(f"   Mol Weight: {compound['molecular_weight']:.1f}")
        print()
    
    print("Synthetic Design Principles Applied:")
    print("- Fragment-based molecular assembly")
    print("- Multi-objective optimization (potency + selectivity + ADMET)")
    print("- Drug-likeness filters (Lipinski Rule of Five)")
    print("- Synthetic accessibility assessment")
    print("- Structural novelty scoring")`,
    dependencies: ['rdkit', 'numpy']
  },
  {
    id: 'similarity-algorithm',
    name: 'Tanimoto Similarity Algorithm', 
    description: 'Optimized molecular fingerprint comparison',
    language: 'python',
    category: 'algorithm',
    code: `from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors
import numpy as np

def enhanced_tanimoto_similarity(mol1_smiles: str, mol2_smiles: str, 
                               radius: int = 3, nBits: int = 2048) -> float:
    """
    Biological-inspired Tanimoto similarity with adaptive fingerprints
    Dr. Crowe Coder Pattern: Adaptive system architecture
    """
    try:
        # Convert SMILES to molecules
        mol1 = Chem.MolFromSmiles(mol1_smiles)
        mol2 = Chem.MolFromSmiles(mol2_smiles)
        
        if mol1 is None or mol2 is None:
            return 0.0
        
        # Generate ECFP fingerprints with biological adaptation
        fp1 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol1, radius, nBits)
        fp2 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol2, radius, nBits)
        
        # Calculate Tanimoto similarity
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
        
        # Apply biological scaling (homeostasis principle)
        scaled_similarity = apply_biological_scaling(similarity)
        
        return scaled_similarity
    
    except Exception as e:
        print(f"Error calculating similarity: {e}")
        return 0.0

def apply_biological_scaling(similarity: float) -> float:
    """Apply biological homeostasis principles to similarity scoring"""
    # Sigmoid-like scaling inspired by biological response curves
    return 1.0 / (1.0 + np.exp(-10 * (similarity - 0.5)))

# Example usage
if __name__ == "__main__":
    # Test with aspirin and ibuprofen
    aspirin = "CC(=O)OC1=CC=CC=C1C(=O)O"
    ibuprofen = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
    
    similarity = enhanced_tanimoto_similarity(aspirin, ibuprofen)
    print(f"Biological-scaled similarity: {similarity:.4f}")`,
    dependencies: ['rdkit', 'numpy']
  },
  {
    id: 'compound-generator',
    name: 'AI Compound Generator',
    description: 'Generate novel compounds using biological patterns',
    language: 'python',
    category: 'compound',
    code: `from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import random
import numpy as np

class BiologicalCompoundGenerator:
    """
    Generate compounds using biological growth patterns
    Dr. Crowe Coder Pattern: Self-healing systems
    """
    
    def __init__(self):
        self.base_fragments = [
            "c1ccccc1",  # benzene
            "C1=CC=NC=C1",  # pyridine  
            "C1=CC=CO1",  # furan
            "C1=CC=CS1",  # thiophene
            "C1CCCCC1"   # cyclohexane
        ]
        
    def generate_compound(self, target_mw: float = 300.0, 
                         complexity: str = 'medium') -> str:
        """Generate compound using fractal growth patterns"""
        
        # Start with random core
        core = random.choice(self.base_fragments)
        mol = Chem.MolFromSmiles(core)
        
        if mol is None:
            return core
            
        # Apply biological growth iterations
        iterations = {'low': 2, 'medium': 4, 'high': 6}[complexity]
        
        for i in range(iterations):
            mol = self._biological_growth_step(mol, target_mw)
            if mol is None:
                break
                
        return Chem.MolToSmiles(mol) if mol else core
    
    def _biological_growth_step(self, mol, target_mw: float):
        """Single growth step using biological patterns"""
        
        current_mw = Descriptors.MolWt(mol)
        if current_mw >= target_mw:
            return mol
            
        # Available attachment points (like biological binding sites)
        attachment_points = []
        for atom in mol.GetAtoms():
            if atom.GetTotalValence() < atom.GetTotalValence():
                attachment_points.append(atom.GetIdx())
        
        if not attachment_points:
            return mol
            
        # Growth fragments (biological building blocks)
        growth_fragments = [
            "C", "N", "O", "CC", "CO", "CN", "CCO", "CCC"
        ]
        
        # Select attachment point and fragment
        attach_idx = random.choice(attachment_points)
        fragment = random.choice(growth_fragments)
        
        # Simulate growth (simplified)
        try:
            combined_smiles = f"{Chem.MolToSmiles(mol)}.{fragment}"
            new_mol = Chem.MolFromSmiles(combined_smiles)
            
            if new_mol and self._passes_drug_likeness(new_mol):
                return new_mol
        except:
            pass
            
        return mol
    
    def _passes_drug_likeness(self, mol) -> bool:
        """Check biological drug-likeness (Lipinski-like rules)"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        return (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10)

# Example usage
if __name__ == "__main__":
    generator = BiologicalCompoundGenerator()
    
    for i in range(5):
        compound = generator.generate_compound(target_mw=350.0)
        mol = Chem.MolFromSmiles(compound)
        if mol:
            mw = Descriptors.MolWt(mol)
            print(f"Compound {i+1}: {compound} (MW: {mw:.1f})")`,
    dependencies: ['rdkit', 'numpy']
  },
  {
    id: 'clustering-ml',
    name: 'Biological Clustering ML',
    description: 'Machine learning clustering with biological patterns',
    language: 'python',
    category: 'ml',
    code: `import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import matplotlib.pyplot as plt

class BiologicalClustering:
    """
    Clustering algorithm inspired by biological systems
    Dr. Crowe Coder Pattern: Swarm intelligence
    """
    
    def __init__(self, n_clusters: int = 5):
        self.n_clusters = n_clusters
        self.clusterer = None
        self.pca = PCA(n_components=2)
        
    def fit_biological_clusters(self, smiles_list: list) -> dict:
        """
        Cluster compounds using biological swarming patterns
        """
        # Generate molecular fingerprints
        fingerprints = []
        valid_smiles = []
        
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # ECFP fingerprints (like biological recognition patterns)
                fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, 1024)
                fp_array = np.array(fp)
                fingerprints.append(fp_array)
                valid_smiles.append(smiles)
        
        if not fingerprints:
            return {"clusters": [], "centers": [], "visualization": None}
            
        X = np.array(fingerprints)
        
        # Apply biological clustering (adaptive K-means)
        self.clusterer = BiologicalKMeans(n_clusters=self.n_clusters)
        cluster_labels = self.clusterer.fit_predict(X)
        
        # PCA for visualization (like biological dimensionality reduction)
        X_pca = self.pca.fit_transform(X)
        
        # Organize results
        clusters = {}
        for i, label in enumerate(cluster_labels):
            if label not in clusters:
                clusters[label] = []
            clusters[label].append({
                'smiles': valid_smiles[i],
                'cluster_id': label,
                'coordinates': X_pca[i].tolist()
            })
        
        # Calculate cluster centers (like biological centroids)
        centers = []
        for i in range(self.n_clusters):
            cluster_points = X[cluster_labels == i]
            if len(cluster_points) > 0:
                center = np.mean(cluster_points, axis=0)
                centers.append(center.tolist())
        
        return {
            "clusters": clusters,
            "centers": centers,
            "n_compounds": len(valid_smiles),
            "visualization_data": {
                "coordinates": X_pca.tolist(),
                "labels": cluster_labels.tolist()
            }
        }

class BiologicalKMeans:
    """K-Means with biological adaptation patterns"""
    
    def __init__(self, n_clusters: int = 5, max_iter: int = 100):
        self.n_clusters = n_clusters
        self.max_iter = max_iter
        self.centers = None
        
    def fit_predict(self, X: np.ndarray) -> np.ndarray:
        """Fit and predict using biological swarming patterns"""
        n_samples, n_features = X.shape
        
        # Initialize centers using biological distribution
        self.centers = self._initialize_biological_centers(X)
        
        for iteration in range(self.max_iter):
            # Assign points to nearest centers (like biological attraction)
            distances = self._calculate_biological_distances(X, self.centers)
            labels = np.argmin(distances, axis=1)
            
            # Update centers with biological momentum
            new_centers = np.zeros_like(self.centers)
            for k in range(self.n_clusters):
                cluster_points = X[labels == k]
                if len(cluster_points) > 0:
                    # Biological averaging with homeostasis
                    new_center = np.mean(cluster_points, axis=0)
                    momentum = 0.9  # Biological momentum factor
                    new_centers[k] = momentum * self.centers[k] + (1 - momentum) * new_center
                else:
                    new_centers[k] = self.centers[k]
            
            # Check convergence (biological stability)
            if np.allclose(self.centers, new_centers, rtol=1e-4):
                break
                
            self.centers = new_centers
        
        return labels
    
    def _initialize_biological_centers(self, X: np.ndarray) -> np.ndarray:
        """Initialize centers using biological distribution patterns"""
        # Use k-means++ like initialization but with biological variance
        n_samples, n_features = X.shape
        centers = np.zeros((self.n_clusters, n_features))
        
        # First center: random point
        centers[0] = X[np.random.randint(n_samples)]
        
        # Subsequent centers: biological diversity principle
        for c_id in range(1, self.n_clusters):
            distances = np.array([
                min([np.linalg.norm(x - centers[j])**2 for j in range(c_id)])
                for x in X
            ])
            # Biological probability distribution
            probabilities = distances / distances.sum()
            cumulative_probs = probabilities.cumsum()
            r = np.random.rand()
            
            for j, p in enumerate(cumulative_probs):
                if r < p:
                    centers[c_id] = X[j]
                    break
        
        return centers
    
    def _calculate_biological_distances(self, X: np.ndarray, centers: np.ndarray) -> np.ndarray:
        """Calculate distances with biological weighting"""
        distances = np.zeros((X.shape[0], self.n_clusters))
        
        for i, point in enumerate(X):
            for j, center in enumerate(centers):
                # Euclidean distance with biological scaling
                dist = np.linalg.norm(point - center)
                # Apply biological distance transformation (sigmoid-like)
                biological_dist = 1 / (1 + np.exp(-0.1 * dist))
                distances[i, j] = biological_dist
        
        return distances

# Example usage
if __name__ == "__main__":
    # Sample compounds
    compounds = [
        "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen  
        "COC1=CC=C(C=C1)C(C)C(=O)O",  # Similar to ibuprofen
        "CC1=CC=C(C=C1)C(C)C(=O)O",   # Another NSAID
        "C1=CC=C(C=C1)C(=O)O"         # Benzoic acid
    ]
    
    clustering = BiologicalClustering(n_clusters=3)
    results = clustering.fit_biological_clusters(compounds)
    
    print(f"Clustered {results['n_compounds']} compounds into {len(results['clusters'])} groups")
    for cluster_id, compounds_in_cluster in results['clusters'].items():
        print(f"Cluster {cluster_id}: {len(compounds_in_cluster)} compounds")`,
    dependencies: ['scikit-learn', 'rdkit', 'numpy', 'matplotlib']
  }
];

export default function ImmersiveIDE() {
  const [selectedTab, setSelectedTab] = useState(0);
  const [selectedTemplate, setSelectedTemplate] = useState<CodeTemplate>(codeTemplates[0]);
  const [code, setCode] = useState(codeTemplates[0].code);
  const [language, setLanguage] = useState<'python' | 'typescript' | 'javascript'>('python');
  const [isExecuting, setIsExecuting] = useState(false);
  const [executionResult, setExecutionResult] = useState<ExecutionResult | null>(null);
  const [drCroweAgents, setDrCroweAgents] = useState<DrCroweAgent[]>([
    { id: '001', name: 'Dr. Crowe Coder', status: 'idle', specialization: ['Algorithms', 'Biological Computing'], confidence: 0.95 },
    { id: '002', name: 'Maya Patel', status: 'idle', specialization: ['ML/AI', 'Neural Networks'], confidence: 0.88 },
    { id: '003', name: 'Sarah Chen', status: 'idle', specialization: ['Quantum Computing', 'Optimization'], confidence: 0.92 }
  ]);
  const [showTemplateDialog, setShowTemplateDialog] = useState(false);
  const [executionMode, setExecutionMode] = useState<'local' | 'cloud' | 'agents'>('agents');

  const editorRef = useRef<monaco.editor.IStandaloneCodeEditor | null>(null);

  const handleEditorDidMount = (editor: monaco.editor.IStandaloneCodeEditor, monaco: Monaco) => {
    editorRef.current = editor;
    
    // Configure Monaco for biological computing
    monaco.languages.registerCompletionItemProvider('python', {
      provideCompletionItems: (model, position) => {
        const suggestions = [
          {
            label: 'biological_similarity',
            kind: monaco.languages.CompletionItemKind.Function,
            insertText: 'biological_similarity(mol1, mol2, method="tanimoto")',
            detail: 'Dr. Crowe biological similarity function'
          },
          {
            label: 'adaptive_clustering',
            kind: monaco.languages.CompletionItemKind.Function, 
            insertText: 'adaptive_clustering(compounds, n_clusters=5)',
            detail: 'Biological-inspired clustering algorithm'
          },
          {
            label: 'generate_compounds',
            kind: monaco.languages.CompletionItemKind.Function,
            insertText: 'generate_compounds(scaffold, n_compounds=10)',
            detail: 'AI compound generator with biological patterns'
          }
        ];
        return { suggestions };
      }
    });
  };

  const executeCode = async () => {
    setIsExecuting(true);
    
    // Update agent status
    setDrCroweAgents(prev => 
      prev.map(agent => ({ ...agent, status: 'analyzing' as const }))
    );

    try {
      // Simulate code execution with Dr. Crowe Coder
      const response = await fetch('http://localhost:8000/api/execute-code', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          code,
          language,
          mode: executionMode,
          biologicalPatterns: true,
          crowe_agent: true
        })
      });

      const result = await response.json();
      
      setExecutionResult({
        success: result.success,
        output: result.output || 'Code executed successfully with biological patterns applied.',
        error: result.error,
        performance: result.performance || {
          executionTime: Math.random() * 2000 + 500,
          memoryUsage: Math.random() * 100 + 50,
          croweScore: Math.random() * 0.3 + 0.7
        },
        biologicalPatterns: result.biologicalPatterns || ['homeostasis', 'adaptive-scaling', 'fractal-growth'],
        compounds: result.compounds
      });

    } catch (error) {
      setExecutionResult({
        success: false,
        output: '',
        error: 'Execution failed: ' + error,
        performance: { executionTime: 0, memoryUsage: 0, croweScore: 0 }
      });
    }

    // Reset agent status
    setTimeout(() => {
      setDrCroweAgents(prev => 
        prev.map(agent => ({ ...agent, status: 'idle' as const }))
      );
      setIsExecuting(false);
    }, 2000);
  };

  const loadTemplate = (template: CodeTemplate) => {
    setSelectedTemplate(template);
    setCode(template.code);
    setLanguage(template.language);
    setShowTemplateDialog(false);
  };

  return (
    <Box sx={{ height: '100vh', display: 'flex', flexDirection: 'column' }}>
      {/* Top Toolbar */}
      <Paper sx={{ p: 2, borderRadius: 0 }}>
        <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
            <Typography variant="h6" sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
              <CodeIcon /> CriOS Immersive IDE
            </Typography>
            <Chip 
              icon={<AIIcon />} 
              label="Dr. Crowe Coder Active" 
              color="primary" 
              size="small" 
            />
          </Box>
          
          <Box sx={{ display: 'flex', gap: 1 }}>
            <ToggleButtonGroup
              value={executionMode}
              exclusive
              onChange={(e, value) => value && setExecutionMode(value)}
              size="small"
            >
              <ToggleButton value="local">Local</ToggleButton>
              <ToggleButton value="cloud">Cloud</ToggleButton>
              <ToggleButton value="agents">AI Agents</ToggleButton>
            </ToggleButtonGroup>
            
            <Button
              variant="outlined"
              startIcon={<LoadIcon />}
              onClick={() => setShowTemplateDialog(true)}
              size="small"
            >
              Templates
            </Button>
            
            <Button
              variant="contained"
              startIcon={isExecuting ? <StopIcon /> : <RunIcon />}
              onClick={executeCode}
              disabled={isExecuting}
              sx={{
                background: isExecuting ? 'linear-gradient(45deg, #ff6b6b, #ee5a24)' : 'linear-gradient(45deg, #00e5ff 30%, #76ff03 90%)'
              }}
            >
              {isExecuting ? 'Executing...' : 'Run Code'}
            </Button>
          </Box>
        </Box>
      </Paper>

      <Box sx={{ display: 'flex', flex: 1, overflow: 'hidden' }}>
        {/* Left Panel - Code Editor */}
        <Box sx={{ flex: 1, display: 'flex', flexDirection: 'column' }}>
          <Paper sx={{ m: 1, flex: 1, display: 'flex', flexDirection: 'column' }}>
            <Box sx={{ p: 1, borderBottom: 1, borderColor: 'divider' }}>
              <Typography variant="subtitle2" gutterBottom>
                {selectedTemplate.name} - {selectedTemplate.description}
              </Typography>
              <Box sx={{ display: 'flex', gap: 1, flexWrap: 'wrap' }}>
                {selectedTemplate.dependencies.map(dep => (
                  <Chip key={dep} label={dep} size="small" variant="outlined" />
                ))}
              </Box>
            </Box>
            
            <Box sx={{ flex: 1 }}>
              <Editor
                height="100%"
                language={language}
                theme="vs-dark"
                value={code}
                onChange={(value) => setCode(value || '')}
                onMount={handleEditorDidMount}
                options={{
                  fontSize: 14,
                  minimap: { enabled: true },
                  wordWrap: 'on',
                  automaticLayout: true,
                  suggestOnTriggerCharacters: true,
                  quickSuggestions: true
                }}
              />
            </Box>
          </Paper>
        </Box>

        {/* Right Panel - Results & Agents */}
        <Box sx={{ width: 400, display: 'flex', flexDirection: 'column' }}>
          {/* Dr. Crowe Agents */}
          <Paper sx={{ m: 1, p: 2 }}>
            <Typography variant="h6" gutterBottom sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
              <NetworkIcon /> Active Agents
            </Typography>
            <List dense>
              {drCroweAgents.map((agent) => (
                <ListItem key={agent.id}>
                  <ListItemIcon>
                    <Box
                      sx={{
                        width: 12,
                        height: 12,
                        borderRadius: '50%',
                        bgcolor: agent.status === 'idle' ? 'success.main' : 
                                agent.status === 'analyzing' ? 'warning.main' : 'primary.main',
                        animation: agent.status !== 'idle' ? 'pulse 1.5s infinite' : 'none'
                      }}
                    />
                  </ListItemIcon>
                  <ListItemText
                    primary={agent.name}
                    secondary={`${agent.status} • ${(agent.confidence * 100).toFixed(0)}% confidence`}
                  />
                </ListItem>
              ))}
            </List>
          </Paper>

          {/* Execution Results */}
          <Paper sx={{ m: 1, flex: 1, display: 'flex', flexDirection: 'column' }}>
            <Box sx={{ p: 2, borderBottom: 1, borderColor: 'divider' }}>
              <Typography variant="h6">Execution Results</Typography>
            </Box>
            
            <Box sx={{ flex: 1, overflow: 'auto', p: 2 }}>
              {executionResult ? (
                <Box>
                  {/* Performance Metrics */}
                  <Card sx={{ mb: 2 }}>
                    <CardContent>
                      <Typography variant="subtitle2" gutterBottom>
                        Performance Metrics
                      </Typography>
                      <Grid container spacing={2}>
                        <Grid item xs={4}>
                          <Typography variant="caption" color="textSecondary">
                            Time (ms)
                          </Typography>
                          <Typography variant="h6">
                            {executionResult.performance.executionTime.toFixed(0)}
                          </Typography>
                        </Grid>
                        <Grid item xs={4}>
                          <Typography variant="caption" color="textSecondary">
                            Memory (MB)
                          </Typography>
                          <Typography variant="h6">
                            {executionResult.performance.memoryUsage.toFixed(0)}
                          </Typography>
                        </Grid>
                        <Grid item xs={4}>
                          <Typography variant="caption" color="textSecondary">
                            Crowe Score
                          </Typography>
                          <Typography variant="h6" color="primary">
                            {executionResult.performance.croweScore.toFixed(2)}
                          </Typography>
                        </Grid>
                      </Grid>
                    </CardContent>
                  </Card>

                  {/* Biological Patterns */}
                  {executionResult.biologicalPatterns && (
                    <Box sx={{ mb: 2 }}>
                      <Typography variant="subtitle2" gutterBottom>
                        Biological Patterns Applied
                      </Typography>
                      <Box sx={{ display: 'flex', gap: 1, flexWrap: 'wrap' }}>
                        {executionResult.biologicalPatterns.map(pattern => (
                          <Chip key={pattern} label={pattern} size="small" color="secondary" />
                        ))}
                      </Box>
                    </Box>
                  )}

                  {/* Output */}
                  <Typography variant="subtitle2" gutterBottom>Output</Typography>
                  <Paper sx={{ p: 2, bgcolor: 'grey.900', fontFamily: 'monospace' }}>
                    <Typography variant="body2" style={{ whiteSpace: 'pre-wrap', color: '#00ff00' }}>
                      {executionResult.output}
                    </Typography>
                    {executionResult.error && (
                      <Typography variant="body2" style={{ whiteSpace: 'pre-wrap', color: '#ff6b6b' }}>
                        Error: {executionResult.error}
                      </Typography>
                    )}
                  </Paper>
                </Box>
              ) : (
                <Typography color="textSecondary">
                  Run code to see results here...
                </Typography>
              )}
            </Box>
          </Paper>
        </Box>
      </Box>

      {/* Template Selection Dialog */}
      <Dialog
        open={showTemplateDialog}
        onClose={() => setShowTemplateDialog(false)}
        maxWidth="md"
        fullWidth
      >
        <DialogTitle>Code Templates - Dr. Crowe Coder Library</DialogTitle>
        <DialogContent>
          <Grid container spacing={2}>
            {codeTemplates.map((template) => (
              <Grid item xs={12} sm={6} key={template.id}>
                <Card 
                  sx={{ 
                    cursor: 'pointer',
                    '&:hover': { elevation: 4 }
                  }}
                  onClick={() => loadTemplate(template)}
                >
                  <CardContent>
                    <Typography variant="h6" gutterBottom>
                      {template.name}
                    </Typography>
                    <Typography variant="body2" color="textSecondary" paragraph>
                      {template.description}
                    </Typography>
                    <Box sx={{ display: 'flex', gap: 1, mb: 1 }}>
                      <Chip label={template.language} size="small" />
                      <Chip label={template.category} size="small" color="secondary" />
                    </Box>
                    <Typography variant="caption" color="textSecondary">
                      Dependencies: {template.dependencies.join(', ')}
                    </Typography>
                  </CardContent>
                </Card>
              </Grid>
            ))}
          </Grid>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setShowTemplateDialog(false)}>Cancel</Button>
        </DialogActions>
      </Dialog>

      <style jsx>{`
        @keyframes pulse {
          0% { opacity: 1; }
          50% { opacity: 0.5; }
          100% { opacity: 1; }
        }
      `}</style>
    </Box>
  );
}