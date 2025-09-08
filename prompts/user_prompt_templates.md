# Reusable User Prompt Templates

## A) Blueprint a Full Discovery Program
Goal: <program goal>
Inputs: <targets, SMILES, assay table, constraints>
Constraints: <MW<500, cLogP 1–3, patent exclusions>
Deliverable: Crowe Logic (OBSERVE→VALIDATE), options matrix, scoring rubric, next-actions. Include Python+R code to: validate SMILES, compute descriptors, baseline model, rank shortlist.

## B) QSAR/QSPR Baseline (Scaffold Split)
Dataset: <name>
Features: MorganFP(2,2048), MACCS, Mordred.  
Models: RF, XGBoost, MLP.  
Metrics: AUC/PRC or RMSE with CIs (bootstrap).  
Deliverable: Comparison table, Python+R code, feature importance snapshot.

## C) Generative Library Design + Filters
Objective: Generate 500 candidates; 3 scaffolds + 12 fragments.  
Constraints: Ro5, SA≤6, no PAINS.  
Novelty reference: <reference.smi>.  
Deliverable: Ranked list, novelty stats, minimal docking rationale.

## D) PBPK/PD Modeling Sketch
Compound: <compound>  
Indication: <disease>  
Deliverable: Coupled ODE system (Python+R), identifiability flags, simulation grid.

## E) IND/IMPD Module Scaffolds
Program: <program>  
Deliverable: Module 2 & 3 outlines, risk register (CMC/nonclinical/clinical), mitigations, document checklist.

## F) Quantum Add-On (Optional)
If value-add: run laptop-safe conformer scoring (quantum simulated), derive q-feature, integrate into ranking.  
Deliverable: Feature definition, integration rationale, resource caveats.

## Kickoff Prompt Template
Project: <program or target>
Scope: IDE + API + CLI + pipelines (validate/filter/score/generate/benchmark)
Inputs: <SMILES/FASTA/CSV>
Constraints: <Ro5, SA≤6, timeline N weeks>
Deliverable: Architecture blueprint, prioritized backlog, pipeline specs & baseline code (Python+R), CLI/API contracts, benchmark plan, next-actions checklist.  
House rules: Dual-language code, mapping table, scoring rubric, assumptions/limits, safety disclaimer, optional quantum example.

## Disclaimer
**Research planning only. No clinical advice or efficacy claims; verification and regulatory review required before any real-world use.**
