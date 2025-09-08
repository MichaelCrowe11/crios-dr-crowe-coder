# Developer Playbook

## Monorepo Conventions
- TypeScript: web (Vite+React), desktop (Electron), API (Express+SSE).  
- Python: `crios` package (cheminformatics core, pipelines, CLI).  
- Testing: Vitest + Playwright (TS); Pytest (coverage ≥80%) for Python.  
- Lint/Type: ESLint/Prettier (TS), ruff + mypy (`--strict-optional`) for Python.  
- Security: rate limiting, schema validation, redacted structured logs, sandboxed execution.  
- Docs: README (quickstart), `/docs` (architecture, API, agents, examples), SECURITY.md, CONTRIBUTING.md.  

## CLI Layer
- Node: `crowe assist|complete|run|bench`.  
- Python: `crios validate|search|cluster|process|score|discover|benchmark|web`.  
- Consistent `--json` output; stable exit codes.  

## Scientific Pipeline Pattern
Validate → Filter → Describe → Search/Cluster → Rank → (Dock/MD optional) → Report.  
Filters: Lipinski, PAINS, reactive groups, SA score, novelty bounds.  
Ranking rubric: potency proxy, ADMET flags, QED, novelty, synthesizability, route risk.  

## Benchmarks
MoleculeNet (BACE, BBBP, Tox21, ClinTox), FS-Mol few-shot, docking sanity sets.  
Use fixed seeds, scaffold splits, JSON metrics logs.  

## Quality Gates (PR)
- Tests green; coverage delta ≥ 0.  
- Lints pass.  
- Docs & changelog updated for CLI/API changes.  
- Security checklist completed (inputs validated; no secrets leakage).  
- ADR for material design shifts.  

## Performance & Repro
- Deterministic seeds.  
- Log configuration + dataset hash.  
- Cache heavy descriptors/fingerprints.  

## Observability
- Structured JSON logs (event, component, duration_ms, status).  
- Tracing for long-running screens (span: generate / dock / score).  

## Release
- Semantic versioning.  
- Pre-release channel for experimental models.  

## Decommission / Deprecation
- Mark deprecated APIs with sunset date.  
- Migration guides in `/docs/migrations`.  
