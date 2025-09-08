# Dr. Crowe Coder Orchestration Prompts

This directory contains the master system prompt, developer playbook, reusable user prompt templates, and mini examples implementing the Crowe Logic methodology.

## Files

- `master_system_prompt.md` – System identity, divisions, Crowe Logic method, output contracts, safety & compliance, style.
- `developer_playbook.md` – Monorepo, pipelines, benchmarks, quality gates, PR standards.
- `user_prompt_templates.md` – Reusable user prompts (Blueprint program, QSAR baseline, Generative design, PBPK/PD sketch, IND scaffolds, Quantum add-on).
- `mini_examples.md` – Few-shot paired Python/R snippets.
 - `response_skeleton.md` – Canonical output section ordering and labels.
 - `acceptance_criteria.md` – Checklist for verifying deliverables.

A CLI helper can emit these prompts for integration with external orchestrators.

## Crowe Logic Sections
1. OBSERVE  
2. DECOMPOSE  
3. CONNECT  
4. SYNTHESIZE  
5. VALIDATE  

All assistant outputs should include a Next-Actions Checklist, assumptions, limits, selection rubric (when comparing options), dual-language code (unless user restricts), and the research-planning disclaimer.

## Disclaimer
**Research planning only. No clinical advice or efficacy claims; verification and regulatory review required before any real-world use.**
