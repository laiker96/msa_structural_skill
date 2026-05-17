---
name: structural-evo-analysis
description: Use modular protein MSA and structure-analysis scripts for diversity-first homolog search, optional OGT-aware clade annotation, conservation scoring, structure scoring, and vulnerable-residue ranking. Use the full workflow only when the user needs end-to-end vulnerability analysis; otherwise run the smallest relevant module recipe.
---

# Protein MSA And Structure Analysis

Use `scripts/structural_evo_analysis/` as the maintained module set. The full
pipeline is one recipe, not the only way to use this skill.

## Quick Orientation

- Read `REFERENCE.md` for module recipes, inputs, outputs, and environment
  variables.
- Treat the directory containing this `SKILL.md` as the skill root, even when it
  has been copied under `~/.codex/skills/`.
- Runtime data defaults outside the skill root:
  `SEA_WORK_DIR=~/structural_evo_analysis`, with results in
  `~/structural_evo_analysis/results` and structures in
  `~/structural_evo_analysis/structures`.
- Use `test/photoHymenobact.fa` as the bundled example query.
- Keep committed smoke-test inputs under `test/`. Keep generated structures in
  the runtime `SEA_STRUCTURE_DIR`.
- Use `data/ogt_taxid_summary.tsv` only when OGT-aware mode is requested.
- Install the environment from the skill root when missing:
  `bash setup_envs.sh`. Use `<skill-root>/envs/structural_evo/bin/python`
  after setup.

## Choose The Smallest Recipe

- Diversity MSA only: run step 01 and step 02, or use
  `SEA_PIPELINE_MODE=msa`. No query PDB is required.
- Conservation only: run step 05 from an existing alignment, or use
  `SEA_PIPELINE_MODE=conservation` to run search, alignment, optional clade
  annotation, and conservation. No query PDB is required.
- OGT context: use `SEA_OGT_AWARE=1` or step 01 `--join-ogt`, then run tree and
  clade annotation. Non-OGT runs must not join OGT during search.
- Structure scoring only: run step 07 with query/AFDB PDBs. Use
  `--skip-aggrescan3d` if the optional Aggrescan3D environment is unavailable.
- Vulnerability ranking: combine conservation and structure scores with steps
  08 and 09. A real query PDB is required.
- Full workflow: use `run_pipeline.sh` with default `SEA_PIPELINE_MODE=full`.
  This requires a query PDB and may run long searches, tree inference, AFDB
  downloads, and structure scoring.

## Workflow Rules

1. Confirm the user’s actual goal before launching a full run. Do not require a
   PDB for MSA-only, OGT-only, or conservation-only tasks.
2. Ask whether to use the default work directory
   `~/structural_evo_analysis` or a user-selected directory. Set
   `SEA_WORK_DIR`, or set `SEA_OUT_DIR` and `SEA_STRUCTURE_DIR` separately.
3. Ask for the sequence database location before step 01. Accept exact
   `SEA_DB_FASTA`/`SEA_DB_MMSEQS` paths or `UNIREF_DIR` plus `SEA_DB`
   (`50`, `90`, `100`, `uniref50`, `uniref90`, `uniref100`). If absent, ask
   before downloading UniRef with `prepare_uniref_database.py --download`.
4. Check whether `<skill-root>/envs/structural_evo/bin/python` exists. If not,
   run `bash setup_envs.sh` in `tmux` from the skill root and log it. For
   structure scoring with Aggrescan3D, install the optional environment with
   `bash setup_envs.sh --skip-main --with-aggrescan3d` when
   `<skill-root>/envs/aggrescan3d/bin/aggrescan` is missing.
5. Run long searches, IQ-TREE jobs, AFDB downloads, and structure scoring in
   `tmux` with a log.
6. Report exact commands, key outputs, and any steps not run.

## Important Constraints

- Do not silently change search thresholds, tree settings, clade thresholds,
  filters, scoring parameters, or representative-set size.
- Step 01 uses a diverse identity-stratified subset by default; set
  `SEA_MAX_REPSET_SEQS=0` only when the user explicitly wants all filtered hits.
- Non-OGT runs should not join OGT metadata during search. Enable OGT mode with
  `SEA_OGT_AWARE=1` or step 01 `--join-ogt`.
- The primary full-workflow endpoint is
  `vulnerability/top_vulnerable_positions.tsv`. Treat it as prioritization, not
  experimental validation.
- Keep generated outputs under `~/structural_evo_analysis/results` and
  downloaded/generated structures under `~/structural_evo_analysis/structures`
  unless the user chooses another location.
- Do not download UniRef or large AFDB sets without making the storage location
  and expected command explicit to the user.
- Treat `scripts/msa_OGT/` as legacy/reference code. Its scorer entrypoints are
  compatibility wrappers around `scripts/structural_evo_analysis/`.
