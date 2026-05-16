---
name: structural-evo-analysis
description: Run and adapt a reproducible single-protein structural evolution pipeline. Use when Codex needs to analyze a query protein with homolog search, MSA, phylogeny, metadata-driven clade annotation, conservation, AlphaFold DB structure download, or optional CamSol-style/Aggrescan3D structure scoring; also use when building a skill or workflow that reuses parts of scripts/structural_evo_analysis.
---

# Structural Evolution Analysis

Use `scripts/structural_evo_analysis/` as the maintained pipeline.

## Quick Orientation

- Read `scripts/structural_evo_analysis/README.md` for step-level commands.
- Use `sequences/photoHymenobact.fa` as the bundled example query.
- Use `data/growth_temp_dataset_OGTFinder.tsv` as taxonomy-keyed metadata for
  OGT/regime clade annotation.
- Use `./envs/structural_evo/bin/python` after running `bash setup_envs.sh`.

## Workflow

1. Confirm there is exactly one query protein sequence.
2. Confirm the local protein database path through `SEA_DB_FASTA` and
   `SEA_DB_MMSEQS`, or through `UNIREF_DIR` and `SEA_DB`.
3. Run the numbered scripts in order, or use `run_pipeline.sh`.
4. For OGT clade calls, pass `repset_metadata.tsv` to step 04 with
   `--trait-column ogt --low-threshold 20 --high-threshold 45`.
5. Run full searches, IQ-TREE, downloads, and scoring inside `tmux` with a log.
6. Report exact commands, outputs, and any steps not run.

## Important Constraints

- Do not silently change search thresholds, tree settings, clade thresholds,
  filters, or scoring parameters.
- Keep generated outputs under `results/` and downloaded structures under
  `structures/`.
- Treat `scripts/msa_OGT/` as legacy/reference code except for the scorer
  modules reused by step 07.
- Use `--skip-aggrescan3d` for structure scoring if the optional
  Aggrescan3D environment has not been installed.
