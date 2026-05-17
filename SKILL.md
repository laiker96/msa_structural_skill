---
name: structural-evo-analysis
description: Run and adapt a reproducible protein vulnerability-analysis pipeline using a diverse MMseqs2 homolog subset, OGT-aware MSA metadata, conservation, CamSol-style solubility scores, Aggrescan3D aggregation scores, and residue-level vulnerability ranking. Use when Codex needs to identify vulnerable enzyme positions compared with homologs and available query/AFDB structures; optional OGT clade enrichment is supported but not required.
---

# Protein Vulnerability Analysis

Use `scripts/structural_evo_analysis/` as the maintained pipeline.

## Quick Orientation

- Read `scripts/structural_evo_analysis/README.md` for step-level commands.
- Treat the directory containing this `SKILL.md` as the skill root, even when
  it has been copied under `~/.codex/skills/`.
- Runtime data defaults outside the skill root:
  `SEA_WORK_DIR=~/structural_evo_analysis`, with results in
  `~/structural_evo_analysis/results` and structures in
  `~/structural_evo_analysis/structures`.
- Use `sequences/photoHymenobact.fa` as the bundled example query.
- Use `data/ogt_taxid_summary.tsv` as the only bundled OGT metadata file.
- Use `<skill-root>/envs/structural_evo/bin/python` after running
  `<skill-root>/setup_envs.sh`.

## Workflow

1. Confirm the user provided exactly one query protein FASTA and one query PDB
   structure. The PDB can be crystal, AlphaFold, AF3, or another query model,
   but it is required for vulnerability ranking.
2. Ask whether to use the default work directory
   `~/structural_evo_analysis` or a user-selected directory. Set
   `SEA_WORK_DIR` for the default directory choice, or set `SEA_OUT_DIR` and
   `SEA_STRUCTURE_DIR` when results and structures should live separately.
3. Ask for the sequence database location. Accept either exact
   `SEA_DB_FASTA`/`SEA_DB_MMSEQS` paths or `UNIREF_DIR` plus `SEA_DB`
   (`50`, `90`, `100`, `uniref50`, `uniref90`, `uniref100`). If the database is
   absent, ask whether the user wants UniRef downloaded before running
   `prepare_uniref_database.py --download`.
4. Check whether `<skill-root>/envs/structural_evo/bin/python` exists. If not,
   run setup in `tmux` from the skill root and log it, for example
   `tmux new -s structural-evo-setup 'cd <skill-root> && bash setup_envs.sh 2>&1 | tee logs/setup_envs.tmux.log'`.
5. If Aggrescan3D scoring is required, check for
   `<skill-root>/envs/aggrescan3d/bin/aggrescan`; if missing, run
   `bash setup_envs.sh --skip-main --with-aggrescan3d` in `tmux`.
6. Run the numbered scripts in order, or use `run_pipeline.sh`. Step 01 uses a
   diverse identity-stratified subset by default; set `SEA_MAX_REPSET_SEQS=0`
   only when the user explicitly wants all filtered hits.
7. Set `SEA_QUERY_PDB=/path/to/query.pdb` or place the query PDB at
   `~/structural_evo_analysis/structures/query.pdb` when using the default work
   directory.
8. OGT-aware mode is optional and off by default. Enable it with
   `SEA_OGT_AWARE=1`; this runs OGT clade calling and makes logos/viewers
   clade-aware. For manual OGT clade calls, pass `repset_metadata.tsv` to step 04 with
   `--trait-column ogt --low-threshold 20 --high-threshold 45`.
9. Run full searches, IQ-TREE, downloads, scoring, and vulnerability analysis
   inside `tmux` with a log.
10. Report exact commands, outputs, and any steps not run.

## Important Constraints

- Do not silently change search thresholds, tree settings, clade thresholds,
  filters, or scoring parameters.
- The agent may adjust MMseqs2/search/subset parameters when needed to obtain a
  usable homolog set, but must keep changes bounded, document exact values, and
  explain why they were changed.
- The primary output is `vulnerability/top_vulnerable_positions.tsv`, which
  combines conservation, CamSol-style solubility, and Aggrescan3D aggregation
  scores. Treat it as prioritization, not experimental validation.
- OGT enrichment/clade analysis is optional context, not the default endpoint.
- Viewer outputs are self-contained HTML files under `viewers/`; sequence-logo
  images are under `logos/`.
- Keep generated outputs under `~/structural_evo_analysis/results` and
  downloaded/generated structures under `~/structural_evo_analysis/structures`
  unless the user chooses another `SEA_WORK_DIR`, `SEA_OUT_DIR`, or
  `SEA_STRUCTURE_DIR`.
- Do not download UniRef or large AFDB sets without making the storage location
  and expected command explicit to the user. Use `tmux` and logs for downloads.
- Treat `scripts/msa_OGT/` as legacy/reference code. Its scorer entrypoints are
  compatibility wrappers around `scripts/structural_evo_analysis/`.
- Use `--skip-aggrescan3d` for structure scoring if the optional
  Aggrescan3D environment has not been installed.
