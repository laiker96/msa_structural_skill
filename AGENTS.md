# Protein Vulnerability Skill Guide

## Scope

This repository is a reusable skill/pipeline scaffold for single-protein
vulnerability analysis. The main maintained workflow is:

```text
scripts/structural_evo_analysis/
```

Use it to search homologs, build a diverse OGT-aware MSA, compute conserved
positions, score available query/AFDB structures for solubility/aggregability,
and rank vulnerable query-enzyme positions. OGT enrichment/clade analysis is
optional context, not the default endpoint.

The ANKros/photolyase sequence in `sequences/photoHymenobact.fa` is the bundled
test query. Do not hard-code ANKros-specific residue ranges, domain boundaries,
thresholds, or biological interpretations into the generalized pipeline.

## Operating Priorities

1. Correctness
2. Reproducibility
3. Safe handling of long-running tasks
4. Clarity
5. Minimal changes
6. Honest uncertainty

## Main Inputs

- Query FASTA: `sequences/photoHymenobact.fa` for the example, or a user
  supplied single-protein FASTA.
- Query PDB: required for vulnerability ranking. The user should provide a
  crystal, AF3, AlphaFold, or other query structure via `SEA_QUERY_PDB` or
  `~/structural_evo_analysis/structures/query.pdb`.
- Search database: a local protein FASTA such as UniRef, configured with
  `SEA_DB_FASTA`/`SEA_DB_MMSEQS` or `UNIREF_DIR`/`SEA_DB`. UniRef shorthand
  values `50`, `90`, `100`, `uniref50`, `uniref90`, and `uniref100` are
  accepted. Ask for the database location before running; if it is absent, ask
  whether to download UniRef and use `prepare_uniref_database.py --download`
  only after approval.
- Metadata: `data/ogt_taxid_summary.tsv` is the only bundled data file and the
  default taxonomy-keyed OGT source. Step 01 propagates `TaxID=` from
  UniRef-style headers and joins OGT/regime values into `repset_metadata.tsv`.
- Structures: AlphaFold DB PDB files downloaded by step 06 under
  `~/structural_evo_analysis/structures/afdb/` by default; the required query
  PDB can be provided with `SEA_QUERY_PDB` or placed at
  `~/structural_evo_analysis/structures/query.pdb`.

## Runtime Storage

The default work directory is:

```text
~/structural_evo_analysis/
```

This keeps generated run data outside the skill root, which may live under
`~/.codex/skills/`. Ask whether the user wants this default or a different
directory before launching a run. Set `SEA_WORK_DIR` for the base directory, or
set `SEA_OUT_DIR` and `SEA_STRUCTURE_DIR` separately.

## Environment

`setup_envs.sh` is the source of truth for environment creation.

```bash
bash setup_envs.sh
```

This creates the main environment at:

```text
envs/structural_evo/
```

Use its Python explicitly:

```bash
./envs/structural_evo/bin/python scripts/structural_evo_analysis/01_mmseqs_search.py --help
```

Aggrescan3D is optional and separate:

```bash
bash setup_envs.sh --with-aggrescan3d
```

Do not add Amber, MD, PyRosetta, ThermoMPNN, OpenMM, CUDA, notebooks, or other
large stacks unless a maintained structural-evolution step actually requires
them.

## Long-Running Tasks

Run full searches, IQ-TREE jobs, AFDB downloads, and structure scoring in
`tmux`, with logs:

```bash
mkdir -p logs
tmux new -s structural-evo \
  'SEA_QUERY_PDB=/path/to/query.pdb SEA_WORK_DIR=~/structural_evo_analysis bash scripts/structural_evo_analysis/run_pipeline.sh sequences/photoHymenobact.fa 2>&1 | tee logs/photoHymenobact_example.log'
```

Record the command, working directory, inputs, outputs/logs, relevant
environment variables, and how to monitor or stop the task.

## Pipeline Rules

- Preserve the numbered pipeline order unless a change is explicitly requested.
- Treat `vulnerability/top_vulnerable_positions.tsv` as the primary maintained
  endpoint. OGT/regime clade enrichment is optional context and enabled with
  `SEA_OGT_AWARE=1`.
- Do not silently change defaults, search thresholds, clade thresholds, filters,
  or scoring parameters.
- Keep new outputs under `~/structural_evo_analysis/results/` or a
  user-selected `SEA_WORK_DIR`/`SEA_OUT_DIR`.
- Keep downloaded/generated structures under
  `~/structural_evo_analysis/structures/` or a user-selected
  `SEA_WORK_DIR`/`SEA_STRUCTURE_DIR`.
- Treat `scripts/msa_OGT/` as legacy/reference code. Its scorer entrypoints are
  compatibility wrappers around `scripts/structural_evo_analysis/`.
- Update documentation when behavior, inputs, outputs, parameters,
  dependencies, or assumptions change.

## Validation

After edits, run the smallest relevant checks available:

```bash
./envs/structural_evo/bin/python -m py_compile scripts/structural_evo_analysis/*.py
bash setup_envs.sh --help
```

For full-pipeline validation, use the bundled query sequence and a local
UniRef-compatible database. Do not claim the pipeline completed unless the logs
and expected output files were inspected.
