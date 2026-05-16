# Structural Evolution Skill Guide

## Scope

This repository is a reusable skill/pipeline scaffold for single-protein
structural evolutionary analysis. The main maintained workflow is:

```text
scripts/structural_evo_analysis/
```

Use it to search homologs, build an MSA/tree, annotate clades from metadata,
compute conserved positions, download AlphaFold DB structures, and optionally
score structures for solubility/aggregability.

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
- Search database: a local protein FASTA such as UniRef, configured with
  `SEA_DB_FASTA`/`SEA_DB_MMSEQS` or `UNIREF_DIR`/`SEA_DB`.
- Metadata: `data/growth_temp_dataset_OGTFinder.tsv` is the default
  taxonomy-keyed growth-temperature source. Step 01 propagates `TaxID=` from
  UniRef-style headers and joins this table into `repset_metadata.tsv`.
- Structures: AlphaFold DB PDB files downloaded by step 06 under
  `structures/structural_evo_analysis/afdb/`; an optional query PDB can be
  placed at `structures/structural_evo_analysis/query.pdb`.

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
  'SEA_OUT_DIR=results/photoHymenobact_example SEA_LOW_THRESHOLD=20 SEA_HIGH_THRESHOLD=45 bash scripts/structural_evo_analysis/run_pipeline.sh sequences/photoHymenobact.fa results/photoHymenobact_example/repset_metadata.tsv ogt 2>&1 | tee logs/photoHymenobact_example.log'
```

Record the command, working directory, inputs, outputs/logs, relevant
environment variables, and how to monitor or stop the task.

## Pipeline Rules

- Preserve the numbered pipeline order unless a change is explicitly requested.
- Do not silently change defaults, search thresholds, clade thresholds, filters,
  or scoring parameters.
- Keep new outputs under `results/structural_evo_analysis/` or a
  user-selected `SEA_OUT_DIR`.
- Keep downloaded/generated structures under
  `structures/structural_evo_analysis/` or a user-selected
  `SEA_STRUCTURE_DIR`.
- Treat `scripts/msa_OGT/` as legacy/reference code except for the scorer
  modules currently reused by step 07.
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
