# ANKROS Enzyme Project Guide

## Scope

This repository centers on ANKros-CPD, a class I CPD photolyase from
*Hymenobacter* sp. UV11. Align changes with the existing pipeline layout.

## Biology Reference

- `ANKros` = *Photohymenobacter* photolyase.
- Protein: ANKros-CPD / ANW48624.1.
- Domains:
  - `1-130`: antenna domain.
  - `131-205`: interdomain linker.
  - `206-437`: catalytic/helical domain.
- Cofactors: FAD and FMN.

## Info Directory

Use `info/INDEX.md` before opening PDFs or raw Design-Expert files for extra
biological context.

## Environments

`setup_envs.sh` defines how environments are created. This file only documents
how to use them.

### Main environment

Use `./envs/ankros/bin/python` for MSA, structural work, and most experimental
scripts.

```bash
./envs/ankros/bin/python scripts/msa_OGT/01_mmseqs_search.py
./envs/ankros/bin/python scripts/structural/orchestrate_holo_model.py
```

### Specialized environments

- Amber CUDA: `./envs/amber_cuda/bin/pmemd.cuda`
- Aggrescan3D: `./envs/aggrescan3d/bin/aggrescan`
- ThermoMPNN: `./envs/thermompnn/bin/python`

Examples:

```bash
PMEMD=./envs/amber_cuda/bin/pmemd.cuda bash scripts/experiments/fad_fmn_md_modeling/04_run_amber_cuda_juanma_ramp.sh
./envs/aggrescan3d/bin/aggrescan --help
./envs/thermompnn/bin/python -c "import torch, pandas, omegaconf, Bio"
```

ThermoMPNN is cloned/configured under `external/ThermoMPNN/`.

`setup_envs.sh` is idempotent and logs to `logs/setup_envs.log`.

## Experimental Pipelines

Experimental work belongs under:

```text
scripts/experiments/<topic>/
```

Outputs belong under:

```text
results/experiments/<topic>/
```

Recommended layout:

```text
scripts/experiments/<topic>/run_pipeline.sh
scripts/experiments/<topic>/README.md
results/experiments/<topic>/{figures,tables,logs,intermediates}/
```

Use existing shared outputs instead of duplicating them:

- MSA outputs: `results/msa_OGT/`
- Structural outputs: `results/structural/`

Use `./envs/ankros/bin/python` for new experimental pipelines unless a
specialized environment is explicitly required.

## Project-Specific Cautions

- Preserve residue numbering.
- Preserve existing pipeline numbering and order.
- Do not silently change biological thresholds, cutoffs, or residue mappings.
- When adding an environment, place its config under `config/`, wire it into
  `setup_envs.sh`, and document the exact invocation command.
