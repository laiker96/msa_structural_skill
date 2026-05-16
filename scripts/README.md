# `scripts/`

Two coordinated pipelines plus shared utilities. Every script here is
part of a concrete pipeline; experimental / ad-hoc scripts live in
`scripts/experiments/` until they mature into pipeline steps.

```
scripts/
├── msa_OGT/             01–17  UniRef taxid → OGTFinder → profile-HMM MSA
├── structural/          FAD/FMN holo-model entrypoint and structural helpers
├── solubility/                 Legacy CamSol-style ANKros/AFDB scoring
├── aggregability/              Legacy Aggrescan3D ANKros/AFDB scoring
└── experiments/                MD and ad-hoc analyses
```

## Existing dependency chain

```
   ┌─────────── scripts/msa_OGT/ ──────────┐
   │                                      │
   │  UniRef90 → MMseqs2 search →         │
   │  class I filter → master set →      │
   │  HMM align → QC → linker refine →   │
   │  IQ-TREE → regime clades/logos      │
   │                                      │
   │  Output: results/msa_OGT/*           │
   └───────────────────┬──────────────────┘
                       │ classI_confirmed.fa + tree
                       ▼
   ┌─────── scripts/structural/ ──────────┐
   │                                      │
   │  Build the DiffDock-free FAD/FMN     │
   │  holo model from 1TEZ and 2J09       │
   │  donors; write measured FAD/FMN,     │
   │  CPD, and DNA residue annotations.   │
   │                                      │
   │  Output: results/structural/         │
   └──────────────────────────────────────┘
```

Each pipeline subdirectory has its own `README.md` with step-by-step details.
`msa_OGT/14_validate_structural_matches.py` screens AFDB models against the
ANKros fold with Foldseek before structure scoring. `msa_OGT/15_compute_solubility_aggregability.py` runs the pipeline-owned
CamSol and Aggrescan3D scorer utilities and merges their structure scores for
the MSA pipeline. The step 16 HTML reads that merged table when present, and
adds ANKros plus representative-structure score heatmaps for each displayed
clade.

## Common requirements

All scripts assume the `ankros` conda env is active:

```bash
conda activate ./envs/ankros
```

Aggrescan3D and Amber CUDA are optional separate envs — see `setup_envs.sh`.
The structural holo builder stages its required crystals under
`results/structural/inputs/`. AFDB homolog models for MSA representatives are
fetched with `msa_OGT/13_download_afdb.py`.
Organism OGT data is under `data/` (BacDive cache, TEMPURA CSV).
