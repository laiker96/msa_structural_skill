# ANKros-CPD

Structural and biophysical analysis of **ANKros-CPD**, a class I CPD photolyase
(EC 4.1.99.3) from the Antarctic bacterium *Hymenobacter* sp. UV11 (Bacteroidota).
GenBank: KX118295.1 · UniProt/NCBI: ANW48624.1 · 437 aa · Tm ≈ 49 °C.

## Repository layout

```
scripts/                  # Active pipelines and experiments
├── README.md                 top-level pipeline map
├── msa_OGT/             01-16 UniRef -> class-I filter -> HMM/MSA -> tree/views
├── structural/          FAD/FMN holo builder and structural utilities
└── experiments/         MD, cell-membrane, GNM/ANM, and oligomer analyses

config/                   # Conda specs and installers used by setup_envs.sh
├── environment.yml           main ankros environment
├── amber_cuda_environment.yml
├── aggrescan3d_environment.yml
├── thermompnn_environment.yml
└── install_amber_cuda.sh

data/                     # Small committed inputs and parsed metadata
├── 200617_TEMPURA.csv
├── growth_temp_dataset_OGTFinder.tsv
├── bacdive/
└── DoE/

sequences/                # Reference/query FASTAs and AF3 submission FASTAs
├── README.md
└── af3/

structures/               # External/re-downloadable structure pools
├── README.md
├── af3/                      AF3 server predictions
└── afdb/                     AlphaFold DB models for MSA representatives

results/                  # Pipeline outputs retained or regenerated locally
├── msa_OGT/                  MSA handoff files and tree outputs
└── structural/               holo-builder inputs, QC, and FAD/FMN model

info/                     # Reference papers, enzyme notes, and DoE parser
archive/                  # Retired legacy workflows kept for auditability
best_results.md           # Snapshot of notable current results
requirements-lock.txt     # Post-install package-version audit

envs/                     # Local conda environments (generated, gitignored)
external/                 # Local third-party checkouts/builds (generated, gitignored)
logs/                     # Setup/build logs (generated, gitignored)
```

Each of `scripts/`, `scripts/msa_OGT/`, `scripts/structural/`, `structures/`,
and `sequences/` carries its own `README.md` with step-by-step inputs and
outputs. Start there for sub-pipeline detail.

## Environment setup

Use the repository setup script. It is idempotent, logs to
`logs/setup_envs.log`, and creates the main environment at `envs/ankros` from
`config/environment.yml`. The main solve is long-running, so launch it in
`tmux`:

```bash
tmux new -s ankros-setup 'bash setup_envs.sh'
```

Run commands with the environment Python directly:

```bash
./envs/ankros/bin/python scripts/msa_OGT/01_mmseqs_search.py
bash scripts/structural/run_pipeline.sh
```

PyTorch is pinned in the pip block of `config/environment.yml` as the `cu121`
wheel. The conda CUDA stack is retained for AmberTools, OpenFF, and OpenMM.
Exact pinned package versions from the installed environment are also captured
in `requirements-lock.txt` for audit.

Optional environments are separate and only created when requested:

```bash
# Amber/pmemd CUDA; licensed package stays outside git.
tmux new -s amber-cuda \
  'AMBER_CUDA_PACKAGE=~/Downloads/pmemd26.tar.bz2 bash setup_envs.sh --with-amber-cuda --amber-cuda-jobs 2'

# Aggrescan3D scorer.
tmux new -s aggrescan3d 'bash setup_envs.sh --skip-ankros --with-aggrescan3d'

# ThermoMPNN checkout + environment under external/ThermoMPNN and envs/thermompnn.
tmux new -s thermompnn 'bash setup_envs.sh --skip-ankros --with-thermompnn'
```

### Optional Amber/pmemd CUDA

The repo does not track licensed Amber source packages. To build a local
Amber/pmemd CUDA archive, keep the package outside git and pass its path
explicitly:

```bash
tmux new -s amber-cuda 'AMBER_CUDA_PACKAGE=~/Downloads/pmemd26.tar.bz2 bash config/install_amber_cuda.sh'
```

`setup_envs.sh --with-amber-cuda` creates a separate `envs/amber_cuda` build
environment from `config/amber_cuda_environment.yml` because
Amber/pmemd rejects CUDA 12.9 and newer. The installer extracts under ignored
`external/amber_cuda/`, builds in a package/GPU-specific directory such as
`external/amber_cuda/build-pmemd26-sm_61-only`, installs to
`envs/amber_cuda/opt/pmemd26`, writes wrappers such as
`envs/amber_cuda/bin/pmemd.cuda`, and logs to `logs/install_amber_cuda.log`.
It auto-detects local GPU targets with `nvidia-smi`; override with
`ANKROS_AMBER_CUDA_GPU_TARGETS="sm_61 sm_86"` if needed. By default, those GPU
targets replace Amber's broad CUDA architecture list to reduce build memory and
time; set `ANKROS_AMBER_CUDA_ONLY_GPU_TARGETS=0` to keep Amber's defaults too.
The build defaults to two parallel jobs to reduce memory pressure; override
with `ANKROS_AMBER_CUDA_JOBS=4` or `setup_envs.sh --amber-cuda-jobs 4` only on
machines with enough RAM. It can also be run as part of setup:

```bash
tmux new -s ankros-setup 'AMBER_CUDA_PACKAGE=~/Downloads/pmemd26.tar.bz2 bash setup_envs.sh --with-amber-cuda --amber-cuda-jobs 2'
```

## Fresh Clone Quickstart

From a new clone, run commands from the repository root:

```bash
git clone <repo-url> ANKROS_enzyme
cd ANKROS_enzyme
mkdir -p logs
```

Create the main environment plus Amber CUDA. The Amber archive is licensed and
must stay outside git:

```bash
tmux new -s setup-envs 'AMBER_CUDA_PACKAGE=~/Downloads/pmemd26.tar.bz2 bash setup_envs.sh --with-amber-cuda --amber-cuda-jobs 12 2>&1 | tee logs/setup_envs.log'
```

If the holo model is missing, build or refresh it first:

```bash
bash scripts/structural/run_pipeline.sh
```

Prepare the standard pH 6.3, 500 mM NaCl water system:

```bash
tmux new -s md-water-setup 'bash scripts/experiments/fad_fmn_md_modeling/run_md_setup.sh --ph 6.3 2>&1 | tee logs/md_water_setup.log'
```

Prepare the pH 6.3 mixed-solvent system with 500 mM NaCl, 23% w/v sorbitol,
and 32% v/v glycerol:

```bash
tmux new -s md-mixed-setup 'bash scripts/experiments/fad_fmn_md_modeling/run_md_setup.sh --solvent mixed --ph 6.3 2>&1 | tee logs/md_mixed_setup.log'
```

Run quick Amber CUDA diagnostics after setup:

```bash
# Water, constant 400 K diagnostic
tmux new -s water-constant-test 'bash scripts/experiments/fad_fmn_md_modeling/06_run_amber_cuda_constant_temp.sh --diagnostic --temp-k 400 --prod-ns 0.01 --prod-segment-ns 0.01 --ig-seed 630401 --run-name water_pH6p3_constant_400K_quick 2>&1 | tee logs/water_pH6p3_constant_400K_quick.log'

# Water, quick ramp diagnostic
tmux new -s water-ramp-test 'bash scripts/experiments/fad_fmn_md_modeling/04_run_amber_cuda_juanma_ramp.sh --diagnostic --run-name water_pH6p3_quick_ramp --ig-seed 630501 2>&1 | tee logs/water_pH6p3_quick_ramp.log'

# Mixed solvent, constant 400 K diagnostic
tmux new -s mixed-constant-test 'bash scripts/experiments/fad_fmn_md_modeling/06_run_amber_cuda_constant_temp.sh --diagnostic --temp-k 400 --prod-ns 0.01 --prod-segment-ns 0.01 --ig-seed 630402 --run-name mixed_pH6p3_constant_400K_quick --system-dir results/experiments/fad_fmn_md_modeling/amber/systems/ankros_pH6p3_500mM_NaCl_23wv_sorbitol_32vv_glycerol 2>&1 | tee logs/mixed_pH6p3_constant_400K_quick.log'

# Mixed solvent, quick ramp diagnostic
tmux new -s mixed-ramp-test 'bash scripts/experiments/fad_fmn_md_modeling/04_run_amber_cuda_juanma_ramp.sh --diagnostic --run-name mixed_pH6p3_quick_ramp --ig-seed 630502 --system-dir results/experiments/fad_fmn_md_modeling/amber/systems/ankros_pH6p3_500mM_NaCl_23wv_sorbitol_32vv_glycerol 2>&1 | tee logs/mixed_pH6p3_quick_ramp.log'
```

All MD outputs are under:

```text
results/experiments/fad_fmn_md_modeling/
```

For each run, the most useful parseable files are:

- `protocol_manifest.json` for exact topology, coordinates, seeds, and stage plan.
- `restart_manifest.tsv` for each stage, input restart, output restart, trajectory, and status.
- `latest_restart.txt` for resuming from the newest completed checkpoint.
- `stage_metrics.tsv` for temperature, density, volume, energy, CA RMSD, and CA RMSF summaries.
- `metrics/<stage>.json` for per-stage machine-readable metrics.

Rerun the same MD command to resume: completed stages are skipped when their
Amber output contains `Total wall time` and the expected restart exists.

## Downloading external data

The active holo builder stages its required crystal inputs under
`results/structural/inputs/`:

```bash
bash scripts/structural/run_pipeline.sh --no-download
```

Use the MSA utility for AlphaFold DB homolog models. With no arguments it reads
the QC-filtered representative set and downloads every available AFDB model:

```bash
./envs/ankros/bin/python scripts/msa_OGT/13_download_afdb.py
```

Two data sources that require manual steps:

1. **AF3 server predictions.** Not available via public API. Submit the
   relevant FASTAs from `sequences/af3/` to the AlphaFold 3 server. The setup
   script looks for `folds_*.zip` under `structures/af3/` and extracts missing
   files in place. For a manual refresh:
   ```bash
   mkdir -p structures/af3
   unzip -nq structures/af3/folds_*.zip -d structures/af3
   ```
2. **UniRef database** (for `scripts/msa_OGT/01_mmseqs_search.py` only):
   ```bash
   mkdir -p ~/databases/uniref90
   wget -O ~/databases/uniref90/uniref90.fasta.gz \
     https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
   # Then point the pipeline at it:
   export UNIREF_DB=uniref90       # or uniref50 (default: uniref90)
   export UNIREF_DIR=~/databases   # default: ~/databases
   ```

## Running the pipelines

See the subdirectory READMEs for step-by-step instructions. Quick
orchestrators:

```bash
# MSA pipeline (requires UniRef DB)
./envs/ankros/bin/python scripts/msa_OGT/01_mmseqs_search.py
./envs/ankros/bin/python scripts/msa_OGT/02_extract_annotate.py
./envs/ankros/bin/python scripts/msa_OGT/03_classify_classI.py
./envs/ankros/bin/python scripts/msa_OGT/04_build_master_set.py
./envs/ankros/bin/python scripts/msa_OGT/05_align_hmm.py
./envs/ankros/bin/python scripts/msa_OGT/06_alignment_qc.py
./envs/ankros/bin/python scripts/msa_OGT/07_refine_linker.py
./envs/ankros/bin/python scripts/msa_OGT/08_concat_refined_alignment.py
./envs/ankros/bin/python scripts/msa_OGT/09_plot_domain_logos.py
./envs/ankros/bin/python scripts/msa_OGT/10_plot_conservation.py
./envs/ankros/bin/python scripts/msa_OGT/11_build_tree_iqtree.py
./envs/ankros/bin/python scripts/msa_OGT/12_call_regime_clades.py
./envs/ankros/bin/python scripts/msa_OGT/13_download_afdb.py
./envs/ankros/bin/python scripts/msa_OGT/14_validate_structural_matches.py
./envs/ankros/bin/python scripts/msa_OGT/15_compute_solubility_aggregability.py --workers 8
./envs/ankros/bin/python scripts/msa_OGT/16_plot_thermal_clade_logos.py
./envs/ankros/bin/python scripts/msa_OGT/17_classify_helix_extension.py
./envs/ankros/bin/python scripts/msa_OGT/18_plot_interactive_tree.py

# Holo model construction (requires crystals + AF3)
bash scripts/structural/run_pipeline.sh
```

## Further reading

- `AGENTS.md` — project-specific operating rules, domains, environments, and
  pipeline layout.
- `info/INDEX.md` — map of reference papers, enzyme notes, and Design-Expert
  parsing context.
- `scripts/README.md` — current pipeline dependency map.
- `structures/README.md` — structure-pool layout and provenance.
- `best_results.md` — current result snapshot and interpretation notes.
