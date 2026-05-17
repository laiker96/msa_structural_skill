# Protein Vulnerability Analysis Skill

Reusable pipeline scaffold for single-protein vulnerability analysis. The
maintained workflow is `scripts/structural_evo_analysis/`: MMseqs2 homolog
search, diverse OGT-aware MSA, conservation scoring, optional OGT enrichment,
AlphaFold DB download, CamSol-style/Aggrescan3D structure scoring, and
residue-level vulnerability ranking.

The bundled example query is:

```text
sequences/photoHymenobact.fa
```

A query PDB is required. Provide it with `SEA_QUERY_PDB=/path/to/query.pdb` or
place it at `~/structural_evo_analysis/structures/query.pdb`.

The only bundled metadata table is `data/ogt_taxid_summary.tsv`, a compact
taxid-keyed OGT table used when database hit headers include `TaxID=` values,
as UniRef headers do.

## Layout

```text
scripts/structural_evo_analysis/   main self-contained pipeline
scripts/msa_OGT/                   legacy/reference code with scorer wrappers
config/                            minimal conda environment specs
data/                              OGT metadata used for clade annotation
sequences/                         bundled example query FASTA
structures/                        downloaded or user-provided structures
results/                           generated pipeline outputs
logs/                              local run logs
```

By default, runtime outputs are not written into the skill directory. They go
under:

```text
~/structural_evo_analysis/results/
~/structural_evo_analysis/structures/
```

Set `SEA_WORK_DIR` to change that base directory, or set `SEA_OUT_DIR` and
`SEA_STRUCTURE_DIR` independently.

## Setup

Install conda, mamba, or micromamba, then run from the repository root:

```bash
bash setup_envs.sh
```

This creates `envs/structural_evo` from `config/environment.yml` and smoke-tests
Biopython, MMseqs2, MAFFT, IQ-TREE, and matplotlib.

Optional Aggrescan3D support for structure scoring:

```bash
bash setup_envs.sh --with-aggrescan3d
```

Use the environment Python explicitly:

```bash
./envs/structural_evo/bin/python scripts/structural_evo_analysis/01_mmseqs_search.py --help
```

## Skill Installation

This directory can be copied under a Codex skills directory such as
`~/.codex/skills/structural-evo-analysis/`. The scripts keep environments under
the installed skill root:

```text
<skill-root>/envs/structural_evo/
<skill-root>/envs/aggrescan3d/
```

When running as a skill, check for those envs first. If they are missing, run
`bash setup_envs.sh` from the skill root, and use `tmux` plus a log for the
install because conda solves can be slow.

Before running an analysis as a skill, the agent should ask whether to use the
default `~/structural_evo_analysis` work directory or a user-selected directory.

## External Database

The full pipeline needs a local searchable protein FASTA. By default it expects:

```text
${UNIREF_DIR}/${SEA_DB}/${SEA_DB}.fasta.gz
```

with `UNIREF_DIR=~/databases` and `SEA_DB=uniref90`. `SEA_DB` and `--db-name`
also accept UniRef shorthand values such as `50`, `90`, `100`, `uniref50`,
`uniref90`, or `uniref100`.

The agent should ask for the database location before running step 01. If the
database is not already present, ask whether the user wants to download UniRef
and which identity level to use. Downloads can be large and should run in
`tmux` with a log. The helper below validates an existing location, or downloads
from UniProt only when `--download` is explicitly provided:

```bash
./envs/structural_evo/bin/python scripts/structural_evo_analysis/prepare_uniref_database.py \
  --db uniref90 \
  --uniref-dir ~/databases \
  --download \
  --create-mmseqs
```

Examples:

```bash
export SEA_DB=uniref50
./envs/structural_evo/bin/python scripts/structural_evo_analysis/01_mmseqs_search.py \
  --query sequences/photoHymenobact.fa \
  --out-dir results/photoHymenobact_uniref50

./envs/structural_evo/bin/python scripts/structural_evo_analysis/01_mmseqs_search.py \
  --query sequences/photoHymenobact.fa \
  --db-name 100 \
  --out-dir results/photoHymenobact_uniref100
```

Override exact paths with:

```bash
export SEA_DB_FASTA=/path/to/protein_database.fa.gz
export SEA_DB_MMSEQS=/path/to/mmseqs_database_prefix
```

If `SEA_DB_MMSEQS`/`--db-mmseqs` is omitted, step 01 creates a prefix beside
the selected FASTA using `<fasta-stem>_db`.

## Example Run

Full searches, tree inference, downloads, and scoring can be long-running. Use
`tmux` and log the run:

```bash
mkdir -p logs
tmux new -s structural-evo \
  'SEA_QUERY_PDB=/path/to/query.pdb SEA_WORK_DIR=~/structural_evo_analysis bash scripts/structural_evo_analysis/run_pipeline.sh sequences/photoHymenobact.fa 2>&1 | tee logs/photoHymenobact_example.log'
```

That command uses the bundled query, joins the committed OGT metadata during
step 01, builds a diverse MSA subset, scores available structures, and writes
the vulnerability report. OGT enrichment is optional; to call low/mid/high OGT
clades, pass `repset_metadata.tsv ogt` and set thresholds:

```bash
SEA_LOW_THRESHOLD=20 SEA_HIGH_THRESHOLD=45 \
  bash scripts/structural_evo_analysis/run_pipeline.sh \
  sequences/photoHymenobact.fa results/photoHymenobact_example/repset_metadata.tsv ogt
```

By default, step 01 builds the MSA from a diverse query-identity-stratified
subset of filtered MMseqs2 hits, capped by `SEA_MAX_REPSET_SEQS=500` including
the query. Set `SEA_MAX_REPSET_SEQS=0` only when all filtered hits are required.
The agent may adjust MMseqs/search/subset parameters to obtain a usable homolog
set, but those changes must be logged with exact values and rationale.

OGT-aware mode is off by default. Enable clade calling, clade-aware logos, and
tree/alignment clade coloring with `SEA_OGT_AWARE=1`.

## OGT Metadata

Use `data/ogt_taxid_summary.tsv` for agent-facing inspection and default
pipeline joins. It has one resolved OGT row per covered NCBI taxid, with
taxid/name/species/OGT/regime first, followed by provenance fields such as
source table rows, raw temperatures, source IDs, and whether the value came
from exact-taxid or species-level fallback.

No other files should be stored in `data/`. If regenerating the OGT summary from
an external raw OGTFinder table, provide that table explicitly:

```bash
./envs/structural_evo/bin/python scripts/structural_evo_analysis/build_ogt_summary.py \
  --input /path/to/growth_temp_dataset_OGTFinder.tsv \
  --output data/ogt_taxid_summary.tsv
```

Run steps manually when you need tighter control:

```bash
./envs/structural_evo/bin/python scripts/structural_evo_analysis/01_mmseqs_search.py \
  --query sequences/photoHymenobact.fa \
  --out-dir results/photoHymenobact_example

./envs/structural_evo/bin/python scripts/structural_evo_analysis/02_align_mafft.py \
  --out-dir results/photoHymenobact_example

./envs/structural_evo/bin/python scripts/structural_evo_analysis/03_build_tree_iqtree.py \
  --out-dir results/photoHymenobact_example/tree

./envs/structural_evo/bin/python scripts/structural_evo_analysis/04_annotate_clades.py \
  --tree results/photoHymenobact_example/tree/query_msa.treefile \
  --metadata results/photoHymenobact_example/repset_metadata.tsv \
  --trait-column ogt \
  --low-threshold 20 \
  --high-threshold 45
```

## Key Outputs

```text
results/<run>/mmseqs_search_results.tsv
results/<run>/hits_metadata.tsv
results/<run>/repset.fa
results/<run>/repset_metadata.tsv
results/<run>/repset_aligned.fa
results/<run>/query_column_map.tsv
results/<run>/tree/query_msa.treefile
results/<run>/metadata_clades/clade_annotations.tsv
results/<run>/metadata_clades/called_clades.tsv
results/<run>/conservation/conserved_positions.tsv
results/<run>/afdb_downloads/download_manifest.tsv
results/<run>/structure_scores/
results/<run>/vulnerability/top_vulnerable_positions.tsv
results/<run>/vulnerability/group_score_summary.tsv
results/<run>/logos/
results/<run>/viewers/tree.html
results/<run>/viewers/alignment.html
```

## Notes For Skill Authors

- Keep the skill focused on `scripts/structural_evo_analysis/`.
- Load `scripts/structural_evo_analysis/README.md` for step-level details.
- Use `sequences/photoHymenobact.fa` as the reproducible example query.
- Use `data/ogt_taxid_summary.tsv` only as metadata; do not treat it as
  sequence input.
- Prefer explicit environment variables and output directories so runs are
  reproducible and do not overwrite each other.
