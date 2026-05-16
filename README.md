# Structural Evolution Analysis Skill

Reusable pipeline scaffold for single-protein structural evolutionary analysis.
The maintained workflow is `scripts/structural_evo_analysis/`: MMseqs2 homolog
search, MAFFT alignment, IQ-TREE phylogeny, metadata-driven clade annotation,
conservation scoring, AlphaFold DB download, and optional structure scoring.

The bundled example query is:

```text
sequences/photoHymenobact.fa
```

The committed metadata table in `data/growth_temp_dataset_OGTFinder.tsv` is used
to annotate temperature-associated clades when database hit headers include
`TaxID=` values, as UniRef headers do.

## Layout

```text
scripts/structural_evo_analysis/   main pipeline
scripts/msa_OGT/                   legacy/reference code; step 07 reuses scorers
config/                            minimal conda environment specs
data/                              metadata used for clade annotation
sequences/                         bundled example query FASTA
structures/                        downloaded or user-provided structures
results/                           generated pipeline outputs
logs/                              local run logs
```

## Setup

Install conda or mamba, then run from the repository root:

```bash
bash setup_envs.sh
```

This creates `envs/structural_evo` from `config/environment.yml` and smoke-tests
Biopython, MMseqs2, MAFFT, and IQ-TREE.

Optional Aggrescan3D support for structure scoring:

```bash
bash setup_envs.sh --with-aggrescan3d
```

Use the environment Python explicitly:

```bash
./envs/structural_evo/bin/python scripts/structural_evo_analysis/01_mmseqs_search.py --help
```

## External Database

The full pipeline needs a local searchable protein FASTA. By default it expects:

```text
${UNIREF_DIR}/${SEA_DB}/${SEA_DB}.fasta.gz
```

with `UNIREF_DIR=~/databases` and `SEA_DB=uniref90`. Override with:

```bash
export SEA_DB_FASTA=/path/to/protein_database.fa.gz
export SEA_DB_MMSEQS=/path/to/mmseqs_database_prefix
```

Step 01 creates the MMseqs database prefix if it is missing.

## Example Run

Full searches, tree inference, downloads, and scoring can be long-running. Use
`tmux` and log the run:

```bash
mkdir -p logs
tmux new -s structural-evo \
  'SEA_OUT_DIR=results/photoHymenobact_example SEA_LOW_THRESHOLD=20 SEA_HIGH_THRESHOLD=45 bash scripts/structural_evo_analysis/run_pipeline.sh sequences/photoHymenobact.fa results/photoHymenobact_example/repset_metadata.tsv ogt 2>&1 | tee logs/photoHymenobact_example.log'
```

That command uses the bundled query, joins the committed OGT metadata during
step 01, and asks step 04 to call low/mid/high clades from the `ogt` column
using 20 C and 45 C thresholds.

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
```

## Notes For Skill Authors

- Keep the skill focused on `scripts/structural_evo_analysis/`.
- Load `scripts/structural_evo_analysis/README.md` for step-level details.
- Use `sequences/photoHymenobact.fa` as the reproducible example query.
- Use `data/growth_temp_dataset_OGTFinder.tsv` only as metadata; do not treat it
  as sequence input.
- Prefer explicit environment variables and output directories so runs are
  reproducible and do not overwrite each other.
