# structural_evo_analysis

General MSA, tree, conservation, and structure-scoring pipeline for a single
query protein.

This is a generalized branch of the ANKros `scripts/msa_OGT/` workflow. It
keeps the original ANKros pipeline unchanged and writes default outputs under:

```text
results/structural_evo_analysis/
structures/structural_evo_analysis/
```

## Inputs

- A single-query protein FASTA.
- A local searchable sequence database FASTA, defaulting to
  `${UNIREF_DIR}/${SEA_DB}/${SEA_DB}.fasta.gz` with `SEA_DB=uniref90`.
- Optional user metadata TSV for tree clade annotation. The metadata must have
  an ID column matching tree tip IDs, usually `id`, and a trait column such as
  `temperature` or another continuous/categorical value.
- Optional query PDB for structure scoring. AFDB structures are downloaded for
  homolog accessions resolved from `repset_metadata.tsv`.

## Recommended Long Run

Full searches, IQ-TREE, downloads, and structure scoring can run for a long
time. Run the full pipeline inside `tmux` and log it:

```bash
tmux new -s structural-evo
mkdir -p logs
SEA_OUT_DIR=results/structural_evo_analysis \
  bash scripts/structural_evo_analysis/run_pipeline.sh sequences/query.fa \
  2>&1 | tee logs/structural_evo_analysis.log
```

Detach with `Ctrl-b d`. Monitor with:

```bash
tmux attach -t structural-evo
```

## Stepwise Use

```bash
./envs/ankros/bin/python scripts/structural_evo_analysis/01_mmseqs_search.py \
  --query sequences/query.fa \
  --out-dir results/structural_evo_analysis

./envs/ankros/bin/python scripts/structural_evo_analysis/02_align_mafft.py \
  --out-dir results/structural_evo_analysis

./envs/ankros/bin/python scripts/structural_evo_analysis/03_build_tree_iqtree.py \
  --out-dir results/structural_evo_analysis/tree

./envs/ankros/bin/python scripts/structural_evo_analysis/04_annotate_clades.py \
  --tree results/structural_evo_analysis/tree/query_msa.treefile \
  --metadata metadata.tsv \
  --trait-column temperature \
  --low-threshold 20 \
  --high-threshold 45

./envs/ankros/bin/python scripts/structural_evo_analysis/05_conserved_positions.py \
  --alignment results/structural_evo_analysis/repset_aligned.fa

./envs/ankros/bin/python scripts/structural_evo_analysis/06_download_afdb.py \
  --metadata results/structural_evo_analysis/repset_metadata.tsv

./envs/ankros/bin/python scripts/structural_evo_analysis/07_score_structures.py \
  --query-pdb structures/structural_evo_analysis/query.pdb
```

## Outputs

- `mmseqs_search_results.tsv`: raw MMseqs hit table.
- `hits_metadata.tsv`: searchable hit metadata and filtering status.
- `repset.fa`: query plus retained homologs.
- `repset_aligned.fa`: MAFFT MSA.
- `query_column_map.tsv`: alignment columns mapped to query positions.
- `tree/query_msa.treefile`: IQ-TREE output tree.
- `metadata_clades/clade_annotations.tsv`: per-clade trait summaries.
- `metadata_clades/called_clades.tsv`: clades passing category enrichment
  thresholds.
- `conservation/conserved_positions.tsv`: conserved query/alignment columns.
- `afdb_downloads/download_manifest.tsv`: AFDB per-accession status.
- `structure_scores/`: merged CamSol-style and Aggrescan3D score tables.

## Notes

- This pipeline does not use ANKros domain boundaries, class-I photolyase
  filtering, or OGTFinder unless the user provides equivalent metadata.
- Continuous metadata can be summarized without thresholds, or converted into
  low/mid/high categories with `--low-threshold` and `--high-threshold`.
- The CamSol implementation reused here is the repository's local
  CamSol-style scorer, not a claim of identity with the closed CamSol server.
