# structural_evo_analysis

Maintained pipeline for applying MSA, phylogeny, clade annotation,
conservation, AFDB download, and optional structure scoring to one query
protein.

Default outputs:

```text
results/structural_evo_analysis/
structures/structural_evo_analysis/
```

Override with `SEA_OUT_DIR` and `SEA_STRUCTURE_DIR`.

## Inputs

- Single-sequence protein FASTA. The bundled test query is
  `sequences/photoHymenobact.fa`.
- Local protein database FASTA, usually UniRef. Configure with
  `SEA_DB_FASTA`/`SEA_DB_MMSEQS`, or with `UNIREF_DIR` and `SEA_DB`.
- Optional metadata for tree annotation. Step 01 automatically adds
  taxonomy-derived OGT metadata from `data/growth_temp_dataset_OGTFinder.tsv`
  when hit headers contain `TaxID=`.
- Optional query PDB at `structures/structural_evo_analysis/query.pdb`.

## Full Example

Run long jobs in `tmux`:

```bash
mkdir -p logs
tmux new -s structural-evo \
  'SEA_OUT_DIR=results/photoHymenobact_example SEA_LOW_THRESHOLD=20 SEA_HIGH_THRESHOLD=45 bash scripts/structural_evo_analysis/run_pipeline.sh sequences/photoHymenobact.fa results/photoHymenobact_example/repset_metadata.tsv ogt 2>&1 | tee logs/photoHymenobact_example.log'
```

`SEA_LOW_THRESHOLD` and `SEA_HIGH_THRESHOLD` are passed to the clade annotator.
For OGT, the example uses `ogt <= 20` as low, `ogt >= 45` as high, and the
middle range as mid.

## Stepwise Use

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

./envs/structural_evo/bin/python scripts/structural_evo_analysis/05_conserved_positions.py \
  --alignment results/photoHymenobact_example/repset_aligned.fa \
  --out-dir results/photoHymenobact_example/conservation

./envs/structural_evo/bin/python scripts/structural_evo_analysis/06_download_afdb.py \
  --metadata results/photoHymenobact_example/repset_metadata.tsv \
  --dest structures/structural_evo_analysis/afdb \
  --manifest results/photoHymenobact_example/afdb_downloads/download_manifest.tsv

./envs/structural_evo/bin/python scripts/structural_evo_analysis/07_score_structures.py \
  --afdb-dir structures/structural_evo_analysis/afdb \
  --out-dir results/photoHymenobact_example/structure_scores
```

Use `--skip-aggrescan3d` on step 07 when the optional Aggrescan3D environment
has not been installed.

## Outputs

- `mmseqs_search_results.tsv`: raw MMseqs hit table.
- `hits_metadata.tsv`: hit metadata, filters, taxonomy, and OGT join status.
- `repset.fa`: query plus retained homologs.
- `repset_metadata.tsv`: metadata keyed by tree-tip IDs.
- `repset_aligned.fa`: MAFFT MSA.
- `query_column_map.tsv`: alignment columns mapped to query positions.
- `tree/query_msa.treefile`: IQ-TREE output tree.
- `metadata_clades/clade_annotations.tsv`: per-clade trait summaries.
- `metadata_clades/called_clades.tsv`: clades passing category enrichment
  thresholds.
- `conservation/conserved_positions.tsv`: conserved query/alignment columns.
- `afdb_downloads/download_manifest.tsv`: AFDB per-accession status.
- `structure_scores/`: merged CamSol-style and Aggrescan3D score tables.

## Parameters To Treat Carefully

- MMseqs defaults: `SEA_MMSEQS_S`, `SEA_MMSEQS_E`, `SEA_MIN_SEQ_ID`,
  `SEA_COVERAGE`, `SEA_MAX_SEQS`.
- Post-search filters: `SEA_MIN_LENGTH`, `SEA_MAX_LENGTH`,
  `SEA_POST_MIN_QCOV`, `SEA_POST_MIN_IDENTITY`, `SEA_POST_MAX_IDENTITY`.
- Tree thresholds: `--min-column-occupancy`, `--model`, `--bootstrap`,
  `--seed`.
- Clade thresholds: `--low-threshold`, `--high-threshold`,
  `--min-labelled`, `--min-fraction`.

Do not change these silently. Record overrides in logs or run manifests.

## Limitations

- The OGT join depends on usable taxonomy in hit headers. For non-UniRef
  databases, provide a metadata TSV keyed by sequence ID or ensure headers
  include `TaxID=`.
- The CamSol implementation reused by step 07 is a transparent local
  approximation, not the closed CamSol server.
- Aggrescan3D scoring requires `bash setup_envs.sh --with-aggrescan3d`.
