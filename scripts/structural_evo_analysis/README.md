# structural_evo_analysis

Maintained pipeline for protein vulnerability analysis: diverse MMseqs2 homolog
selection, OGT-aware MSA metadata, conservation, optional OGT enrichment,
AFDB/query-structure scoring, and residue-level vulnerability ranking.

Default outputs:

```text
~/structural_evo_analysis/results/
~/structural_evo_analysis/structures/
```

Override with `SEA_WORK_DIR`, or set `SEA_OUT_DIR` and `SEA_STRUCTURE_DIR`
independently.

The pipeline is portable when this repository is installed as a skill under a
directory such as `~/.codex/skills/structural-evo-analysis/`. Environments and
scripts remain under the skill root, while runtime results and structures
default to the user's home directory. `run_pipeline.sh` can be launched from
outside the skill directory.

## Inputs

- Single-sequence protein FASTA. The bundled test query is
  `test/photoHymenobact.fa`.
- Query PDB structure. This is required for vulnerability ranking. Provide it
  with `SEA_QUERY_PDB` or `--query-pdb`, or place it at
  `~/structural_evo_analysis/structures/query.pdb`. If a small real test
  structure is committed, `SEA_QUERY_PDB=test/query.pdb` is suitable for smoke
  tests. MSA-only and conservation-only runs do not require a PDB.
- Local protein database FASTA, usually UniRef. Configure with
  `SEA_DB_FASTA`/`SEA_DB_MMSEQS`, or with `UNIREF_DIR` and `SEA_DB`.
  UniRef shorthand values `50`, `90`, `100`, `uniref50`, `uniref90`, and
  `uniref100` are accepted by `SEA_DB` and step 01 `--db-name`.
  Ask for the database location before running. If it needs downloading, use
  `prepare_uniref_database.py --download` only after user approval.
- OGT metadata for OGT-aware annotations. Step 01 does not join OGT metadata by
  default; use `SEA_OGT_AWARE=1` or step 01 `--join-ogt` to add taxonomy-derived
  OGT metadata from `data/ogt_taxid_summary.tsv` when hit headers contain
  `TaxID=`.

## Full Example

Run long jobs in `tmux`:

```bash
mkdir -p logs
tmux new -s structural-evo \
  'SEA_QUERY_PDB=/path/to/query.pdb SEA_WORK_DIR=~/structural_evo_analysis bash scripts/structural_evo_analysis/run_pipeline.sh test/photoHymenobact.fa 2>&1 | tee logs/photoHymenobact_example.log'
```

For MSA-only work, use:

```bash
SEA_PIPELINE_MODE=msa \
  bash scripts/structural_evo_analysis/run_pipeline.sh test/photoHymenobact.fa
```

OGT enrichment is optional; omit the metadata and trait arguments when the goal
is diversity-first vulnerability ranking with no OGT join during search. Set
`SEA_OGT_AWARE=1` to join OGT metadata during step 01, run OGT clade enrichment,
and make logos/viewers clade-aware. For manual OGT clade enrichment, pass
`repset_metadata.tsv ogt` and set `SEA_LOW_THRESHOLD`/`SEA_HIGH_THRESHOLD`; for
OGT, the example thresholds use `ogt <= 20` as low and `ogt >= 45` as high.

Step 01 selects a diverse identity-stratified MSA subset by default
(`SEA_MAX_REPSET_SEQS=500`, including query). Agents may adjust MMseqs/search
or subset parameters when needed, but must record exact values and rationale.

Prepare or validate the UniRef database explicitly before step 01:

```bash
./envs/structural_evo/bin/python scripts/structural_evo_analysis/prepare_uniref_database.py \
  --db uniref90 \
  --uniref-dir ~/databases \
  --download \
  --create-mmseqs
```

## Stepwise Use

```bash
./envs/structural_evo/bin/python scripts/structural_evo_analysis/01_mmseqs_search.py \
  --query test/photoHymenobact.fa \
  --join-ogt \
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
  --dest ~/structural_evo_analysis/structures/afdb \
  --manifest results/photoHymenobact_example/afdb_downloads/download_manifest.tsv

./envs/structural_evo/bin/python scripts/structural_evo_analysis/07_score_structures.py \
  --afdb-dir ~/structural_evo_analysis/structures/afdb \
  --out-dir results/photoHymenobact_example/structure_scores

./envs/structural_evo/bin/python scripts/structural_evo_analysis/08_vulnerability_analysis.py \
  --conservation results/photoHymenobact_example/conservation/position_conservation.tsv \
  --scores results/photoHymenobact_example/structure_scores/per_residue_scores.tsv \
  --query-pdb ~/structural_evo_analysis/structures/query.pdb \
  --out-dir results/photoHymenobact_example/vulnerability

./envs/structural_evo/bin/python scripts/structural_evo_analysis/09_group_score_summary.py \
  --metadata results/photoHymenobact_example/repset_metadata.tsv \
  --scores results/photoHymenobact_example/structure_scores/global_scores.tsv \
  --query-pdb ~/structural_evo_analysis/structures/query.pdb \
  --out-dir results/photoHymenobact_example/vulnerability

./envs/structural_evo/bin/python scripts/structural_evo_analysis/10_sequence_logos.py \
  --alignment results/photoHymenobact_example/repset_aligned.fa \
  --group-members results/photoHymenobact_example/vulnerability/group_members.tsv \
  --out-dir results/photoHymenobact_example/logos

./envs/structural_evo/bin/python scripts/structural_evo_analysis/11_write_viewers.py \
  --alignment results/photoHymenobact_example/repset_aligned.fa \
  --tree results/photoHymenobact_example/tree/query_msa.treefile \
  --metadata results/photoHymenobact_example/repset_metadata.tsv \
  --out-dir results/photoHymenobact_example/viewers
```

Use `--skip-aggrescan3d` on step 07 when the optional Aggrescan3D environment
has not been installed.

## OGT Summary Table

`data/ogt_taxid_summary.tsv` is the maintained lookup table for OGT joins. It
contains one resolved row per covered taxid and preserves provenance columns.
No other files should be stored under `data/`. Regenerate it from an external
raw OGTFinder table with:

```bash
./envs/structural_evo/bin/python scripts/structural_evo_analysis/build_ogt_summary.py \
  --input /path/to/growth_temp_dataset_OGTFinder.tsv \
  --output data/ogt_taxid_summary.tsv
```

## Outputs

- `mmseqs_search_results.tsv`: raw MMseqs hit table.
- `hits_metadata.tsv`: hit metadata, filters, taxonomy, and optional OGT join
  fields.
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
- `vulnerability/top_vulnerable_positions.tsv`: highest-priority query
  structure residues by conservation/aggregation/solubility rank.
- `vulnerability/group_score_summary.tsv`: average CamSol/Aggrescan3D global
  structure scores for all representatives or OGT clades, compared with query.
- `logos/`: group image logos and consensus sequences.
- `viewers/tree.html`, `viewers/alignment.html`: self-contained HTML viewers.

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

- The optional OGT join depends on usable taxonomy in hit headers. For
  non-UniRef databases, provide a metadata TSV keyed by sequence ID or ensure
  headers include `TaxID=`.
- The CamSol implementation used by step 07 is a transparent local
  approximation, not the closed CamSol server.
- Aggrescan3D scoring requires `bash setup_envs.sh --with-aggrescan3d`.
