# Reference

## Maintained Modules

The maintained scripts live under `scripts/structural_evo_analysis/`.

| Step | Script | Purpose | Main Outputs |
| --- | --- | --- | --- |
| 01 | `01_mmseqs_search.py` | MMseqs search, post-search filters, diversity-first representative set, optional OGT join | `mmseqs_search_results.tsv`, `hits_metadata.tsv`, `repset.fa`, `repset_metadata.tsv` |
| 02 | `02_align_mafft.py` | Align query plus representatives | `repset_aligned.fa`, `query_column_map.tsv`, `alignment_summary.tsv` |
| 03 | `03_build_tree_iqtree.py` | Build IQ-TREE tree from the MSA | `tree/query_msa.treefile`, `tree/run_manifest.txt` |
| 04 | `04_annotate_clades.py` | Annotate trait-enriched clades such as OGT regimes | `metadata_clades/` |
| 05 | `05_conserved_positions.py` | Score conserved query positions from an MSA | `conservation/position_conservation.tsv`, `conservation/conserved_positions.tsv` |
| 06 | `06_download_afdb.py` | Download/cache AFDB PDBs for representative accessions | `afdb_downloads/download_manifest.tsv`, `structures/afdb/` |
| 07 | `07_score_structures.py` | Score query and AFDB structures with CamSol-style and optional Aggrescan3D scores | `structure_scores/` |
| 08 | `08_vulnerability_analysis.py` | Rank query residues from conservation and per-residue structure scores | `vulnerability/top_vulnerable_positions.tsv` |
| 09 | `09_group_score_summary.py` | Summarize global structure scores by representatives or clades | `vulnerability/group_score_summary.tsv` |
| 10 | `10_sequence_logos.py` | Make sequence-logo summaries for groups | `logos/` |
| 11 | `11_write_viewers.py` | Write self-contained tree/alignment HTML viewers | `viewers/` |

## Wrapper Modes

`run_pipeline.sh` supports `SEA_PIPELINE_MODE`:

- `full`: default end-to-end vulnerability workflow. Requires a query PDB.
- `msa`: runs steps 01 and 02 only. Does not require a query PDB.
- `conservation`: runs steps 01 to 05. Does not require a query PDB.

Examples:

```bash
SEA_PIPELINE_MODE=msa \
bash scripts/structural_evo_analysis/run_pipeline.sh test/photoHymenobact.fa
```

```bash
SEA_PIPELINE_MODE=conservation \
bash scripts/structural_evo_analysis/run_pipeline.sh test/photoHymenobact.fa
```

```bash
SEA_QUERY_PDB=/path/to/query.pdb \
bash scripts/structural_evo_analysis/run_pipeline.sh test/photoHymenobact.fa
```

## OGT And Non-OGT Behavior

Non-OGT runs do not join OGT metadata during search. Step 01 prioritizes the
diverse query-identity-stratified representative set.

For OGT clade context without structure scoring, use:

```bash
SEA_PIPELINE_MODE=conservation SEA_OGT_AWARE=1 \
bash scripts/structural_evo_analysis/run_pipeline.sh test/photoHymenobact.fa
```

For a full OGT-aware vulnerability run, also provide `SEA_QUERY_PDB`.

For manual use, pass `--join-ogt` to step 01 before running clade annotation.

## Required Inputs By Recipe

- Diversity MSA: one protein FASTA and a local searchable protein database.
- Conservation: an aligned FASTA, or the diversity MSA inputs if using the
  wrapper.
- OGT context: OGT-aware metadata, usually `data/ogt_taxid_summary.tsv`, and
  hit headers with usable `TaxID=` values.
- Structure scoring: query PDB and optional AFDB/query-related PDB directory.
- Full vulnerability ranking: query FASTA, query PDB, local protein database,
  and structure-scoring dependencies.

## Environment Variables

- `SEA_WORK_DIR`: base runtime directory. Defaults to
  `~/structural_evo_analysis`.
- `SEA_OUT_DIR`: result directory. Defaults to `$SEA_WORK_DIR/results`.
- `SEA_STRUCTURE_DIR`: structure/cache directory. Defaults to
  `$SEA_WORK_DIR/structures`.
- `SEA_PIPELINE_MODE`: `full`, `msa`, or `conservation`.
- `SEA_QUERY_PDB`: query PDB for full vulnerability runs.
- `SEA_DB`, `UNIREF_DIR`, `SEA_DB_FASTA`, `SEA_DB_MMSEQS`: database selection.
- `SEA_OGT_AWARE`: set to `1` to join OGT metadata and run OGT clade context.
- `SEA_MAX_REPSET_SEQS`: representative-set cap including query. Defaults to
  `500`; use `0` only when all filtered hits are required.
- `SEA_THREADS`: thread count for scripts that support it.

## Long-Running Runs

Run full searches, tree inference, AFDB downloads, and structure scoring in
`tmux` with logs:

```bash
mkdir -p logs
tmux new -s structural-evo \
  'SEA_QUERY_PDB=/path/to/query.pdb bash scripts/structural_evo_analysis/run_pipeline.sh test/photoHymenobact.fa 2>&1 | tee logs/structural_evo.log'
```

Record the command, working directory, inputs, outputs/logs, environment
variables, and how to monitor or stop the task.
