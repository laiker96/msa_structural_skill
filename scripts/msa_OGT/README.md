# msa_OGT

OGT-aware UniRef MSA pipeline.

The core contract is:

1. Search UniRef once with MMseqs2.
2. Preserve the UniRef header and `TaxID=` for every match.
3. Join OGTFinder by taxid before downstream filtering.
4. Keep one master homolog FASTA/TSV.
5. Align a representative set to a curated class-I profile HMM.

No direct UniProt top-up is part of this pipeline.

## Steps

```bash
python3 scripts/msa_OGT/01_mmseqs_search.py
python3 scripts/msa_OGT/02_extract_annotate.py
python3 scripts/msa_OGT/03_classify_classI.py
python3 scripts/msa_OGT/04_build_master_set.py
python3 scripts/msa_OGT/05_align_hmm.py
python3 scripts/msa_OGT/06_alignment_qc.py
python3 scripts/msa_OGT/07_refine_linker.py
python3 scripts/msa_OGT/08_concat_refined_alignment.py
python3 scripts/msa_OGT/09_plot_domain_logos.py
python3 scripts/msa_OGT/10_plot_conservation.py
python3 scripts/msa_OGT/11_build_tree_iqtree.py
python3 scripts/msa_OGT/12_call_regime_clades.py
python3 scripts/msa_OGT/13_download_afdb.py
python3 scripts/msa_OGT/14_validate_structural_matches.py
python3 scripts/msa_OGT/15_compute_solubility_aggregability.py
python3 scripts/msa_OGT/16_plot_thermal_clade_logos.py
python3 scripts/msa_OGT/17_classify_helix_extension.py
python3 scripts/msa_OGT/18_plot_interactive_tree.py
python3 scripts/msa_OGT/17_test_helix_extension_adaptation.py
python3 scripts/msa_OGT/18_score_branch_rosetta_stability.py
python3 scripts/msa_OGT/19_visualize_n59_n60_alignment.py
python3 scripts/msa_OGT/20_prepare_dnds_sequence_inputs.py
python3 scripts/msa_OGT/21_fetch_dnds_cds.py
python3 scripts/msa_OGT/22_build_dnds_codon_alignment.py
```

Outputs live under `results/msa_OGT/`.

`13_download_afdb.py` is the MSA pipeline utility for AlphaFold DB homolog
structures. With no arguments it reads
`results/msa_OGT/repset_metadata_qc.tsv`, skips the ANKros query, and downloads
all available representative-set UniProt models into `structures/afdb/`.
Per-accession status is written to
`results/msa_OGT/afdb_downloads/download_manifest.tsv`. Use `--pre-qc-repset`
to target `repset_metadata.tsv` instead.

Step 14 runs Foldseek in TM-align mode to validate that local AFDB structures
match the ANKros class-I fold before structure scoring. It writes
`results/msa_OGT/structure_validation/quick_classI_mismatch_review.tsv`.
Hard `fail_possible_mismatch` rows are excluded from structure scoring by
default; borderline review rows are retained.

Step 15 runs the local CamSol-style and Aggrescan3D scorer utilities in this
directory for ANKros plus structurally validated AFDB structures. Intermediate
scorer outputs and the merged score tables are written under
`results/msa_OGT/structure_scores/`. Use `--workers N` to score independent
structures in parallel.
Step 16 reads that merged per-residue table by default when it exists.
Its HTML and `column_stats.tsv` also include secondary-structure states
(`H`, `E`, `C`) for ANKros and the selected representative structure in each
displayed clade. These use `mkdssp`/DSSP when available, with Biotite P-SEA as
a fallback.

Step 17 classifies AFDB-backed representatives by whether the ANKros core helix
around qpos 211-220 is extended upstream. It uses the refined alignment to map
ANKros qpos 192-210 plus intervening insertion columns as the candidate
extension region, then scores secondary structure from each structure model.
The binary `extension_signal=positive` is based on a substantial upstream helix
run (`extension_run_len >= 6`), so `no_core_helix` entries with that upstream
run are counted as positive. Outputs are written under
`results/msa_OGT/helix_extension_classifier/`. Run Step 18 after Step 17 so it
reads these calls and overlays helix-extension rings on the interactive tree.

Regime-comparison analyses that consume this pipeline's final alignment,
metadata, and tree live under `scripts/experiments/regime_comparison/`.

For a small/cheap UniRef50 smoke test without touching the default output:

```bash
UNIREF_DB=uniref50 MSA_OGT_OUT_DIR=results/msa_OGT_uniref50 \
  python3 scripts/msa_OGT/01_mmseqs_search.py
```

## Search Defaults

The default search is intentionally broad:

```text
-s 7.5
-e 1e-5
--min-seq-id 0.05
-c 0.50
--cov-mode 1
--max-seqs 50000
```

These can be overridden with environment variables:

```text
MSA_OGT_MMSEQS_S
MSA_OGT_MMSEQS_E
MSA_OGT_MIN_SEQ_ID
MSA_OGT_COVERAGE
MSA_OGT_MAX_SEQS
MSA_OGT_OUT_DIR
```

## Representative Selection

Step 05 uses MMseqs clustering to define sequence neighborhoods, but it does
not blindly accept MMseqs' chosen representative. For each cluster, the pipeline
chooses the best member using `master_homologs.tsv` metadata.

Representative priority is:

```text
has OGTFinder OGT
exact taxid match > species taxid > exact name > binomial name
point OGT > range-midpoint OGT
class-I bitscore
query coverage
search bitscore
```

The optional master precluster is controlled by `MSA_OGT_PRECLUSTER_ID`.
Defaults:

```text
UniRef50: disabled
UniRef90: 0.95
```

When enabled, precluster membership is written to
`master_precluster_members.tsv`. Diversity-cluster membership for meso, thermo,
and unlabelled top-up pools is written to `mmseqs_*_members.tsv`.

## OGTFinder Metadata

`hits_metadata.tsv`, `hits_filtered_metadata.tsv`, `master_homologs.tsv`, and
`repset_metadata.tsv` include OGTFinder provenance columns:

```text
taxid
ogt
regime
ogt_match_type
ogt_taxid
ogt_species_id
ogt_species
ogt_raw_temps
ogt_sources
ogt_types
ogt_source_ids
ogt_parse_modes
ogt_has_range
ogt_range_count
ogt_row_count
```

Range-valued temperatures are preserved in `ogt_raw_temps`. The numeric `ogt`
uses the range midpoint for range rows, then the median across all matching
optimum rows. For example, `15.0-35.0` becomes `25.0` for regime assignment,
while `ogt_raw_temps=15.0-35.0` and `ogt_parse_modes=range_midpoint` preserve
the original source form.

## Alignment QC

Step 06 scores the HMM-aligned representative set in ANKros coordinates using
the project domain boundaries:

```text
antenna:   1-130
linker:    131-205
catalytic: 206-437
```

It removes clear sequence outliers iteratively, but uses a permissive linker
coverage threshold because the linker is the indel-rich, hardest-to-align
region. Defaults can be overridden with:

```text
MSA_OGT_QC_MIN_GLOBAL_COV=0.75
MSA_OGT_QC_MIN_ANTENNA_COV=0.70
MSA_OGT_QC_MIN_LINKER_COV=0.35
MSA_OGT_QC_MIN_CATALYTIC_COV=0.80
MSA_OGT_QC_MAX_ITER=5
```

The PNGs show the ANKros query residue under each query-position column. Main
outputs:

```text
alignment_qc_metrics.tsv
alignment_qc_domain_metrics.tsv
alignment_qc_iterations.tsv
alignment_qc_rejected.tsv
repset_hmmalign_matchcols_qc.fa
repset_hmmalign_qc.fa
repset_metadata_qc.tsv
```

## Linker Refinement

Step 07 keeps the HMM/QC alignment as the ANKros coordinate anchor, extracts a
linker-centered window, fixes conserved anchor islands, and realigns only the
variable intervals with MAFFT E-INS-i:

```bash
python3 scripts/msa_OGT/07_refine_linker.py
```

Defaults:

```text
input:   repset_hmmalign_qc.fa
window:  ANKros 120-215
linker:  ANKros 131-205
anchors: auto-discovered from occupancy/conservation
method:  mafft --genafpair --maxiterate 1000
```

The default preserves all residues in the extracted window. For a compact
core-linker view that prevents rare long insertions from widening every block,
use `--robust`; the sequence metrics report any residue loss from that
projection.

Useful scope options:

```bash
python3 scripts/msa_OGT/07_refine_linker.py --scope annotated
python3 scripts/msa_OGT/07_refine_linker.py --scope meso
python3 scripts/msa_OGT/07_refine_linker.py --scope thermo
```

Main outputs under `linker_refined/`:

```text
linker_refined.fa
linker_refined_column_map.tsv
linker_refined_block_stats.tsv
linker_refined_anchor_metrics.tsv
linker_refined_sequence_metrics.tsv
linker_refined_summary.tsv
```

## Refined Full-Domain Alignment

Step 08 concatenates the HMM-anchored non-linker domains with the refined
linker-window alignment:

```text
HMM/QC columns:        ANKros 1-119
refined linker window: ANKros 120-215
HMM/QC columns:        ANKros 216-437
```

The wider 120-215 splice keeps the linker boundary context from step 07 instead
of cutting through insert-rich boundary blocks.

Main outputs:

```text
repset_hmmalign_linker_refined.fa
repset_hmmalign_linker_refined_column_map.tsv
repset_hmmalign_linker_refined_summary.tsv
```

## Domain Logos

Step 09 renders gap-weighted amino-acid sequence logos for the antenna, linker,
and catalytic domains from the refined full-domain alignment:

```bash
python3 scripts/msa_OGT/09_plot_domain_logos.py
```

Each logo includes a grey occupancy track. By default, the logos use only
query-position columns, so the linker is plotted as the 75 ANKros linker
positions. To inspect insertion columns fully contained inside each domain:

```bash
python3 scripts/msa_OGT/09_plot_domain_logos.py --include-contained-inserts \
  --out-dir results/msa_OGT/figures/domain_logos_with_inserts
```

Main outputs:

```text
figures/domain_logos/{antenna,linker,catalytic}_logo.png
figures/domain_logos/{antenna,linker,catalytic}_logo_columns.tsv
figures/domain_logos/domain_logo_summary.tsv
```

## Conservation Profile

Step 10 plots the per-column conservation score on the ANKros coordinate
frame, using only query-position columns from the refined full-domain
alignment:

```bash
python3 scripts/msa_OGT/10_plot_conservation.py
```

The curve is the fraction of non-gap residues matching the column consensus.
A light occupancy fill is drawn underneath for context.

Main outputs:

```text
figures/conservation_profile/conservation_profile.png
figures/conservation_profile/conservation_profile.tsv
```

## Final Tree And Clades

Step 11 builds the final IQ-TREE phylogeny from the refined full-domain
alignment after excluding the indel-rich linker. By default it uses ANKros
query-position columns from the antenna and catalytic domains and runs full ML
mode with `LG+F+R4` and UFBoot:

```bash
python3 scripts/msa_OGT/11_build_tree_iqtree.py
```

Default tree input:

```text
antenna:   ANKros 1-130
catalytic: ANKros 206-437
model:     LG+F+R4
bootstrap: 1000 UFBoot
seed:      42
```

By default, `sequences/outgroup.fa` is added to the tree MSA with
`mafft --add --keeplength`, and the resulting outgroup IDs are passed to
IQ-TREE with `-o` for rooted output. The outgroup panel is restricted to
verified class II CPD photolyases from `classII_known.fa`; the more distant
6-4 photolyase / cryptochrome negatives remain classifier exclusions, not tree
rooting taxa. Use `--no-outgroup` for an unrooted tree.

Use `--fast` for the older quick `LG -fast` topology, or `--prepare-only` to
write the tree MSA without launching IQ-TREE.

Main outputs:

```text
tree_antenna_catalytic/antenna_catalytic.fa
tree_antenna_catalytic/antenna_catalytic_ingroup.fa
tree_antenna_catalytic/antenna_catalytic_column_map.tsv
tree_antenna_catalytic/antenna_catalytic_outgroup_ids.txt
tree_antenna_catalytic/antenna_catalytic_summary.tsv
tree_antenna_catalytic/ankros_antenna_catalytic.treefile
tree_antenna_catalytic/run_manifest.txt
ankros_antenna_catalytic.nw
```

Step 12 calls thermal-regime clades on the final tree. The clade-calling
parameters are command-line options:

```bash
python3 scripts/msa_OGT/12_call_regime_clades.py \
  --min-support 90 \
  --min-labelled 3 \
  --min-target-fraction 1.0 \
  --max-unlabelled-fraction 1.0
```

Missing branch supports are allowed by default so clades can still be called
from fast or support-free treefiles. Add `--require-support` when using a full
IQ-TREE output where support filtering should be enforced.

The regime call is based only on labelled descendants. Unlabelled/unknown tips
inside the same subtree are retained as phylogenetic context but do not dilute
`target_fraction`; control them with `--max-unlabelled-fraction`. For example,
a subtree with 3 meso labelled tips and 9 unknown tips is a meso clade by
default because the labelled fraction is 3/3.

When local AFDB structures are present, step 12 also appends within-clade
domain geometry summaries to `regime_clades.tsv`: number of usable structures,
pair count, and mean/median C-alpha RMSD for antenna, linker, and catalytic
domains. Use `--skip-domain-geometry` to omit this structural pass.

Main outputs:

```text
regime_clades/regime_clades.tsv
regime_clades/regime_clades_tip_metadata.tsv
regime_clades/regime_clades_summary.txt
regime_clades/regime_clades_rooted_ingroup.nwk
```

Step 18 writes a rectangular interactive HTML tree. Tips are colored by regime,
all internal nodes are named as `N*`, and called clades are marked at their
MRCA nodes. Run this after Step 17 when helix-extension calls are available so
the tree includes the helix-extension overlay. The radial layout from the
exploratory viewer has been removed.

```bash
python3 scripts/msa_OGT/18_plot_interactive_tree.py
```

Main outputs:

```text
interactive_tree/ogt_msa_tree.html
interactive_tree/ogt_msa_tree_node_summary.tsv
interactive_tree/ogt_msa_tree_named_nodes.nwk
```

## Thermal Clade Logos

Step 16 generates ANKros-positioned consensus logos for thermal regime clades
from the refined full-domain MSA and `regime_clades/regime_clades.tsv`.
Residue classes use ANKros structure evidence and local AFDB structures when
available.

```bash
python3 scripts/msa_OGT/16_plot_thermal_clade_logos.py
```

By default it plots the two largest labelled psychro, meso, and thermo clades.
Use `--max-clades-per-regime 0` to include all clades, `--html` for the
interactive single-table view, or `--merge-regime-clades` to build one
consensus/logo per thermal regime from the union of all MSA sequences assigned
to that regime.

In HTML mode, insert columns compact dynamically based on the currently visible
clades and expand again when a visible clade has residues in that column. The
HTML also includes a top-row total-MSA conservation track aligned to the ANKros
positioned sequence columns.

Main outputs:

```text
figures/thermal_clade_consensus_logos_representatives/by_clade/<clade>/<domain>_logo.png
figures/thermal_clade_consensus_logos_representatives/by_domain/<domain>_clades.png
figures/thermal_clade_consensus_logos_representatives/by_domain/<domain>_hydropathy_volume_tracks.png
figures/thermal_clade_consensus_logos_representatives/thermal_clade_consensus_logos.html
figures/thermal_clade_consensus_logos_representatives/column_stats.tsv
figures/thermal_clade_consensus_logos_representatives/total_msa_conservation.tsv
figures/thermal_clade_consensus_logos_representatives/run_summary.tsv
figures/thermal_clade_consensus_logos_representatives/core_metric_summary.tsv
```

Older exploratory ASR, linker-path, and AF3-preparation scripts have been moved
to `scripts/msa_OGT/archive/`. Their old derived outputs have been moved to
`results/msa_OGT/archive/`.
