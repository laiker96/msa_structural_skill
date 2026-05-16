# scripts/

The maintained workflow for this repository is:

```text
scripts/structural_evo_analysis/
```

It provides a reusable single-query protein pipeline:

```text
query FASTA
  -> MMseqs2 homolog search
  -> MAFFT MSA
  -> IQ-TREE phylogeny
  -> metadata-driven clade annotation
  -> conserved-position tables
  -> optional AFDB download and structure scoring
```

Use `./envs/structural_evo/bin/python` for these scripts.

## Directory Map

```text
scripts/
├── structural_evo_analysis/   maintained generalized pipeline
└── msa_OGT/                   legacy ANKros-specific reference code
```

`scripts/msa_OGT/15_compute_solubility_aggregability.py` and its scorer modules
are still reused by `structural_evo_analysis/07_score_structures.py`. Other
`msa_OGT` steps are retained as reference material and should not drive new
skill behavior unless they are deliberately generalized.

## Long Runs

MMseqs searches, IQ-TREE, AFDB downloads, and structure scoring can be
long-running. Use `tmux` and log commands under `logs/`.
