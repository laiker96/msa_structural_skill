# sequences/

Bundled example query sequences.

```text
photoHymenobact.fa
```

Use this FASTA as the reproducible test query for the generalized structural
evolution pipeline:

```bash
bash scripts/structural_evo_analysis/run_pipeline.sh \
  sequences/photoHymenobact.fa \
  results/photoHymenobact_example/repset_metadata.tsv \
  ogt
```

The pipeline expects exactly one protein sequence in the query FASTA.
