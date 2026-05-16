# data/

Committed metadata used by the structural-evolution pipeline.

```text
growth_temp_dataset_OGTFinder.tsv
```

This table is taxonomy-keyed. Step 01 reads `TaxID=` from UniRef-style hit
headers and joins optimum growth temperature metadata into
`repset_metadata.tsv`. The clade annotator can then use:

```text
--metadata results/<run>/repset_metadata.tsv --trait-column ogt
```

For OGT-based clade calls, use explicit thresholds such as:

```text
--low-threshold 20 --high-threshold 45
```

Do not use this file as a sequence database.
