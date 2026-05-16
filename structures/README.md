# structures/

External structural data pools used by the active pipelines.

Large structure files are gitignored. The active FAD/FMN holo builder does not
read crystals or derived holo models from this directory; it stages its own
crystal inputs and writes its own outputs under `results/structural/`.

## Layout

```text
structures/
├── README.md
├── af3/        AlphaFold 3 server predictions, extracted from submitted jobs
└── afdb/       AlphaFold DB models for MSA representatives
```

Retired legacy crystal, ligand-template, and derived-holo folders were
intentionally removed. Current crystal provenance for the holo builder lives in
`results/structural/inputs/crystals/crystal_inputs.tsv`, and current holo model
outputs live in `results/structural/docked_holo/`.

## Reproducing Current Contents

AFDB representative structures are downloaded by the MSA pipeline utility:

```bash
python scripts/msa_OGT/13_download_afdb.py
```

AF3 server predictions require a manual server submission. After downloading a
fresh `folds_*.zip` archive, extract it here:

```bash
rm -rf structures/af3
unzip -q structures/folds_*.zip -d structures/af3
```

The active structural holo builder stages its required RCSB crystals into the
results tree:

```bash
bash scripts/structural/run_pipeline.sh
```
