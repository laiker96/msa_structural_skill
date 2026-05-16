# structures/

Local structure inputs and downloaded model caches.

Large generated/downloaded structures should stay out of git. The generalized
pipeline writes AlphaFold DB downloads here by default:

```text
structures/structural_evo_analysis/afdb/
```

Optional user-provided query structure:

```text
structures/structural_evo_analysis/query.pdb
```

Step 06 downloads AFDB PDB models from accessions in `repset_metadata.tsv`:

```bash
./envs/structural_evo/bin/python scripts/structural_evo_analysis/06_download_afdb.py \
  --metadata results/photoHymenobact_example/repset_metadata.tsv \
  --dest structures/structural_evo_analysis/afdb
```

Step 07 scores structures from this directory. Use `--skip-aggrescan3d` if the
optional Aggrescan3D environment is not installed.
