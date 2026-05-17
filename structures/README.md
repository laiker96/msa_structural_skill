# structures/

Local structure inputs and downloaded model caches.

Large generated/downloaded structures should stay out of git. This directory is
only a placeholder for repository-local examples. The maintained pipeline now
writes AlphaFold DB downloads under the home work directory by default:

```text
~/structural_evo_analysis/structures/afdb/
```

Optional user-provided query structure:

```text
~/structural_evo_analysis/structures/query.pdb
```

Step 06 downloads AFDB PDB models from accessions in `repset_metadata.tsv`:

```bash
./envs/structural_evo/bin/python scripts/structural_evo_analysis/06_download_afdb.py \
  --metadata results/photoHymenobact_example/repset_metadata.tsv \
  --dest ~/structural_evo_analysis/structures/afdb
```

Step 07 scores structures from this directory. Use `--skip-aggrescan3d` if the
optional Aggrescan3D environment is not installed.
