# test/

Small committed example inputs for smoke tests and documentation.

```text
photoHymenobact.fa
query.pdb   optional real query structure, when available
```

`photoHymenobact.fa` is the bundled single-sequence protein FASTA. A matching
`query.pdb` can be placed here for local smoke tests, but no synthetic or
placeholder PDB is committed because vulnerability ranking requires a real
query structure.

Runtime outputs and downloaded AlphaFold DB structures should not be written
here. They default to:

```text
~/structural_evo_analysis/results/
~/structural_evo_analysis/structures/
```
