# msa_OGT

Legacy ANKros-specific pipeline retained as reference material.

The reusable skill workflow is `scripts/structural_evo_analysis/`, not this
directory. Do not add new generalized behavior here unless the goal is
explicitly to preserve or migrate an ANKros-specific method.

Currently reused by the main pipeline:

```text
15_compute_solubility_aggregability.py
structure_score_camsol.py
structure_score_aggrescan3d.py
```

These modules provide the local CamSol-style scorer, Aggrescan3D wrapper, and
score-merging logic used by:

```text
scripts/structural_evo_analysis/07_score_structures.py
```

Historical scripts in this directory include ANKros domain boundaries,
class-I photolyase filters, and OGT/regime assumptions. Treat those as
project-specific examples, not defaults for the generalized skill.
