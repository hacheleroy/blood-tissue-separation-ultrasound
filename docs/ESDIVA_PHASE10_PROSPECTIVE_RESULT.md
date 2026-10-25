# Phase 10: prospective seed 2333 result

25 October 2026. The [prospective contract](ESDIVA_PHASE10_PROSPECTIVE_FREEZE.md)
was committed as `86ce030` **before generating RF for seed 2333**. The
[fixed reviewer decision](../runs/esdiva_phase10_prospective/decision.json) is
**prospective development passed; reserves unopened**. All four required
scenes completed, with no missing scenes or failures. Their [source reports,
RF replay audits, and scores](../runs/esdiva_phase10_prospective/fresh/2333/)
are bound by a [SHA-256 manifest](../runs/esdiva_phase10_prospective/manifest.json).
Large RF checkpoints remain local and are not in Git.

The 40-pixel core and 192-pixel wide anatomical tissue supports passed the
evaluator-only source-purity check: maximum blood/tissue power was **0.05623**
in the core and **0.01517** in the wide support (each limit 0.100). Every
prior-calibration gate passed. Minimum wide and core validation coherence
were **0.93410** (limit 0.700) and **0.74948** (limit 0.500, with the
train-ratio condition). RF source additivity error was at most
**1.45×10⁻⁷**, and scored-core IQ replay error was zero. Each scene selected
hybrid motion correction using the predeclared rule.

| Scene | SVD blood Δ / tissue ratio | Axis 1 blood Δ / tissue ratio |
|:--|:--|:--|
| Stationary fast | +0.0003 / 0.9984 | +0.0020 / 0.9977 |
| Rigid fast | +0.0123 / 0.8212 | −0.0614 / 0.8581 |
| Smooth slow | −0.0288 / 0.7332 | −0.0337 / 0.8009 |
| Smooth zero | −0.0318 / 0.7367 | −0.0347 / 0.8002 |

Blood Δ is the selected-minus-uncorrected paired added-blood IQ NRMSE; a
positive value means worse blood preservation. Tissue ratio compares the
selected and uncorrected background leakage, so a smaller value means more
suppression. SVD cut 1 and Axis 1: Motion-compensated virtual long ensemble
cut 1 passed their **≤0.020 blood Δ in every scene** gate. The largest
increase was **+0.01234**. Both methods passed the **≤0.900** tissue ratio
gate in rigid fast and smooth slow; the largest required ratio was **0.85814**.
The stationary and zero-flow ratios are reported but have no moving-tissue
reduction requirement. Both variant roundtrips passed 1e−10.

The uncorrected smooth-slow Axis 2: Motion-demodulated multi-lag separation
rank 4/cap 0.4 control had blood NRMSE **0.63068** and tissue leakage
**0.16330**, passing its **0.700/0.180** limits. The archived shared method
rows also retain Baranger SSM, Stephen Lee AER, and Axis 3: Local ultrasound
STAP, alongside the SVD and three-axis comparators; this paired motion gate
does not imply superiority over those comparators.

## Interpretation and next step

This is **one prospective scatterer cloud** with four correlated scenes and
16 scored frames per scene. It confirms the frozen phase-10 candidate on
software RF simulation after the two-cloud development pass. It does not
validate a physical phantom, tissue heterogeneity beyond this generator, or
the reference-mode path on a fresh cloud: all four 2333 scenes chose hybrid.

Seeds **2401/2411 remain unopened**. The next comparison needs its own
predeclared, hash-bound reserve protocol, including a gate that can exercise
reference-mode support on independent tissue if that claim is intended.
Physical phantom acquisition then needs independent displacement, blood
preservation and tissue-leakage measurements with stationary, rigid,
nonrigid, slow-flow and zero-flow controls.
