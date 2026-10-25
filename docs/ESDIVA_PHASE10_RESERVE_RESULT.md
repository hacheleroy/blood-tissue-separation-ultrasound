# Phase 10: reserved scatterer clouds 2401 and 2411

25 October 2026. The [reserve protocol](ESDIVA_PHASE10_RESERVE_FREEZE.md)
was committed as `67fd0e6` before either cloud's RF was generated. The
[fixed eight-scene decision](../runs/esdiva_phase10_reserve/decision.json) is
**reserved bank passed; physical bridge eligible for acquisition preparation**.
There are no missing scenes or gate failures. All [source reports, wide-IQ
replay audits and scores](../runs/esdiva_phase10_reserve/fresh/) are bound by
the [SHA-256 manifest](../runs/esdiva_phase10_reserve/manifest.json). The RF
checkpoints remain local and are not part of the Git archive.

All eight scenes passed the 40-pixel core and 192-pixel wide source-purity
gates. Maximum evaluator-only blood/tissue power was **0.00603** in the core
and **0.00551** in the wide support, each below 0.100. Every prior-calibration
gate passed. Minimum wide and core validation coherence were **0.89193** and
**0.90350**, above their 0.700 and 0.500 limits, respectively; the core
train-relative requirement also passed. RF source additivity error was at
most **1.46 × 10⁻⁷**; the scored-core IQ replay errors were zero. Each scene
selected **hybrid** correction by the frozen mixture/prior-calibration rule.

| Cloud and scene | SVD blood Δ / tissue ratio | Axis 1 blood Δ / tissue ratio |
|:--|:--|:--|
| 2401 stationary fast | +0.0001 / 0.9993 | +0.0001 / 0.9996 |
| 2401 rigid fast | −0.0002 / 0.7869 | −0.0047 / 0.8188 |
| 2401 smooth slow | −0.0446 / 0.6471 | −0.0392 / 0.6925 |
| 2401 smooth zero | −0.0379 / 0.7426 | −0.0330 / 0.7697 |
| 2411 stationary fast | +0.0041 / 0.7071 | +0.0034 / 0.7753 |
| 2411 rigid fast | +0.0013 / 0.5477 | −0.0049 / 0.5787 |
| 2411 smooth slow | +0.0055 / 0.6165 | −0.0369 / 0.6432 |
| 2411 smooth zero | −0.0058 / 0.6234 | −0.0211 / 0.6471 |

Blood Δ is the selected-minus-uncorrected paired added-blood IQ NRMSE;
tissue ratio is the selected/uncorrected refitted-background leakage. SVD
cut 1 and Axis 1: Motion-compensated virtual long ensemble cut 1 passed the
**≤0.020** blood-increase gate in every scene; the maximum was **+0.00547**.
Their maximum *required moving-scene* tissue ratio was **0.81876**, below
0.900. Stationary and zero-flow ratios are reported without that moving
suppression requirement. Both variant roundtrips passed 1e−10.

The uncorrected smooth-slow Axis 2: Motion-demodulated multi-lag separation
rank 4/cap 0.4 control had blood NRMSE/tissue leakage **0.56081/0.15438**
on 2401 and **0.55842/0.16449** on 2411, passing the 0.700/0.180 limits.
The archived scores also include Baranger SSM, Stephen Lee AER and Axis 3:
Local ultrasound STAP; this paired gate does not establish superiority over
those comparators.

The reviewer reproduced the decision byte for byte on a second run, all 29
manifest entries verified, and the 206-test local numerical suite passed.
The [paired gate figure](figures/esdiva_phase10_reserve_gates.svg) shows both
moving scenes and the frozen bounds.

## Interpretation and physical handoff

These are two independently declared scatterer clouds with four correlated
scenes per cloud and 16 scored frames per scene. Together with the earlier
2333 prospective result, they support the frozen hybrid path on this **RF
simulator**. They do not validate a physical phantom or an independent
reference-mode selection: none of these eight scenes selected reference.

The [phase-10 phantom bridge](PHANTOM_PHASE10_BRIDGE.md) now binds this
decision and archive. Its geometry-aware port has been checked on already
seen simulated IQ, including a reference selection, but it has no scanner
measurement. The [incomplete acquisition card and 180-slot capture
inventory](../runs/phantom_phase10_plan/) specify the matched development
and fresh-remount physical captures. Their scanner, structural geometry,
beamformer, optical, stage, flowmeter/PIV and raw RF fields must be measured
and sealed before the first physical capture; no physical validation result
exists yet.
