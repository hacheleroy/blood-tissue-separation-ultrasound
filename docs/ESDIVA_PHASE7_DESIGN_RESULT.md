# Phase-7 boundary-aware complex transport: development result

25 October 2026. The [phase-7 split](../configs/esdiva_phase7_split.json)
declared seeds 2503 and 2521 before their RF was generated. The
[single diagnostic source and config](ESDIVA_PHASE7_FROZEN_DIAGNOSTIC.md) were
committed on `main` before opening either cloud. This is a development study,
not a fresh held-out or physical phantom validation. The fixed
[all-scene decision](../runs/esdiva_phase7_design/decision.json) is **no
qualified candidate; prospective and reserve unopened**. All sixteen
old/new scene scores are present. The eight new acquisition reports, wider
source-replay audits and scores are archived under
[`runs/esdiva_phase7_design`](../runs/esdiva_phase7_design/).

## Source and input eligibility

The pinned Example 22 source resolved tissue and blood RF separately for four
scenes on each new seed. Source-additivity errors were 1.41–1.44×10⁻⁷, the
wider replay verified sources and noise, and its central eight IQ rows matched
the original scoring tile exactly in all eight scenes. The regular filter
received mixture IQ, prior tissue calibration and structural geometry; source
components and prescribed motion were evaluation-only.

| New seed and scene | Wide coherence train | Wide coherence validation | Frozen input gate ≥0.700 |
|:--|--:|--:|:--|
| 2503 stationary fast | 0.700056 | 0.699431 | Fail, unscored |
| 2503 rigid fast | 0.700020 | 0.699395 | Fail, unscored |
| 2503 smooth slow | 0.696424 | 0.700591 | Fail, unscored |
| 2503 smooth zero | 0.696483 | 0.701045 | Fail, unscored |
| 2521 stationary fast | 0.7927 | 0.7890 | Pass |
| 2521 rigid fast | 0.7927 | 0.7889 | Pass |
| 2521 smooth slow | 0.7875 | 0.7929 | Pass |
| 2521 smooth zero | 0.7945 | 0.7969 | Pass |

The near-threshold seed 2503 is a calibration-robustness finding, not a
paired filter failure or success. Its four scores deliberately contain no
filter variants. Lowering the 0.700 gate after seeing this seed would break
the fixed decision; no such score is used for selection.

## Paired blood and tissue response

The candidate interpolates only the *prior complex calibration reference*
within its measured support, masks invalid/weak regions, and applies a
phase-only correction to the exact core. This avoids periodic wrap and has an
exact inverse. The table shows change relative to the same scene's
no-correction control. Both SVD cutoff 1 and **Axis 1:
Motion-compensated virtual long ensemble** cutoff 1 require blood NRMSE
increase ≤0.020 on every scene and tissue leakage ratio ≤0.900 in rigid fast
and smooth slow.

| Seed / moving scene | SVD blood Δ | SVD tissue ratio | Axis 1 blood Δ | Axis 1 tissue ratio |
|:--|--:|--:|--:|--:|
| 2309 rigid fast | +0.0180 | 0.8494 | +0.0149 | 0.8979 |
| 2309 smooth slow | −0.0324 | 1.0480 | −0.0257 | 1.0539 |
| 2311 rigid fast | +0.0210 | 0.8829 | +0.0186 | 0.8844 |
| 2311 smooth slow | +0.1117 | 0.8401 | +0.0600 | 0.8509 |
| 2521 rigid fast | +0.1760 | 2.3576 | +0.1669 | 2.4304 |
| 2521 smooth slow | −0.0813 | 2.5838 | −0.0491 | 2.5244 |

![Paired blood error and tissue leakage on eligible moving development scenes.
Dashed lines are the frozen limits; seed 2503 is absent because its input gate
failed.](figures/esdiva_phase7_development_gates.svg)

Stationary fast on seed 2521 was unchanged by the correction. Smooth zero
reduced paired blood NRMSE by 0.0847 (SVD) and 0.0516 (Axis 1), while its
tissue ratios increased to 2.3616 and 2.1932. The fixed no-correction
**Axis 2: Motion-demodulated multi-lag separation** slow-flow control on
2521 remained within its separate bounds: blood NRMSE 0.6261 ≤0.700 and
tissue leakage 0.1479 ≤0.180. The no-correction records also retain
Baranger SSM, Stephen Lee AER and **Axis 3: Local ultrasound STAP**. Maximum
paired roundtrip error over scored scenes and both variants was 5.22×10⁻¹⁶,
well below 1×10⁻¹⁰. Algebraic reversibility does not establish useful
separation.

The 2521 rigid case contradicts generalization from the favorable 2309 rigid
response, despite similar gross fitted motion and measured in-bounds
reference support. On 2521 rigid, only about 51% of core pixels met the
weak-reference support test on average, and mixture/background fits accepted
94/95 view-frames. This is a diagnostic association, not proof of a single
failure mechanism. The tissue leakage more than doubles even as the
roundtrip remains exact. The prescribed simulator displacement is not used
as independent physical motion truth.

## Decision and next experiment

No passing candidate was frozen for prospective seed **2333**. That seed and
fresh reserve **2401/2411** remain unopened. Seed 2503 is now seen and must
not be relabeled held out or rescued by a post-hoc gate. The conditional
[phase-5-to-phantom geometry bridge](PHANTOM_PHASE5_BRIDGE.md) remains a
preacquisition contract; no physical data or scanner-specific card exists.

The next model should address both calibration support across clouds and the
seed-dependent blood/tissue response. In particular, a measured spatial
phase-consistency or nonrigid residual criterion could decide when to use a
local phase correction versus reference transport, but gross translation
magnitude alone does not distinguish the old rigid and smooth cases. Develop
that rule only on now-seen development data, declare another fresh
confirmation bank and freeze its code, input gate and paired decision before
opening that bank. A full pass there is needed before a separate prospective
freeze for 2333; a later 2333 pass is needed before any fresh reserve or
physical phantom claim.
