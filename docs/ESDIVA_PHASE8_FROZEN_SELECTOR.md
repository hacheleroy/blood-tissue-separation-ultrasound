# Phase-8 guarded calibration and measured phase selector

25 October 2026. This is one **development candidate** chosen using the
already seen 2309, 2311, 2503 and 2521 clouds. It does not rescue the
historically frozen phase-7 result. The [new split](../configs/esdiva_phase8_split.json)
declares 2707 and 2729 as fresh development confirmation clouds before their
RF is generated. The [candidate configuration](../configs/esdiva_phase8_candidate.json)
hashes the [selector](../examples/esdiva_phase8_selector.py),
[score](../examples/esdiva_phase8_score.py) and
[review](../examples/esdiva_phase8_review.py). All source-resolved evaluations
use the original phase-5 numerical source and the same four scenes, source
replay checks and paired limits.

The selected wide-tile tissue guard is 4 or 5 mm, chosen by maximum minimum
angular coherence on the first eight prior calibration frames. The remaining
eight frames must independently reach 0.700 coherence, bounded noise power
gain at most 4, and at least 128 tissue pixels. This repairs seed 2503's
phase-7 near-threshold guard-5 failure without lowering its coherence gate.
Only the **selected** guard is validated; a guard cannot be switched after its
validation result is known. The core's phase-7 wide-supported diagnostic rule
is retained, including its separate validation floor and train/validation
ratio.

Guard 4 uses accepted local tissue phase at fraction 0.5 and smoothing width
4 pixels. Guard 5 uses the same local phase plus a half-strength rigid shift
in an exactly reversible, **nonperiodic** skew-derivative operator on the
original core. One narrow exception selects the earlier in-bounds complex
reference phase transport at strength 1.1: weighted phase discordance between
the locally corrected mixture and its prior calibration must exceed 0.28,
and the reference mode's additional discordance must be below 0.18. The
discordances use only tissue check pixels and frames 0–63, before the scored
packet 64–79. The mode is selected from mixture data once and applied to
the paired evaluator's isolated background. Source truth, prescribed motion
and scene labels cannot influence the selector.

The rule is intentionally a **hypothesis**, not a demonstrated physical
criterion. On the seen 2309 rigid scene it selects reference despite a worse
raw phase loss, because that mode passes the paired blood and leakage gates;
the simple rule that minimizes phase loss or pre-score filtered tissue power
selected a failing mode. The cutoffs were chosen after inspecting those
seen outcomes. The new clouds can therefore refute the rule and cannot be
used to retune it while still being called confirmation.

The [partial all-scene review](../runs/esdiva_phase8_design/partial_decision.json)
has zero failures on all sixteen seen scenes and lists the eight new scenes
as missing. The moving-scene numbers below are **retrospective**. Each pair
is SVD / Axis 1, with blood NRMSE change and background tissue leakage ratio
relative to its own uncorrected control.

| Seen cloud / scene | Selected | Blood NRMSE Δ, SVD / Axis 1 | Tissue ratio, SVD / Axis 1 |
|:--|:--|--:|--:|
| 2309 rigid fast | Reference | +0.0180 / +0.0149 | 0.8494 / 0.8979 |
| 2309 smooth slow | Hybrid | −0.0424 / −0.0316 | 0.8350 / 0.8592 |
| 2311 rigid fast | Hybrid | +0.0136 / +0.0090 | 0.6808 / 0.7195 |
| 2311 smooth slow | Hybrid | −0.0496 / −0.0447 | 0.5739 / 0.6265 |
| 2503 rigid fast | Local | −0.0009 / −0.0017 | 0.8749 / 0.8824 |
| 2503 smooth slow | Local | −0.0105 / −0.0132 | 0.8119 / 0.8350 |
| 2521 rigid fast | Hybrid | −0.0008 / −0.0095 | 0.8222 / 0.8746 |
| 2521 smooth slow | Hybrid | −0.0916 / −0.0674 | 0.8132 / 0.8734 |

The all-scene reviewer requires source and input eligibility, paired
roundtrip error at most 1e-10, blood NRMSE increase at most 0.020 versus
the same-scene no-motion control for SVD cutoff 1 and **Axis 1:
Motion-compensated virtual long ensemble** cutoff 1 in all scenes, and
tissue leakage ratio at most 0.900 for both on rigid-fast and smooth-slow.
The no-motion slow-flow **Axis 2: Motion-demodulated multi-lag separation**
control retains its blood 0.700 and tissue 0.180 limits. Baseline records
also retain Baranger SSM, Stephen Lee AER and **Axis 3: Local ultrasound
STAP**. A missing or failed scene blocks the path to prospective 2333;
2401/2411 remain a separate unopened reserve. Physical phantom acquisition
still requires a passing reserve and a scanner-specific geometry card.
