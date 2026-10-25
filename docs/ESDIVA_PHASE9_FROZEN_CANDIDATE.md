# Phase 9: tissue support repair, frozen development comparison

25 October 2026. The [split](../configs/esdiva_phase9_split.json),
[candidate](../configs/esdiva_phase9_candidate.json), acquisition, scorer, and
reviewer are sealed **before RF generation for seeds 3001 and 3019**. Seeds
2309, 2311, 2503, 2521, 2707 and 2729 are development data already seen.
Seed 2729's previously generated RF checkpoints were [replayed](../examples/esdiva_phase9_replay_seen.py)
for the three scenes that phase 8 rejected before writing IQ. The completed
zero-flow IQ replay matches the old saved arrays exactly. The [24-scene
retrospective screen](../runs/esdiva_phase9_design/retrospective.json) has zero
failures under the revised rules, but cannot establish new-cloud performance.

## Frozen source and input eligibility

The known Example 22 tilted-vessel geometry fixes a 5 mm guard on the tissue
mask, on both the 32 × 16 support tile (192 pixels) and the eight-row scoring
tile (48 pixels). This is chosen by geometry without evaluating isolated
blood. The acquisition **archives** source-resolved RF and IQ even if the
evaluator-only blood-to-tissue power check exceeds 0.100; the scorer records
an input failure and performs no filtering in that case. Both wide and core
regions must satisfy the check. The earlier lateral-edge core mask was
contaminated on seed 2729; its replacement is part of the new acquisition
contract, not a retrospective reinterpretation of phase-8 eligibility.

For the stationary, noisy tissue-only prior calibration, the first eight
frames fit per-view spatial phase relative to the central angle. The
cross-view products are masked to the fixed tissue region, smoothed at 0.4
pixel with zero extension and normalized by the similarly smoothed mask.
The phase map is applied to calibration frames **only to measure eligibility**:
the last eight frames must retain wide coherence ≥0.700, core coherence
≥0.500 and core validation coherence ≥0.85 times core train coherence.
The corresponding train minima and bounded noise power ≤4 also apply.
The regular filters receive the original IQ and calibration, with no mapped
IQ substituted. No blood truth, prescribed motion, or evaluation IQ enters
the tissue mask, phase map, mode selection, or filter fit. The two calibration
halves share a stationary scatterer anatomy, so this gate alone cannot prove
robustness to new anatomy.

## Motion and paired rule

The [mixture-only phase-consistency selector](../examples/esdiva_phase8_selector.py)
chooses between reference transport and a hybrid on the fixed 5 mm guard,
using only frames 0–63. The reference parameters stay as in phase 8. The
hybrid combines accepted local phase at fraction **0.6**, spatial sigma 4,
and half-strength nonperiodic unitary IQ transport. Frames 64–79 are scored.
On both the mixture and a separately refitted tissue-plus-noise background,
SVD cut 1 and Axis 1: Motion-compensated virtual long ensemble cut 1 must
have paired blood NRMSE increase ≤0.020 versus their no-motion controls on
**every scene**. Their moving-scene tissue leakage ratio must be ≤0.900 in
rigid fast and smooth slow. Paired roundtrip is ≤1e-10. The uncorrected slow
Axis 2: Motion-demodulated multi-lag separation rank 4/cap 0.4 must have
blood NRMSE ≤0.700 and leakage ≤0.180. The scorer also retains the SVD,
Baranger, Stephen Lee AER and Axis 3: Local ultrasound STAP comparators.

Both 3001 and 3019 must satisfy **all four scenes**, with no missing scene
or source failure. No threshold, selector, motion parameter, acquisition
rule or numerical source may be changed after new RF without declaring a
different candidate and fresh seeds. A complete pass permits a separate
prospective freeze for seed 2333. Seeds 2333, 2401 and 2411 are unopened;
no physical phantom evidence is claimed here.
