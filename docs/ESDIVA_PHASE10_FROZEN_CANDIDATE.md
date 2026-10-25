# Phase 10: core support and reference-mode selection freeze

25 October 2026. This [split](../configs/esdiva_phase10_split.json),
[candidate](../configs/esdiva_phase10_candidate.json), acquisition, scorer,
and all-scene reviewer are fixed **before generating RF for seeds 3307 and
3323**. Seeds 2309, 2311, 2503, 2521, 2707, 2729, 3001, and 3019 are
already seen development clouds. The [32-scene retrospective screen](../runs/esdiva_phase10_design/retrospective.json)
passes the revised source, calibration, paired blood and tissue gates, but
was used to choose these parameters and is not an independent confirmation.

## Anatomical tissue support

The prior-known tilted vessel geometry fixes a **5.5 mm core guard** (40
tissue pixels in the eight-row score tile) and a **5 mm wide guard** (192
pixels in the 32-row motion tile). The extra half millimetre in the core
reduces the now-seen 3001 fast-flow blood/tissue source power from 0.1139
and 0.1126 to 0.080 and 0.079. These are *evaluator-only* measurements;
neither isolated source entered the mask rule. The scorer requires source
power ≤0.100 in both regions and archives an ineligible scene without
filtering. The acquisition report also records the older lateral-edge
annotation as a diagnostic.

The prior stationary tissue calibration uses the phase-9 masked angular
phase map only for eligibility. Its first eight frames fit the map, its last
eight validate wide coherence ≥0.700 and core coherence ≥0.500, core
validation ≥0.85 of core train, and bounded noise gain ≤4. All filters
receive **unmodified mixture IQ** and the original calibration, and the
phase-9 numerical method families remain unchanged. A stationary calibration
on the same anatomy across both halves is a temporal check, not a test on
new anatomy.

## Reference transport decision and paired gate

The phase-8 pre-score mixture-only phase-loss rule first requests reference
or hybrid. Reference is used only if the prior calibration reference has
adequate amplitude at **≥90% of the 5.5 mm core tissue pixels in every
view**, relative to 20% of that view's median wide-tissue reference
amplitude. Otherwise the hybrid applies accepted local phase fraction
**0.75**, spatial sigma 4 and half-strength nonperiodic unitary IQ shift.
The reference option retains transport strength 1.1 and its existing
boundary handling. The new support check uses calibration alone; selection
uses only calibration and pre-score mixture frames 0–63. It retains
reference for seen 2309 rigid fast and rejects reference for seen 3001,
where reference had failed the smooth-slow tissue gate.

The frozen paired rule remains: SVD cut 1 and Axis 1: Motion-compensated
virtual long ensemble cut 1 must increase paired blood NRMSE by ≤0.020
against no motion in **every** scene. Their tissue leakage ratio must be
≤0.900 in rigid fast and smooth slow. Roundtrip NRMSE is ≤1e−10. The
uncorrected slow-flow Axis 2: Motion-demodulated multi-lag separation
rank 4/cap 0.4 requires blood NRMSE ≤0.700 and tissue leakage ≤0.180.
The baseline records SVD, Baranger SSM, Stephen Lee AER and Axis 3:
Local ultrasound STAP as comparators.

Both new clouds, **3307 and 3323**, must pass all four scenes. A missing
or failed scene disqualifies the candidate; there is no threshold or
strength adjustment after these RF data are opened. A complete new-bank
pass would permit a **separate prospective freeze** for seed 2333. Seed
2333 and reserves 2401/2411 remain unopened, and no physical phantom
claim follows from this simulation freeze.
