# Phase 10: reserved-bank freeze

25 October 2026. The [prospective 2333 decision](../runs/esdiva_phase10_prospective/decision.json)
passed all four scenes after its own prior freeze. The [split](../configs/esdiva_phase10_reserve_split.json),
[candidate](../configs/esdiva_phase10_reserve_candidate.json),
[scorer](../examples/esdiva_phase10_reserve_score.py) and
[reviewer](../examples/esdiva_phase10_reserve_review.py) are sealed **before
generating RF for reserved scatterer clouds 2401 and 2411**. The split and
candidate bind the prospective decision and archive manifest by SHA-256.
The reviewer additionally checks that the numerical settings, scene matrix,
acquisition and thresholds equal the prospective configuration. The reserve
scorer changes only the bound parent result and split, admitted seeds, and
diagnostic status label. It retains the acquisition, source gate, 5 mm wide
and 5.5 mm core support, prior calibration eligibility, reference safeguard,
hybrid correction, method settings and score window unchanged.

Each cloud has stationary fast, rigid fast, smooth slow and smooth zero scenes
with 16 prior calibration frames and 80 evaluation frames. Score frames
64–79. RF source and noise identities, broad IQ replay, the 40-pixel core
and 192-pixel wide source-purity limit (blood/tissue power ≤0.100), and
prior calibration gates are required for every scene. No failed or missing
scene can be discarded. Both selected and no-motion roundtrips must be
≤1e−10. SVD cut 1 and Axis 1: Motion-compensated virtual long ensemble cut 1
must each have paired blood NRMSE increase ≤0.020 in every scene, and tissue
leakage ratio ≤0.900 in both moving scenes. The uncorrected smooth-slow Axis
2: Motion-demodulated multi-lag separation rank 4/cap 0.4 control must have
blood NRMSE ≤0.700 and tissue leakage ≤0.180. Report stationary and zero-flow
leakage without imposing a moving-tissue suppression threshold on them.

All eight scenes and all source and paired gates must pass for a reserve-bank
pass. The reference-mode safeguard remains a fixed selection rule; an absence
of reference-mode choices is reported as limited coverage, not retuned after
the fact. Source truth and known motion stay evaluator-only; the filtering
and mode choice use regular IQ and prior calibration. Archived comparators
include SVD, Baranger SSM, Stephen Lee AER, Axis 2 and Axis 3: Local
ultrasound STAP, but the paired motion gate is specifically SVD and Axis 1.

A pass permits **preparation** for independent physical phantom acquisition,
not an experimental validation claim. Physical data require a measured
acquisition card and new preacquisition geometry bridge. No experimental RF
has been supplied. A reserve failure remains a failure and blocks that
candidate from a phantom performance claim; the failed clouds become seen.
