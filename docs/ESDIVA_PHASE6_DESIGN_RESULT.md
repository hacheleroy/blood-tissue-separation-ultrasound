# Phase-6 calibration and rigid-motion development result

25 October 2026. The phase-5 reserved comparison failed on seeds 2081 and
2093, which are permanently seen. [The phase-6 split](../configs/esdiva_phase6_split.json)
was published before source-resolved RF on design seeds 2309 and 2311. The
prospective development seed 2333 and fresh reserved seeds 2401 and 2411 have
not been opened. This document reports *development*, not another held-out
validation or a physical phantom result.

The [trial rules](../configs/esdiva_phase6_trial_rules.json) declared twelve
fixed configurations: the original phase-5 core input gate; eight local-phase
fractions and smoothing widths under wider tissue support; and three rigid
translation controls. The previous paired release limits remain: blood NRMSE
increase at most 0.020 for SVD and Axis 1, moving-tissue leakage ratio at most
0.900, exact roundtrip at most 1e-10, plus separately fixed slow-flow Axis 2
blood NRMSE at most 0.700 and leakage at most 0.180. Every complete trial
must pass all four scenes on *both* design seeds. The complete machine-readable
[decision](../runs/esdiva_phase6_design/decision.json) checks trial presence,
source hashes and every scene before declaring eligibility.

## Mechanistic findings

The original small core calibration fails its cross-angle coherence check on
seed 2309, while it passes on seed 2311's stationary, rigid and zero-flow
scenes. A vessel-guarded 32-row tissue support passes its independently
measured coherence check on both seeds: approximately 0.789 on seed 2309
stationary and rigid, but only 0.713–0.717 on seed 2311 stationary and rigid.
This is enough to score the fixed filters using a bounded, temporally checked
core ROI; it does not itself establish blood preservation or motion correction.
Source-resolved replay exactly reproduces the original central eight rows.
Known scatterer displacement remains evaluator-only.

On seed 2309 rigid fast, no local-phase fraction/width meets the moving-tissue
bound. The lowest SVD leakage ratio among the eight phase-only variants is
0.944 versus the required 0.900. At that setting, SVD blood NRMSE increases
0.033 and Axis 1 increases 0.025, both above 0.020. Gentler corrections
reduce blood damage but retain too much tissue. The cyclic eight-row shift
and its phase combinations increase leakage above 1.24 for both filters.
An additional matched full-width control avoids the narrow periodic edge yet
still raises leakage for SVD (0.353 versus 0.347) and Axis 1 (0.334 versus
0.298). The fitted magnitude translation tracks simulator displacement on
that scene (axial correlation 0.991 and RMSE 0.060 pixel; lateral correlation
0.997 and RMSE 0.021 pixel), measured **only after** the data-only fits. This
supports a diagnosis of a deficient complex tissue correction / filter
interaction rather than a simple inability to estimate gross translation;
the simulator trajectory is not independent physical ground truth.

On seed 2309 smooth slow, all eight phase-only variants satisfy moving-tissue
leakage and blood bounds; the fixed no-motion Axis 2 slow-flow control has
blood NRMSE 0.591 and tissue leakage 0.162. This favorable scene cannot
override the rigid-scene failure. The original phase-5 gate control is
blocked by input calibration, as it should be.

The seed 2311 rigid scene demonstrates the complementary failure: phase
fraction 0.50 at width 4 meets the blood and leakage limits (SVD increase
0.006 and leakage ratio 0.893; Axis 1 increase 0.004 and ratio 0.884), but
the same setting on seed 2309 has leakage ratios 0.991 and 0.985. The more
aggressive fraction 0.80 at width 2 reduces seed 2311 leakage to 0.839 and
0.825 while raising blood NRMSE by 0.094 and 0.081. Thus a single-scene
selection would hide the seed-dependent tradeoff.

![Paired blood-response increase versus tissue leakage for eight phase-only
configurations on two development seeds and moving scenes. The shaded region
passes both individual bounds.](figures/esdiva_phase6_design_gates.svg)

On the corrected seed 2311 slow-flow scene, the fixed no-motion Axis 2
control has blood NRMSE 0.544 and tissue leakage 0.156. Its zero-flow scene
is also source-resolved and scored. Some aggressive local-phase settings on
seed 2311 slow flow improve tissue leakage but exceed the blood increment
bound; the reviewer retains these failures rather than selecting by tissue
alone.

One seed 2311 smooth-motion acquisition mistakenly used flow multiplier 0.18,
whereas the declared `smooth_slow` scene uses 0.02. The scorer rejected this
as an unknown scene before producing a score. It was excluded, without tuning
or substituting its observations; the specified 0.02 scene was reacquired. Its
[exclusion audit](../runs/esdiva_phase6_design/excluded_acquisition.json)
retains the report and wider-replay hashes and records zero selection scores.

## Release and physical bridge

The all-scene [decision](../runs/esdiva_phase6_design/decision.json) is
**no design candidate passes**: all twelve declared trials have at least one
failure on the eight source-resolved scenes. The original core gate control
blocks all four seed 2309 scenes. Every wide-supported phase-only trial fails
at least one seed 2309 rigid leakage bound; narrow width 2 also incurs blood
damage on one or both rigid clouds. Rigid translation controls either worsen
seed 2309 leakage or damage blood on seed 2311. Exact roundtrip and the fixed
slow Axis 2 bounds do not rescue a failed moving-tissue/blood response.

The revised **design space and fresh reserve were frozen**, but no numerical
candidate qualified for a reserve freeze. Seed 2333 and reserved seeds 2401
and 2411 therefore remain unopened. Neither the failed phase-5 candidate nor
a phase-6 trial is eligible for phantom validation. The next model should
target complex speckle transport and boundary handling while measuring its
paired blood response, then be checked on additional declared development
clouds before freezing a single candidate for seed 2333. Do not retune from
the earlier seen 2081/2093 held-out results or infer physical validation from
the simulator's motion labels.

The [phase-5-to-phantom bridge](PHANTOM_PHASE5_BRIDGE.md) and its
[machine-readable contract](../configs/phantom_phase5_bridge.json) were
preregistered before physical capture. They fix the measured three-view
geometry, 32×16 wide tile and exact central eight-row core, first 16 of 64
stationary calibration frames, structural tissue/flow masks, independent
motion and flow references, 256-frame packet scoring, pump-on/off endpoints
and confidence-bound gates. The checker passes the conditional contract and
rejects an incomplete acquisition card. There are no physical acquisitions,
no scanner-specific card, and no confirmed simulation-qualified candidate.
Before capture, seal the actual geometry, masks, port source/configuration
hashes, optical/flow calibration and successful fresh-reserve decision in the
card. Source-resolved simulated blood NRMSE is unavailable from independent
physical pump-on/off captures; the bridge uses a separately specified pump
response endpoint.
