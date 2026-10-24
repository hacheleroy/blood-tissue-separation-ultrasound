# eSDIva broad-support local phase registration

24 October 2026. The [pre-new-cloud manifest](../configs/esdiva_phase5_design_freeze.json)
fixes the [candidate](../configs/esdiva_phase5_candidate.json), source code,
decision rule and all twelve prior development records before new scatterer
cloud 2111 is acquired. Reserved seeds 2081 and 2093 remain unopened.

The pinned Example 22 RF checkpoints reconstruct a 32-row tissue tile at
14.4–17.5 mm depth. The central eight rows exactly replay the published IQ
score tile, sample for sample. A 5 mm guard around the known tilted vessel
selects 192 tissue pixels; 127 fit the phase field, and 65 disjoint pixels
check each frame and view. Only the measured mixture, stationary disjoint-time
tissue calibration and geometry reach the filter. Isolated blood and tissue
and known scatterer motion are evaluator-only. The fit smooths a normalized
complex tissue-reference cross product at four-pixel scale, unwraps phase
along depth and applies half the phase if its disjoint-pixel improvement is
above 0.01. This is local phase registration, **not** a measured geometric
displacement or spatial resampling. Its unit-magnitude factor keeps the
frozen complex-blood roundtrip exact. Separate mixture and tissue-plus-noise
fits implement the paired source check.

The original eight-row input gate and a fixed wider-ROI gate both require
calibration coherence at least 0.7 and bounded input noise-power gain at most
4. The wider ROI is checked for blood contamination by the evaluator (at most
0.1). Both gates pass the three prior clouds in every scene; wider ROI
contamination spans 0.0037–0.0203. The full 24-file phase-2 frozen history
and numerical source hash remain intact.

Values below are **local-phase minus no-motion paired added-blood IQ NRMSE /
local-phase divided by no-motion refitted tissue leakage**; the blood change
must be at most +0.020 in all scenes, and the leakage ratio at most 0.900 in
both moving scenes. Both SVD discard 1 and Axis 1: Motion-compensated virtual
long ensemble must pass. All three previous clouds meet these development
conditions:

| Seed | Rigid SVD | Rigid Axis 1 | Smooth slow SVD | Smooth slow Axis 1 |
|---:|---:|---:|---:|---:|
| 2083 | −0.002 / 0.757 | −0.001 / 0.748 | −0.017 / 0.664 | −0.016 / 0.673 |
| 2097 | +0.009 / 0.826 | +0.005 / 0.838 | +0.002 / 0.863 | −0.004 / 0.867 |
| 2101 | +0.004 / 0.845 | −0.003 / 0.856 | −0.079 / 0.763 | −0.056 / 0.792 |

Stationary and moving zero-flow controls also meet the blood-increase bound;
paired roundtrip NRMSE is numerically zero in every scene. The separate slow
Axis 2: Motion-demodulated multi-lag separation rank-4, gain-cap-0.40
setting meets its existing ≤0.700 blood NRMSE and ≤0.180 leakage bounds on
the three previous clouds. Global SVD, Baranger SSM, Stephen Lee AER and
Axis 3: Local ultrasound STAP retain their phase-4 control settings. The
local-phase experiment does not claim to solve the slow or zero-flow
identifiability problem.

The predeclared [decision program](../examples/esdiva_phase5_decide.py) requires
all four scenes on new cloud 2111 to pass the same gates. Any failure keeps
2081/2093 unopened. A pass permits the complete score and source freeze to
`main` before the reserve is used. These are software RF simulation results,
not physical phantom or in vivo validation.
