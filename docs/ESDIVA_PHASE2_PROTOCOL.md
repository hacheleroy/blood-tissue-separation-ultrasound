# eSDIva longer-ensemble motion and flow protocol

24 October 2026. The source-resolved RF adapter uses pinned eSDIva Example 22
at `ab7eaef6d23428762ab5b9079ab03a5aa99af085`. Five matched scenes are
run on development seed 2026. Reserved eSDIva scatterer seeds 2031 and 2053
must remain unopened until the source, reconstruction, filters, gates and
settings are committed with the [freeze manifest](../configs/esdiva_phase2_freeze.json).
No experimental or ULMShare result is part of this tier.

| Scene | Tissue | Blood speed multiplier |
|:--|:--|--:|
| stationary_fast | Stationary | 1.00 |
| rigid_fast | Rigid axial/lateral oscillation | 1.00 |
| smooth_fast | Smoothly varying displacement | 1.00 |
| smooth_slow | Same smooth displacement | 0.02 |
| smooth_zero | Same smooth displacement | 0.00 |

The underlying vessel-axis peak is 0.45 m/s, so 0.02 corresponds to a
9 mm/s peak and 0.00 to stationary *blood scatterers*. The zero-flow scene
still contains blood scatterers; the separate tissue-plus-noise acquisition is
the no-blood counterfactual. Rigid tissue moves up to 0.06 mm axially and
0.30 mm laterally. Smooth motion varies those shifts with position. Tissue is
stationary during the disjoint 16-frame calibration segment; the 80 evaluation
frames span the motion cycle. Every frame fires the same `0°, -5°, +5°` events
at 3.7 kHz frame PRF and 11.1 kHz emission PRF. Only 10% of the upstream
particle cloud is retained: 456 tissue and 618 blood scatterers per seed.
**This reduced-density simulation cannot establish full-density or phantom
performance.**

The adapter acquires isolated tissue and blood channel RF using the same
trajectories and transducer events, verifies one independent combined RF event,
then adds a fixed 6 dB receiver-noise realization. Beamforming preserves
three uncompounded views. The carrier is measured on calibration RF, followed
by identical IQ demodulation and a common 8×16 tile near 15.6 mm. Each method
gets mixture IQ, tissue-plus-noise calibration IQ, the geometric outer
`|x| ≥ 6 mm` training ROI and acquisition geometry. The isolated sources and
analytic scatterer displacement stay with the evaluator. The latter is **not
an exact beamformed IQ warp**, so its pixel RMSE is a geometric diagnostic.

The RF/beam/IQ independent-sum relative errors must each be at most 2×10⁻⁴;
the training ROI blood/tissue power ratio must be at most 0.1. The input gate
requires coherence at least 0.70 on both halves of tissue calibration and
scalar-normalization noise power gain at most 4 for each view. Checkpoints and
raw RF/IQ live outside Git; compact reports record checksums and gates.

All filters score only frames 65–80 (zero-based indices 64–79), in a 16-frame
packet with a four-packet, 64-frame maximum history. Settings are in the
[shared config](../configs/esdiva_phase2_development.json): identity and
roundtrip controls, global SVD discard 1, Baranger SSM, Axis 1:
Motion-compensated virtual long ensemble, Axis 2: Motion-demodulated multi-lag
separation, Axis 3: Local ultrasound STAP (orthogonal and guarded loaded),
and Stephen Lee AER with the previously supplied temporal ε=26. These
operating points are fixed comparisons, not individually optimized claims.

Loaded STAP has 48 angle-time coefficients. Its spatial guard excludes a
one-pixel Chebyshev radius around each test pixel. The data-only gate requires
at least 96 nominal guarded samples and maximum operator gain at most 1.5;
it substitutes orthogonal STAP on failure. Four packets supply a nominal
minimum of 124 samples; spatial and temporal correlation mean **124 is not
124 independent observations**. Relative covariance loading is 3.0. The
motion gate uses disjoint spatial registration and validation pixels and a
0.12 maximum validation residual; rejected estimates produce exactly zero
displacement. A development diagnostic with a 0.30 RF-scale residual ceiling
is recorded separately and is excluded from the frozen comparison because it
increased paired blood-response error on the rigid scene.

The primary two outcomes are (1) IQ NRMSE of the separately refitted mixture
minus separately refitted tissue-plus-identical-noise output against isolated
blood, and (2) tissue-plus-noise output power divided by input tissue power.
The first penalizes attenuation and distortion of the added-blood response;
the second measures background leakage. Frozen source attribution, noise
power, operator gain and displacement error are diagnostic, not a substitute
for this pair. Do not promote an adaptive STAP or motion setting that worsens
either primary outcome on development. Publish both passed and failed strata.
On the reserved seeds, report the same locked settings and all outcomes;
do not tune on a reserved output.

Reproduce development acquisitions with
[`esdiva_phase2_acquire.py`](../examples/esdiva_phase2_acquire.py) for each
scene, then score with [`esdiva_phase2_score.py`](../examples/esdiva_phase2_score.py).
After the manifest is committed, the
[`esdiva_phase2_run_heldout.py`](../examples/esdiva_phase2_run_heldout.py)
entry point checks all frozen file hashes, the reserved seed and scene before
acquiring RF. It checks the acquisition provenance again before scoring.
