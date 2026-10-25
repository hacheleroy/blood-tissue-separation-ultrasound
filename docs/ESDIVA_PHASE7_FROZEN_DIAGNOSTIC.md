# Phase-7 reference-transport diagnostic, frozen before new RF

25 October 2026. This is a **single fixed diagnostic candidate**, not a
qualified replacement for the failed phase-6 comparison. The declared new
development clouds are 2503 and 2521; prospective 2333 and fresh reserve
2401/2411 remain unopened. The [split](../configs/esdiva_phase7_split.json),
[configuration](../configs/esdiva_phase7_candidate.json),
[transport source](../examples/esdiva_phase7_reference_transport.py),
[source-resolved scorer](../examples/esdiva_phase7_score.py) and
[all-scene reviewer](../examples/esdiva_phase7_review.py) are fixed before
2503/2521 RF generation. The config locks their SHA-256 hashes.

The model fits a bounded translation on disjoint tissue magnitude pixels of
the 32×16 wide tile. It samples only the stationary **complex calibration
reference** at translated coordinates, masks outside-the-measured-tile and
weak-reference locations, and validates phase gain on disjoint tissue pixels.
The complex IQ observations stay on the original grid. The proposed phase is
multiplied by the fixed strength 1.10 and applied as a unit-modulus factor to
the exact central 8×16 core. Inverse multiplication makes a roundtrip exact;
neither the simulator's displacement nor the isolated sources enter the fit.
The boundary is explicit, with zero correction where the reference lacks
support. This diagnostic does not claim a physical edge model.

## Seen development response

The fixed candidate was scored on all four scenes of the already seen 2309
and 2311 clouds before new RF. Its paired response relative to the unchanged
no-motion filter is archived in [eight compact scores](../runs/esdiva_phase7_design/).
The [partial review](../runs/esdiva_phase7_design/partial_decision.json)
records the expected eight missing new-development scenes
and the old-development failures. The principal moving-scene responses are:

| Seed / scene | SVD blood NRMSE increase | SVD tissue ratio | Axis 1 blood increase | Axis 1 tissue ratio |
|:--|--:|--:|--:|--:|
| 2309 rigid fast | 0.0180 | 0.8494 | 0.0149 | 0.8979 |
| 2309 smooth slow | −0.0324 | 1.0480 | −0.0257 | 1.0539 |
| 2311 rigid fast | 0.0210 | 0.8829 | 0.0186 | 0.8844 |
| 2311 smooth slow | 0.1117 | 0.8401 | 0.0600 | 0.8509 |

Thus the rigid tissue benefit is real on these two development clouds, but
the candidate already fails the preserved ≤0.020 blood increment and ≤0.900
moving-tissue ratio gates. A calibration phase-validation gain alone cannot
license it. Alternative residual-phase, affine-magnitude, tissue-only and
noncyclic unitary transport probes on these seen clouds did not provide a
complete paired pass. The frozen diagnostic is used only to measure
generalization of the failure mode on 2503/2521, without retuning.

The reviewer requires the unchanged paired gates on **all sixteen** old and
new development scenes. It distinguishes a missing scene from a negative
complete study. A failed old or new scene blocks the path to 2333; neither
the phase-7 diagnostic nor the phantom geometry bridge can release a reserve.
For slow flow, the fixed Axis 2: Motion-demodulated multi-lag separation
control stays at rank 4 and gain cap 0.4. The no-motion records retain
Baranger SSM, Stephen Lee AER and Axis 3: Local ultrasound STAP.
