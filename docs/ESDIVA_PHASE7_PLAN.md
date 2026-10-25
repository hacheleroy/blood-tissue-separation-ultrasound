# Phase 7: boundary-aware complex tissue transport

25 October 2026. [Phase 6](ESDIVA_PHASE6_DESIGN_RESULT.md) found no candidate
that passed both development seeds. The original small core gate changed from
fail to pass with seed, and cyclic IQ translation increased tissue leakage on
seed 2309 despite accurate gross motion tracking. The negative decision and
source files remain intact.

The [phase-7 seed split](../configs/esdiva_phase7_split.json) declares new
source-resolved development confirmation clouds **2503 and 2521** before
their RF is generated. Only the seen 2309/2311 development data may inform
the initial model/setting choice. Seeds 2333, 2401 and 2411 remain unopened.

## Transport model and leakage checks

Fit a tissue translation independently for each view/frame using observed
wide-tile magnitude and disjoint fit/validation tissue pixels. Predict the
*complex* calibration reference at translated coordinates inside its measured
32×16 support, with an explicit validity mask at the axial/lateral boundary.
The phase between the transported and resting references supplies a
speckle-transport term. Apply a validated unit-magnitude phase factor to the
original eight-row IQ core and invert it exactly after filtering. No cyclic
wrap, IQ interpolation of blood, isolated source, or known displacement enters
the fit. The calibration reference may be interpolated only at in-bounds
coordinates; invalid/low-reference support receives zero correction. Report
how often this fallback is used. A residual-phase alternative was explored
on seen development data and rejected before the fixed diagnostic was chosen.

Use the same paired source-resolved SVD cutoff-1 and **Axis 1:
Motion-compensated virtual long ensemble** cutoff-1 responses, exact
roundtrip and no-motion controls. Preserve the fixed **Axis 2:
Motion-demodulated multi-lag separation** slow-flow gate and **Axis 3: Local
ultrasound STAP**, Baranger SSM and Stephen Lee AER comparison records. The
phase-6 numerical limits remain: blood NRMSE increment ≤0.020 in every scene
for SVD and Axis 1; tissue leakage ratio ≤0.900 in rigid fast and smooth slow;
roundtrip ≤1e-10; fixed slow Axis 2 blood NRMSE ≤0.700 and leakage ≤0.180.
Input support, evaluator-only ROI blood contamination and source hashes must
also pass. Neither calibration coherence nor motion-displacement agreement
can replace the paired blood criterion.

On 2309/2311, compare transport against the phase-6 no-motion and local-phase
controls without selecting from one scene. Freeze *one* complete code/config
and decision rule on `main` before opening 2503 or 2521. Run all four scenes
on both new clouds, preserving any negative result. Only a complete pass
permits a separate freeze before the still-unopened prospective 2333; only a
later 2333 pass permits using the fresh 2401/2411 reserve. The conditional
[phantom geometry bridge](PHANTOM_PHASE5_BRIDGE.md) remains a preacquisition
contract and does not authorize physical validation in this phase.

The selected diagnostic and its already observed old-development failures are
recorded in [the frozen diagnostic](ESDIVA_PHASE7_FROZEN_DIAGNOSTIC.md).
