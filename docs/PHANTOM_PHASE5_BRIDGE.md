# Conditional phase-5 local-phase to physical phantom bridge

**Preregistered analysis mapping, 25 October 2026; no physical acquisition.**
The machine-readable [bridge contract](../configs/phantom_phase5_bridge.json)
is separate from the frozen [v07 physical protocol](PHANTOM_PROTOCOL_V07.md).
The latter evaluates a v06 geometric motion estimator on 8×16 IQ. Phase 5
uses a 32×16 tissue phase fit, an exact central 8×16 score tile, 16
stationary calibration frames and a unit-magnitude correction. The phase-5
held-out comparison failed; this bridge **does not release it for phantom
validation**. A phase-6 candidate must first pass a fresh simulation reserve
and its code/configuration must be hashed into a physical acquisition card
before the first physical capture. If that never happens, this protocol
remains a readiness specification, not a result.

| Property | eSDIva phase 5 | Physical bridge, fixed before capture |
|:--|:--|:--|
| Raw input | Three source-resolved simulated RF views | Measured raw RF/channel-IQ retained with hashes; unfiltered complex beamformed three-view IQ enters the filter |
| Geometry | Transmits 0°, −5°, +5°; fixed 0.1 mm axial / 1 mm lateral grid | Measured transmit/receive angles, offsets, same-view PRF, wave speed, IQ carrier and grid entered in the card; transmit angles within 0.5° of 0°, −5°, +5° |
| Spatial support | Wider 32×16, central rows 12–19 score as 8×16 | One contiguous 32×16 tile; core is the exact array slice `wide[:,:,12:20,:]`, with no independent interpolation or crop drift |
| Calibration | Sixteen prior stationary tissue frames, first eight fit and last eight validate | Acquire 64 pump-off, motion-rest frames per v07; **use only frames 0–15** for phase-5 fit and validation; retain all 64 for v07's separate frozen analysis |
| Evaluation | 80 frames, 64 history, last 16 scored | 256 frames, 64 warm-up, score frames [64,256) in twelve 16-frame packets; the physical port must support this longer sequence without altering method settings |
| Tissue support | Hard-coded tilted-vessel guard 5 mm | Preacquisition structural image and measured PSF determine one fixed `registration_roi`, disjoint `tissue_score_roi`, `flow_roi` and one-pixel flow guard; hash masks before development capture |

The 32×16 physical tile must contain ≥128 registration tissue pixels, ≥40
core tissue pixels, and ≥80/40 disjoint fit/check pixels under the frozen
checkerboard assignment `(z_index+2*x_index)%3`. Score tissue and flow ROIs
must each have ≥8 pixels and must not overlap registration pixels. The flow
ROI has at least a one-pixel guard. The optical-to-ultrasound transform, probe
pose, axial/lateral pixel centers in **millimetres**, channel cross-section,
PSF measurement, vessel mask and chosen crop are entered before capture.
The exact tissue region cannot be inferred from a filter output, blood source
map or an optical motion estimate.

The local phase port replaces only Example 22's hard-coded vessel/grid check
with these fixed physical coordinates and structural masks. It keeps
first-half calibration, the disjoint spatial fit/check split, the selected
phase-6 fraction/width/gain criterion, unit-magnitude IQ correction and
independent mixture/background fits. Smoothing widths in pixels refer to the
measured physical grid and their resulting axial/lateral millimetres are
recorded; they are **not** silently converted from the eSDIva grid. Before
capture, a numerical equivalence test must reproduce the phase-6 eSDIva IQ
core exactly when given its published grid and masks, and the port source
hash must be sealed in the card. If the physical scanner cannot meet the
three-view, tile or reference requirements, make a new preregistration
version before acquiring the held-out physical session.

## Physical responses and independent truth

Use v07's independently acquired 3 motion conditions × 2 mean flow speeds
(3 and 20 mm/s) × 5 repeats in each of development and fresh-remount held-out
sessions. Each repeat retains tissue-only, seeded-channel pump-off and
seeded-channel pump-on captures with balanced order, independent flowmeter
and optical/PIV reference, plus optical fiducials/DIC for tissue motion.
Match the stage program and speed within each triplet. Apply each filter to
each capture **independently**; never subtract complex IQ from different
captures. Reference fields `(frame,2,32,16)` and their valid support are
evaluation-only. They stratify and align outcomes and diagnose residual
tissue motion. A phase-only correction is not a geometric displacement
estimate, so the v07 displacement-RMSE release criterion is not transferred
to phase 6.

In the fixed structural tissue ROI report tissue-only output/input power
leakage. In the flow ROI report v07's pump-associated Doppler-band excess
retention `(P_on,out-P_off,out)/(P_on,in-P_off,in)` in the preregistered slow
and nominal frequency bands. It is an *operational flow response*, not exact
blood IQ preservation. Source-resolved blood NRMSE and the simulated paired
roundtrip cannot be computed from these independent physical captures;
report the available unit-modulus and interpolation-only controls instead.

The physical phase-6 candidate, if simulation-qualified, must satisfy all
six held-out motion×speed cells. Whole-triplet 2,000-resample 95% confidence
bounds require mean pump response retention lower bound ≥0.8 and mean tissue
leakage upper bound ≤0.1. In moving cells, local phase versus no phase must
also have tissue leakage ratio upper bound ≤0.9 and pump retention difference
lower bound ≥−0.02. Treat a cell as inconclusive when input pump excess is
≤0.5 of pump-off band power or fewer than six FFT bins fall in its band;
inconclusive, ineligible or missing cells block validation. Report every
repeat, method, control and static/rigid/nonrigid stratum. SVD, Axis 1,
Baranger SSM, Stephen Lee AER, Axis 2 and **Axis 3: Local ultrasound STAP**
remain fixed comparators, not a selection menu after physical capture.

## Remaining preregistration fields

Before **any** physical capture, seal a dated card with measured scanner,
probe, waveform, receive aperture, actual geometry and beamformer hash,
structural masks and their hashes, PSF/crop, optical calibration/uncertainty,
stage program, flowmeter calibration, medium, clocks and synchronization.
The card must also contain the successful phase-6 simulation reserve decision
hash and the exact port source/configuration hashes. The frozen v07 protocol
and its original 180-slot scaffold remain intact. The phase-5 bridge
preflight refuses a missing or unconfirmed candidate/card; no synthetic
substitute makes a physical cell eligible.
