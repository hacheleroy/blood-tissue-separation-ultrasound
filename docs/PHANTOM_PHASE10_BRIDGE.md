# Phase-10 simulation to physical phantom bridge

**Acquisition preparation; no physical measurements yet.** The
[frozen reserve comparison](ESDIVA_PHASE10_RESERVE_RESULT.md) passed all
eight scenes, allowing this preparation to proceed. Its
[machine-readable analysis contract](../configs/phantom_phase10_bridge.json)
binds the existing numerical candidate, physical response thresholds and
control matrix. This bridge
supplements the independent-control [v07 physical protocol](PHANTOM_PROTOCOL_V07.md).
The [phase-10 reserve decision](../runs/esdiva_phase10_reserve/decision.json)
passed under its [prior freeze](ESDIVA_PHASE10_RESERVE_FREEZE.md) before
this candidate was described as ready for physical transfer. The older
[phase-5 bridge](PHANTOM_PHASE5_BRIDGE.md) and the v07 analysis source are
historical, separately frozen protocols; neither is a phase-10 physical
implementation.

## What must be ported and checked before acquisition

The current phase-10 numerical code calls Example 22's `tissue_mask`, which
requires a 32×16 grid with axial positions 14.4–17.5 mm at 0.1 mm spacing,
lateral positions −8…7 mm at 1 mm spacing, and a known 60° vessel axis.
`phase10_support`, `fit_local_phase`, `fit_shifts`, `select_mode`, the
reference transport and scored-core checks inherit that geometry. Applying
these functions directly to a physical scanner crop will fail or impose a
false vessel model. The simulated guard is **lateral distance from the vessel
axis**, not distance from the lumen wall. The 5 mm wide and 5.5 mm core guards
must use an independently measured structural centerline on the physical
grid; the measured wall/boundary must additionally confirm that those pixels
exclude the lumen and flow channel. Their physical support and fit/check counts
must be recorded, not silently forced to Example 22's 192/40 pixels.

The separate [geometry-aware port](../examples/phantom_phase10_port.py) and
[equivalence runner](../scripts/phantom_phase10_equivalence.py) implement the
same candidate with a measured vessel-axis argument. The [engineering record](../runs/phantom_phase10_port/equivalence.json)
uses two previously seen 2333 cases and one previously seen 2309 case: the
structural masks, calibration and selector audits matched, the hybrid and
reference choices were both exercised, and all 54 separately fitted
mixture/background roundtrip/SVD/Axis 1 outputs had zero recorded relative
difference. The paired added-blood errors also matched to 1e−12. This is an equivalence
check on simulated IQ, not a physical test. Before the **first physical
capture**, freeze the port and its analysis configuration in the measured
card. On the published eSDIva grid, vessel axis, IQ and calibration, the
port must continue to reproduce the
phase-10 core IQ, masks, mode choice, paired roundtrip and method outputs to
the documented numerical tolerances. The port must use only unfiltered
complex IQ, measured scanner geometry, fixed structural masks and prior
calibration; no optical motion, source truth or pump state can enter filtering
or mode selection. Hash its source and configuration in the physical
acquisition card. A physical geometry that cannot provide three true views,
the required wide/core support and disjoint fit/check pixels requires a new
preacquisition protocol, not a post hoc replacement.

| Item | Simulation gate | Physical acquisition mapping |
|:--|:--|:--|
| Views | Three `0°, −5°, +5°` nominal views | Three measured angles within 0.5° of nominal in that order, delays and same-view revisit rate; record actual deviations and receive aperture |
| IQ | 32×16 wide, exact central rows 12–19 core | Beamform three unfiltered complex wide IQ views; crop core by the identical array slice without a second resample |
| Prior calibration | 16 stationary frames; first 8 fit, last 8 validate | Acquire v07's 64 pump-off motion-rest frames, use 0–15 for phase 10, preserve all 64 for separate v07 analysis |
| Time | 80 evaluation frames; first 64 select, last 16 score | Acquire 256 frames; apply phase 10 independently to twelve 80-frame windows starting at 0, 16, …, 176 and retain only each last 16-frame output |
| Tissue masks | Fixed simulated vessel-axis guards | Preacquisition structural centerline and wall, physical-grid wide and core lateral guards, disjoint fit/check and independent tissue/flow score ROIs; hash the mask and crop |
| Source scores | Isolated simulated blood/tissue paired IQ and source purity | Tissue-only leakage and paired pump-associated Doppler response; never label these exact blood IQ retention |

## Matched physical observations

Retain v07's independent development and fresh-remount held-out sessions:
static, rigid and nonrigid tissue at independently measured mean flow targets
3 and 20 mm/s, five repeats per motion×speed cell. Each repeat has
tissue-only, seeded-channel pump-off and seeded-channel pump-on captures,
with randomized or balanced capture order. Store raw RF/channel IQ and its
hash, the unfiltered beamformed complex views, scanner and beamformer version,
acquisition clocks, actual angle/offset/grid, temperature/medium, channel
area and flowmeter calibration. Measure tissue displacement with stage
encoder plus optical fiducials for rigid motion and dense optical DIC for
nonrigid motion. Use optical/PIV flow and a calibrated flowmeter as independent
fluid references. Keep optical fields and flow truth evaluation-only.

The optical-to-ultrasound coordinate transform, measured PSF, wide crop,
vessel centerline, wall and scored ROIs must be fixed and hashed from a separate
structural capture before the development session. The flow and tissue ROIs
are disjoint with at least a one-pixel guard, and the tissue score ROI must
not supply motion training pixels. Require at least 192 wide and 40 core
tissue pixels, 80/40 disjoint wide fit/check pixels, and 8 tissue-score and
8 flow pixels. Confirm that the
same physical mask and beamformer remain valid after remount, or preregister
a geometry transformation and its uncertainty before the held-out session.

Each capture receives its own phase-10 fit and filter output. In the fixed
tissue ROI, report tissue-only output/input power. In the flow ROI, report
v07's independently filtered pump-on minus pump-off Doppler-band excess
retention in the predeclared slow/nominal bands. This is an operational flow
response and may differ from source-resolved blood preservation. Report
roundtrip error, mode decisions, calibration coherence, optical displacement
agreement, flowmeter stability and individual repeat outputs. The candidate
is the frozen phase-10 selected-mode Axis 1 cut 1, with same-scene no-motion
Axis 1 and SVD cut 1 controls. Retain the common SVD, Baranger SSM, Stephen
Lee AER, Axis 2: Motion-demodulated multi-lag separation, and Axis 3: Local
ultrasound STAP comparisons. No method, mode or threshold is selected from
held-out physical outputs.

Transfer v07's whole-triplet bootstrap and six-cell requirement: 2,000
fixed-seed resamples over five complete triplets per cell, lower 95% mean
pump-response retention ≥0.8, upper 95% mean tissue leakage ≤0.1. Also
predeclare moving selected-versus-no-motion leakage ratio upper 95% ≤0.9
and pump-retention difference lower 95% ≥−0.02. Static suppression is
reported without a moving-tissue improvement condition. Inconclusive pump
excess (≤0.5×off-band power), fewer than six band bins, missing independent
motion/flow references or incomplete paired captures block a validation
claim. A development-session change requires a new frozen port and a fresh
held-out session before testing it.

## Acquisition card and handoff

The card must be dated and sealed **before first capture** with the successful
reserve decision/manifest hashes, physical port and analysis hashes, scanner,
probe, transmit waveform, receive aperture, measured three-view geometry,
same-view PRF, wave speed and carrier, beamformer hash, wide/core crop and
structural centerline/wall/mask hashes, PSF, optical transform and uncertainty,
stage program, pump, phantom medium, channel area, flowmeter/PIV calibration,
clock synchronization, randomization, capture IDs and the planned two
independent sessions. A self-reported timestamp is not independent proof of
preregistration; retain the dated card in a lab record. This repository has
neither that measured card nor the independent physical RF and reference
captures. Until they exist, the bridge is an acquisition specification, not
a measured phantom result.

After a complete reserve pass, generate an **incomplete** card and a 180-slot
capture inventory (two sessions, 30 triplets per session, three matched
captures per triplet):

```bash
python scripts/phantom_phase10_scaffold.py --output /path/to/phase10_phantom_plan
```

Every measurement, RF/IQ/reference hash and physical port hash is initially
null. The command deliberately refuses a failed or incomplete reserve bank.
Filling its placeholders and independently checking the geometry-aware port
is a separate step; generating the inventory does not authorize capture or
establish a phantom outcome.
