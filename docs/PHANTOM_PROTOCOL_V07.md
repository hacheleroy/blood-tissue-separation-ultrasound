# Preregistered phantom validation: motion and blood/tissue separation

**Status: analysis frozen, acquisition not yet performed.** This protocol
follows the [v06 synthetic software gate](RESULTS_MOTION_REPAIR_V06.md). The
numerical filter source stays at SHA256
`8587a156b3283238920777e96cd059d72203d5a9638a61f5cc43876324cd4e13`.
`configs/phantom_v07.json`, `scripts/phantom_v07.py` and this text are sealed by
`PHANTOM_FREEZE_V07.json`. The hardware acquisition card is the remaining
physical preregistration: enter actual scanner geometry, stage trajectory,
optical-to-ultrasound calibration, medium and pump specifications and hash the
card **before the first capture**. Do not substitute a convenient public
dataset for any of the controls below.

## Physical experiment and independent references

Use beamformed **unfiltered complex angular IQ**, three measured views at
their true same-view revisit rate, one fixed 8×16 tile per capture. Keep the
transmit/receive angles, waveform, PRF, gain, reconstruction grid, probe pose,
tile position, scatterer concentration and beamformer fixed across all arms.
Retain the raw RF/channel-IQ SHA256 for every capture. A separately measured
flowmeter targets **mean lumen speeds of 3 and 20 mm/s** from independently
measured flow rate and channel area (value within ±20%). Record the area and
calibration; optical PIV centerline speed is a separate check. No contrast microbubbles in this
first test. A transparent acoustic/optical phantom with external camera/PIV
access is suitable; optical PIV can provide independent fluid-velocity
information while stage encoders or optical fiducial/DIC tracking provide
**tissue displacement**. They are distinct references and must not be inferred
from the same ultrasound IQ used by the filter. The literature supports
dual-modality phantom feasibility and also shows that ultrasound Doppler and
optical flow can disagree substantially
([dual-modality phantom study](https://pubmed.ncbi.nlm.nih.gov/39752855/),
[transparent PVA phantom study](https://pubmed.ncbi.nlm.nih.gov/36878829/)).

Acquire a 3×2×5 factorial set independently in a **development** session and
a subsequent **held-out** session with a fresh phantom remount/refill or
independent setup. For each repeat, use the same programmed stage motion and
speed in three captures:

| Motion condition | Physical tissue reference | Captures at each speed |
|:--|:--|:--|
| Static | Optically verified near zero | tissue-only; seeded channel, pump off; seeded channel, pump on |
| Rigid | Subpixel translation via stage encoder plus optical fiducials | same three captures |
| Nonrigid | Deformable layer with independently measured dense optical/DIC field | same three captures |

Record the three captures in randomized order where the plumbing permits; the
pump-on/off order is balanced over repeats. Record the ordering and any refill
or temperature drift. Each capture has **64 prior pump-off calibration frames**
while tissue is in its undeformed reference state, followed after stage/pump
settling by **256 acquisition frames**. The analysis warms up on frames
0–63 and scores frames **[64,256)** in 16-frame packets with 64-frame history.
Calibration may be separated by a settling gap, but must never overlap or
follow evaluation. A common-reference tile/masks and the same calibration
policy apply to all methods. A mounted flow channel can contribute wall/tube
echoes; the pump-off capture measures that baseline. The tissue-only capture
has no flowing scatterers and provides the leakage negative control.

An independent optical/stage reference field has shape `(frame,2,z,x)` and
means **reference pixel → native/lab pixel**. Sync it to ultrasound timestamps,
measure optical-to-ultrasound calibration independently, and report an upper
uncertainty ≤0.1 pixel on the evaluated interior. Rigid/nonrigid fields must
have RMS ≥0.2 pixel and peak ≤1.25 pixels; nonrigid spatial variation RMS
≥0.05 pixel. The static optical RMS must be ≤0.1 pixel. These are eligibility
checks, not a means to select filters. Reconstruct reference fields only where
fiducial/optical support exists; no extrapolated displacement is scored.

Draw `registration_roi`, held-out `tissue_score_roi`, and `flow_roi` from a
separate preacquisition structural image in **reference coordinates**. Fix
them before the development session, hash the mask file, keep them across
sessions, and exclude at least one pixel around the flow ROI from tissue
ROIs. Require at least 25, 8 and 8 pixels respectively. No blood/source mask
or optical field enters `PublicInput`; the tissue score mask is disjoint from
registration and STAP training. Optical fields align **both input and output**
only when evaluating scores in the common reference frame. Filter output is
also reported in its native coordinates and the interpolation-only roundtrip
control remains visible.

## Fixed methods and measurement endpoints

The 27 fixed entries in the JSON configuration are v06's 64-frame uniform
SVD, weighted SVD, **Axis 1: Motion-compensated virtual long ensemble**,
Baranger SSM, Stephen Lee AER (temporal and spatial), **Axis 2:
Motion-demodulated multi-lag separation**, and **Axis 3: Local ultrasound
STAP**, each with no registration, gated rigid and gated smooth registration,
plus identity/roundtrip controls. The v06 `MotionGate` parameters are copied
exactly. Do not retune rank, AER epsilon, loading, gain cap, ROI, mode or
frequency band on the held-out session. The candidate gate is **Axis 1 gated
smooth, cutoff 1**. SVD and the other families are matched comparators,
not a menu for post hoc selection. No known-warp filtering path is run on
phantom data; optical motion is evaluation-only.

For a spatial ROI $R$, take the mean across views and reference-aligned
score frames. Define tissue leakage

$$L=\frac{\langle|F(X_{\mathrm{tissue\ only}})|^2\rangle_R}
{\langle|X_{\mathrm{tissue\ only}}|^2\rangle_R}$$

For the **predeclared absolute Doppler-frequency band** expressed as cycles
per same-view frame, slow `[0.005,0.05]`, nominal `[0.02,0.25]`, let $P_B$
be the Hann-windowed, demeaned complex-IQ spectral power integrated over
both Doppler signs. In the flow ROI, define the operational pump-associated
response retention

$$R_B=\frac{P_B(F(X_{\rm on}))-P_B(F(X_{\rm off}))}
{P_B(X_{\rm on})-P_B(X_{\rm off})}$$

The two captures are **independently filtered**. Never subtract their complex
IQ or call $R_B$ exact blood-signal retention: independent scatterer phases,
tube/background changes and pump-induced motion prevent that interpretation.
If input excess is ≤0.5 times the pump-off band power, or if fewer than six
FFT bins fall in the registered band, mark the endpoint **inconclusive**. Do
not tune bands on the held-out acquisition. Alongside both ratios, report
per-frame gate counts, interpolation roundtrip loss, and estimated-versus-
independent-reference displacement RMSE on the valid interior; compare
no compensation, rigid and smooth in the same capture. Optical/PIV velocity
agreement and flowmeter stability are secondary task-level checks, not source
IQ truth. Report every individual repeat and speed/motion stratum before any
pooled result.

## Inclusion, decision and data contract

The `phantom_v07.py preflight` command requires both independent sessions,
every 3×2×5 paired triplet and all three capture files, hashes, clocks and
ROI/reference support. Missing captures, ROI leakage/overlap, clock mismatch,
out-of-range motion, excessive optical uncertainty, changed geometry or
beamformer, and flowmeter mismatch fail eligibility. Do not fill gaps with
synthetic data. Quality exclusion is made from inputs and independent
references before examining method outputs. An ineligible stratum stays
inconclusive rather than being silently dropped.

For each held-out motion×speed cell (five paired repeats), use a fixed-seed
2,000-resample percentile bootstrap over **whole triplets**, not frames. The
Axis 1 gate passes only if all six cells are identifiable, the lower 95% bound
of mean $R_B$ is ≥0.8, the upper bound of mean $L$ is ≤0.1, and the moving
cells have upper mean displacement RMSE ≤0.25 pixel and mean RMSE at least
20% below the no-compensation field. The static cells require the upper mean
accepted-frame fraction ≤1%. All criteria must pass; a failed or incomplete
cell blocks an experimental validation claim. The development session is for
hardware troubleshooting and a separate report; changes to frozen analysis
require a new protocol version **before** acquiring a fresh held-out session.

Manifest and capture paths stay local and are ignored by Git. A manifest has
`protocol_sha256`, an `acquisition_card` `{path,sha256}`, and exactly two
`sessions` with distinct `id`, `split`, card-matching `geometry`, a hashed
`masks` NPZ and 30 `triplets` each. Every triplet specifies `motion`, `speed`,
`replicate` (1–5), and `files` for `tissue_only`, `flow_off`, `flow_on`, each
with `{path,sha256,rf_sha256,beamformer_sha256,reference_source,
calibration_pump_off,calibration_motion_rest}`; `flow_on`
also includes `flowmeter_mps`. Each capture NPZ (`allow_pickle=False`) has
`views` `(256,3,8,16)` and `calibration_views` `(64,3,8,16)` complex,
uniform `timestamps_s`/`calibration_timestamps_s`, independent
`reference_displacement` `(256,2,8,16)`, Boolean `reference_valid`
`(256,8,16)` and scalar `reference_sigma_pixels`. The card includes the
measured `geometry`, beamformer hash, preregistration/capture times, and
nonempty `scanner`, `probe`, `transmit_waveform`, `receive_aperture`,
`optical_registration`, `stage_program`, `pump_medium`, and
`flowmeter_calibration` records. Capture hashes cannot be reused across
repeats. Nonrigid fields require dense independent optical DIC, rather than
a stage encoder alone. The optical reference is **not** a model
input. The output report contains per-triplet metrics, eligibility evidence,
source/config hashes and a transparent gate decision, but no raw acquisition.

```bash
python scripts/phantom_v07.py check-protocol
python scripts/phantom_v07.py scaffold --output /path/to/phantom_v07
python scripts/phantom_v07.py preflight --manifest /path/to/phantom/manifest.json
python scripts/phantom_v07.py run --manifest /path/to/phantom/manifest.json --output /path/to/phantom/report_v07.json
```

The scaffold writes 60 declared triplets (180 capture slots) and a deliberately
incomplete card. Fill the card and its SHA256, copy its measured geometry into
both sessions, populate each capture's checksums/reference/flowmeter fields,
and hash the fixed mask NPZ before preflight. A null placeholder is never an
eligible datum. Preserve the preacquisition card and its dated revision in an
independent lab record; a self-reported timestamp or hash alone cannot prove
when the decision was made.

The existing [Open H-RF inventory](DATASETS.md) has no verified bubble-free
phantom set with this complete independent source/control contract. Until
one is acquired, this is a frozen **plan and executable analysis**, not an
experimental result. The objective is clutter rejection with preservation of
the flow response; Doppler, fUSi, ULM and WSS remain downstream applications.
