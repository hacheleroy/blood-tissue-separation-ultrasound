# Shared benchmark protocol v0.4.0 / study v05

Frozen design, 23 October 2026. The numerical implementation and JSON configuration
are frozen before the 40-case run. This protocol qualifies the new shared runner;
older studies retain their historical implementations and source fingerprints.

## Question and scope

Compare the three research directions on one common observation bank:

1. **Axis 1: Motion-compensated virtual long ensemble**
2. **Axis 2: Motion-demodulated multi-lag separation**
3. **Axis 3: Local ultrasound STAP**

The primary question is whether the actual complex output preserves the prescribed
blood field while rejecting background, at declared observation support and
training-only operating points. The experiment uses statistical angular IQ with
an imposed image-coordinate warp. Doppler and envelope transport are specified
separately. It is not RF, particle-scattering, RBC, acoustic or experimental
validation. No new empirical acquisition or source PDF is part of this archive.

## Shared information and clock

Every arm receives the same three uncompounded views, 8×16 tile, revisit PRF
1 kHz, centre frequency 5 MHz, transmit angles −18°, 0°, +18°, receive angle 0°,
and sequential offsets 0, 1/(3 PRF), 2/(3 PRF). The view dimension is retained in
every primary output. A coherent mean is a secondary companion only.

Each case has 64 calibration frames at times −64/PRF through −1/PRF, followed by
128 evaluation-stream frames. Calibration and evaluation have disjoint times
and independent noise and stochastic blood innovations, but **share the fixed
tissue texture**. The low-rank and identical-statistics cases also share a fixed
blood spatial texture; the low-rank temporal innovations remain independent.
Calibration is unwarped and lacks the evaluation-only drift/acceleration.
This is a deliberately favourable reusable reference, not an independent tissue
realisation. It must not be presented as an in-vivo calibration transfer result.

The common calibration estimates bulk tissue velocity from a fixed annotated
outer tissue ROI and per-view complex gain relative to view 0. All arms receive
this ROI and the same calibration packet. It demodulates the evaluation clock
and equalises relative gains, then restores laboratory complex phase/gain after
filtering. There is no truth-fitted gain, phase, alignment or amplitude correction.
The annotation is an acquisition input; source arrays are evaluator-only.

Packets contain 16 frames. The score interval is [64,128), after 64 warm-up frames.
Long-history arms use at most the latest 64 frames; short controls use 16.
Fits may use the complete current packet and release it at the packet end:
latency is bounded by 16 frames, with no future packet or full-sequence temporal
fit. Baseline SVD folds all views into the spatial dimension. Each algorithm can
choose a subset of the common observations; covariance dimensions and training
counts differ and are logged rather than asserted equal.

## Motion controls

All principal long-history families run with no spatial compensation and with
the same data-only smooth registration. This estimates each frame independently
against the calibration envelope, using robust Gauss–Newton with eight fixed
coefficients (translation, linear z/x terms and a sinusoidal product in each
displacement coordinate), eight iterations, ridge and coefficient bounds fixed
in code. It is a small model of smooth deformation, not a general nonrigid solver.

Additional fixed-setting arms use estimated rigid translation or a known-warp
evaluation-only control. A separate API argument is required for known motion;
it is excluded from operating-point selection. Identity and rigid/smooth/known
round-trip controls expose resampling and boundary losses. Complex bilinear
sampling has zero-filled image boundaries; inverse coordinates use 40 fixed-point
iterations. No periodic wrap, interior-only primary scoring or loss correction
is applied. Truth warp RMSE (whole tile and interior) is computed after fitting.
The field may be known exactly while the resampled complex signal remains lossy.

## Methods and controls

The configuration contains 73 settings. The long none/smooth arms include SVD
cutoffs 1, 4 and an explicit empty-output endpoint, a centred cutoff-1 arm,
Baranger SSM, weighted SVD, the three axes, and Lee AER. The short arm and
rigid/known controls are smaller declared subsets, not equally tuned families.
Different candidate counts limit claims of optimality; all settings are archived.

**Axis 1** maintains packet raw moments with exponential factor 0.75 and exact
expiration after four packets. Direct weighted SVD uses the same finite history,
packet weights, centring and cutoff. Weighted mean subtraction is tested against
the expiring moments. Their output equality is a required control, not evidence
for a distinct separation mechanism. The 16-frame version has no usable history
beyond the current packet. This controlled implementation does not include the
older prototype's quality gates, event resets or Gaussian patch blend.

**Axis 2** uses lags 1, 2, 4, 8, retained rank 4 or 8, covariance threshold 10⁻⁶,
ridge 10⁻³ of the top eigenvalue, and phase-width 0.1 or 0.3 radians/frame.
Finite-history lag covariances use pairs within each packet, with matching
endpoint means if centred. The component heuristic examines contiguous lag
pairs over the retained history, including packet boundaries. Hermitian real
and imaginary lag parts enter complex unitary AJD (40 sweeps, tolerance 10⁻⁶).
The heuristic combines near-DC coherence and spatial extent; it is not a trained
source classifier. The fitted spatial operator has spectral norm capped at 1.
This does not bound sensitivity of the full refitted nonlinear pipeline.
Convergence flags, off-diagonal residuals, signature separation, rank and gains
are retained for actual, background and weak-injection fits.

After the actual-mixture fit, evaluation-only least squares chooses real component
weights in [0,1] in that *uncapped fixed demixing basis*, separately per scored
packet. Its output error versus the uncapped heuristic separates limitations
of the basis from limitations of identification. It does not establish physical
identifiability, and is not an upper bound on every spectrally capped operator.
It never selects a primary setting or reaches a regular filtering API.

**Axis 3** uses a 16×3 angle-time snapshot and a nominal constant tissue steering
vector in the calibrated frame. Orthogonal and loaded oblique tissue-null
projectors return the full view-time IQ. Guarded spatial training uses the
shared annotated outer ROI, excludes Chebyshev distance ≤1 around each target,
and pools packet snapshots within the declared 16- or 64-frame history.
Diagonal loading is 0.01 or 0.1 of mean covariance trace. Fewer than four
training columns trigger the orthogonal fallback. The record includes K/d,
counts below 2d, fallback count and operator gain. Correlated columns are not
assumed independent. It is a small local ultrasound STAP prototype; no receive
aperture, near-field calibration or LCMV detection output is silently substituted
for preserved IQ. Joint-view SVD has access to the same calibration and views.

**Stephen Lee AER** is included explicitly, using the already audited centred
ZCA algorithms: temporal covariance denominator P (384 features), spatial
denominator H−1, inverse-square-root regularisation. Temporal absolute ε is
1, 26 or 1000, with a relative 0.001 option; spatial ε is 26 or relative 0.001.
These absolute values are in this synthetic, commonly calibrated IQ scale and
are not transferred paper-optimal values. The temporal unit-gain extension is a
separately labelled family. Shared-runner outputs match the full-history
published covariance implementation on scored columns. Frozen diagnostics use
the declared full spatial operator extension; equivalent temporal representations
can have different isolated-source attributions. This raw-IQ study does not
replicate Lee's contrast-enhanced ULM tracking experiment.

## Conditions, split and selection

Eight fixed conditions cross slow overlap (3 mm/s, 40 dB clutter), a more separated
low-rank source (20 mm/s), subpixel translation, modelled smooth deformation,
unmodelled deformation, calibration phase drift plus tissue acceleration,
opposed flow with reversed phase convention, and identical source statistics.
Except where overridden: 30 dB clutter, 20 dB SNR and blood AR coefficient 0.985.
Exact generator values are in `configs/shared_axes_v05.json`.

Each condition uses training seeds 101 and 211, and held-out seeds 307, 419 and
523: 16 training and 24 evaluation cases. The same seeds pair conditions to
reduce some stochastic differences. Conditions are not held out; evaluation
therefore tests new draws within the declared conditions, not a new domain.
A separate seed-991 pilot checks runtime, array shapes and algebra; it is not
used to choose scientific settings. No grid is changed after the frozen run.

Selection pools **all training conditions and seeds**. One setting is selected
per family × motion arm × support by lowest arithmetic mean actual-output IQ
NRMSE, subject to mean refitted background/tissue power ≤10⁻² or ≤10⁻⁴. An
unconstrained training-best setting is also reported. Infeasible groups remain
explicit; no target relaxation or held-out replacement occurs. Empty SVD output
is permitted and has blood error 1, so a strict rejection target can expose a
useless solution. A mean training constraint is not a per-case guarantee.
Each chosen ID is then applied unchanged to every evaluation condition/seed;
held-out pass fractions and errors are reported. This is truth-assisted *training*
calibration of a research comparison, not a deployment-ready data-only selector.

## Measurements and archive

Primary measurements use all scored complex view samples with physical amplitude
and phase: actual ‖F(X)−B‖/‖B‖; refitted background power
‖F(T+N)‖²/‖T‖²; and refitted incremental errors for
[F(T+N+αB)−F(T+N)]/α, α=0.1 and 1. All three fits reuse the common calibration
but re-estimate registration/covariances from their own input. Background-only
power is not the tissue attribution inside F(X); injection is not an additive
source decomposition. Actual amplitude, spectral shape, lag covariance, weak
support and pixel-tail errors accompany the aggregate. Invalid spectral shape
for an empty output is recorded, not replaced with a successful score.

Frozen source errors/retentions and linear closure are secondary, conditional
diagnostics. No reference-based adjustment improves a primary error. No velocity,
WSS or ULM claim is inferred from raw IQ alone. Figures show per-condition
held-out means/ranges over three seeds, with no population confidence interval.

Expected archive: 2,920 case-setting evaluations and 8,760 independently fitted
outputs including counterfactuals. The frozen config, all case records, source
hash, input/calibration/source hashes, fit diagnostics, selected IDs and candidate
training losses are retained in lossless XZ parts with SHA256 checks. There are
no raw acquisitions, generated IQ arrays or opaque third-party binaries.
Old milestone hashes are historical; `shared_v05_manifest.json` is the current
source-identity check. CI verifies both old archives and the new one.
Four independent case workers each use one BLAS thread. Checkpoints retain the
source/config identity, and final records are written in configuration order.

```bash
btbench shared-study --config configs/shared_axes_v05.json --output runs/my_shared_study
python scripts/plot_shared_benchmark.py runs/my_shared_study
python scripts/verify_results.py --manifest runs/shared_v05_manifest.json --require-current-code
```

This controlled study is the gate before calibrated bubble-free phantom/RF data,
independent longer in-vivo sequences and geometry-specific downstream validation.
