# Frozen follow-up protocol, 23 October 2026

Axis labels follow the [names in the research note](REFERENCES.md#collaboration-note-and-proposed-research-axes).

## Scope and correction to the previous interpretation

The objective remains tissue rejection while preserving blood information.
The v0.2.0 frozen-component scores depend on which fitted operator is frozen.
For X=U S V^H, U_b U_b^H X and X V_b V_b^H give the same filtered mixture,
but generally act differently on an isolated blood component B. The same issue
arises for spatial versus temporal AER, even when covariance normalization is
matched. These conditional attributions cannot establish an intrinsic advantage
of one implementation of an otherwise identical adaptive filter.

The new primary endpoints are the actual output-to-blood IQ error and the output
power of a separately refitted background-only run F(T+N), divided by input
tissue power. The latter is a background-response diagnostic, not measured
tissue leakage inside F(X). Refitted blood injection responses
[F(T+N+alpha B)-F(T+N)]/alpha, alpha=0.1 and 1, measure incremental behavior;
they do not define an additive decomposition. Frozen component scores remain
secondary, with explicit projection-side/equivalence controls.

## Synthetic operating-point calibration and held-out evaluation

- Configuration: `configs/fidelity_followup_v04.json`, written before execution.
- Seven conditions: slow spectral overlap, accelerating tissue, a texture jump,
  opposite flow direction, padded integer image translation, deforming PSF
  transport, and an equal-static-source failure control.
- Each input has 256 frames, 8 by 12 pixels, 1-kHz revisit and 5-MHz center.
- Source calibration seeds: 11 and 23. Held-out evaluation seeds: 61, 73 and 89.
- Grids compare short/long/global/local SVD, the PDF-audited literature methods,
  Axis 1 memory/registration/centering ablations, original and guarded Axis 2,
  and Lee temporal/spatial centered ZCA. Absolute Lee regularizers 1, 26 and
  1000 retain their input units; relative regularizers and unit-operator scaling
  are separately labeled extensions.
- Target background/tissue power ratios: 1e-2, 1e-4 and 1e-6. Within each family
  and condition, select the lowest mean calibration output-to-blood error among
  settings that meet the target on BOTH calibration seeds. Apply that exact
  setting to all three held-out seeds. Report infeasibility and held-out target
  violations; do not interpolate, fit gains to truth, or replace failed points.
- This selection uses source truth on calibration simulations. It is explicitly
  calibration-assisted, not a deployable data-only operating-point selector.
- Every adaptive filter sees only its own input mixture. Source masks and known
  components are evaluation-only. No evaluation truth selects parameters.
- Also record the former frozen-attribution matching as a secondary analysis,
  to expose its dependence on temporal versus spatial operator representation.

Axis 1 without centering is invariant to global phase demodulation/restoration
because scalar phase cancels in each spatial outer product. Its no-memory
version equals local SVD at matched support. Registration, covariance memory
and mean subtraction must therefore be assessed separately. The translated
case uses padded crops with no wrapping or hidden evaluation-boundary removal.

Axis 2 varies rank truncation and Tikhonov whitening regularization. A separate
extension projects the reduced reconstruction operator onto a spectral-norm
ball of radius one. This bounds that local frozen operator, not the full
adaptive response or every centered/motion-compensated pipeline. Rejection,
preservation, conditioning, original gain, capped gain and lag signatures are
recorded independently. No source-identification guarantee follows from capping.

## Longer experimental pilot

`configs/openhrf_followup_v04.json` fixes two Waterloo CCA baseline acquisitions:
Acq0, anonymous subject 1, and Acq12, anonymous subject 3. Metadata were inspected
before processing; no experimental filtered outputs were inspected to select
these files or operating points. Both identify non-lite us4r and full 128-channel
data. Each file's transmit delays, angles, sample rate, event clock and aperture
are checked independently against the reconstruction contract.

Decode frames 0:1024 in 16-frame RF batches, reconstruct on the original grid
subsampled by two, and evaluate every second reconstructed grid point. Frames
0:256 estimate a bright-region nuisance ROI and bulk velocity. Frames 256:1024
are the disjoint 768-frame evaluation segment; absolute event timing and both
views remain intact. At 3 kHz this evaluates about 0.256 seconds per subject,
24 times the earlier 32-frame evaluation duration. It is not a whole-cycle,
long-ensemble ULM, or population-level experiment.

Compare fixed SVD/literature/Lee AER/Axis 1/Axis 2 arms and matched-view STAP
controls. No experimental oracle or source-separated blood/tissue truth is
created. Bright-region power, spectra, coherence and image structure are proxies;
they do not establish blood preservation. Calibration and evaluation packets
from one recording remain physiologically correlated. Contrast-injection status
and quantitative experimental phase calibration remain unconfirmed.

Raw RF and generated IQ arrays stay outside Git. Portable reports retain the
source revision, selected frames, geometry, subjects, input/code/config hashes,
method settings, errors and numerical diagnostics. Analysis scripts and figures
are versioned. Do not infer a universal ranking from pooled conditions.
