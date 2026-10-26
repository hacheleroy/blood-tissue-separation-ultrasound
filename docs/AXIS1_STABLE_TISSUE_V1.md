# Axis 1: stable motion and a blood-sparing tissue reference

26 October 2026. **The simple stabilized candidate fails the frozen joint
criterion.** This is a source-labeled IQ software simulation, not a physical
Resolve stroke flow phantom measurement. The 9401/9402 held-out seeds were run
once after sealing the numerical sources, configuration, seen-bank archive,
development archive and decision rule. There was no held-out retuning.

## Question and implementation

The previous [native-grid study](AXIS1_NATIVE_SUBTRACTION_V1.md) showed two
losses: blood projects onto the fitted temporal tissue mode, and independently
fitting the motion on flow-on versus flow-off mixtures changes the returned
tissue estimate. This study asks whether a **single-frame-width temporal
smoothing** of the automatic B-mode motion trajectory makes the latter fit
more stable, and whether training rank one only on bright tissue-candidate
pixels reduces the former overlap.

For each observed IQ series \(X\), form the mean across views of its
magnitude, log-compress it, and fit translations to frame zero on the brighter
half of the interior checkerboard support with a separate check half. The
existing identity/translation BIC gate and correlation failure fallback apply.
Smooth the **estimated** two-coordinate displacement through time with a
Gaussian width of one frame; subtract the smoothed frame-zero displacement so
that the reference stays anchored. The complex IQ is aligned by the smoothed
maps and the original, data-fitted per-view phase. Fit a joint-view, uncentered
rank-one temporal SVD to the aligned IQ. Return the estimated tissue to the
original grid and subtract it from the original IQ:

\[
  F(X)=X-R_{\tilde M(X)}\!\left[P_1(A_{\tilde M(X)}X)\right].
\]

The selected **primary arm** uses the whole field for the SVD fit. A diagnostic
ablation chooses the top quarter of pixels by median aligned IQ magnitude,
fits the temporal mode there, then applies that mode to every pixel. The
quarter is a fixed operating choice, not an anatomical mask. There is no
flow-off calibration, supplied tissue ROI, source truth, known displacement,
velocity, or pump setting in the fitting API. The one-frame smoothing width,
rank one, 25% diagnostic support, image-shift bound, phase step and 64-frame
offline ensemble are still **parameters**. The original phase estimate is
retained after geometric smoothing, so this does not establish a locally
phase-exact IQ warp.

Previously seen 5201/5202/8301/8302 cases were used to diagnose the idea.
The bright fit barely changed the result; the simpler smoothed full-field
version was designated the primary arm **before** new development 6301/6302
and held-out 9401/9402 cases. The [seen archive](../runs/axis1_stable_tissue_v1/seen_diagnostic.json),
[development archive](../runs/axis1_stable_tissue_v1/development.json) and
[freeze manifest](../runs/axis1_stable_tissue_v1/freeze.json) record that order.
Development failed the joint criterion; the locked candidate was nevertheless
run on the fresh bank to test generalization.

## Held-out results

Medians over six moving cases (translation, smooth deformation, and unmodeled
deformation with correlated tissue-texture turnover; two seeds each). Lower
is better for both measures. **Added-blood IQ error** is the relative complex
vector error of \(F(T+N+B)-F(T+N)\) versus known blood on the same cropped
support; it includes phase, amplitude and refit artifacts and can exceed one.
**Background power remaining** is the power of independently filtered
\(T+N\) divided by input tissue power, with lower values meaning better
tissue suppression. Neither metric is a clinical sensitivity estimate.

| Method | Added-blood IQ error ↓ | Background power remaining ↓ |
|:--|--:|--:|
| Ordinary rank-one SVD on native IQ | **0.589** | 0.597 |
| Registration then SVD, return filtered IQ | 1.003 | **0.042** |
| Previous native-grid tissue subtraction | 0.875 | 0.094 |
| **Smoothed motion, full-field native subtraction** | **0.833** | **0.090** |
| Smoothed motion, bright-quarter training | 0.833 | 0.089 |

Smoothing improved blood error against the previous native-grid method in all
six paired moving cases, but **all six still failed** the predeclared blood
gate of ordinary SVD plus 0.05. The smoothed full-field arm met the tissue
power limit of 0.10 in four of six moving cases. In the two turnover scenes it
left 0.227 and 0.230; one of those had a B-mode registration fallback. The
two stationary cases met the blood margin. Thus the combined gate failed in
development and held-out. Bright-quarter training did not materially change
either median and adds an operating choice, so it is not promoted.

## What caused the remaining error

Over the moving held-out cases, the median difference between motion maps
fitted to \(T+N+B\) and \(T+N\) fell from **0.067 to 0.041 pixels RMS** after
smoothing. When the background is *diagnostically* refitted using the
**mixture's same motion map**, the primary arm's median blood error is 0.699;
the frozen-operator blood error is 0.698. With independent background motion,
the added-blood error is 0.833. The median norm of the background-output
change caused by independent motion fitting is 0.347 in units of input blood
IQ norm. These vector errors do not add as scalar fractions. The same-motion
calculation uses a source-separated counterfactual **only in the evaluator**;
it is not part of the proposed filter. It shows that motion-refit sensitivity
persists, while the ~0.70 frozen blood error demonstrates substantial removal
of blood even with the motion and projector held fixed.

Bright-only training barely alters this projection because blood and tissue
share a slow-time direction in the short simulated ensembles. Restricting
where the rank-one mode is *learned* does not prevent that mode from removing
the overlapping blood component when it is *applied*. An exactly identical
blood and tissue temporal waveform cannot be separated by a temporal
projection from IQ alone. This simulation does not prove that all real
slow-flow waveforms are identical or that all other approaches fail.

The 64-frame, three-view 16 × 24 fields, four-pixel crop, last 32 scored
frames, flow speeds and turnover surrogate are specified in the
[configuration](../configs/axis1_stable_tissue_v1.json). Turnover is a
correlated complex-texture surrogate, **not** elevational acoustic motion.
The background counterfactual is independent reprocessing of known simulated
tissue plus noise; flow-on minus flow-off in a physical phantom is not exact
pixelwise blood IQ truth. The existing Resolve dataset analysis has no local
raw-IQ replay or flow-off source separation for this numerical test.

The candidate is therefore **not a robust calibration-free solution**. A
useful next test would need additional observable information that distinguishes
the overlapping signals, such as independently acquired tissue reference,
longer motion-stable histories, or verified cross-view/angular structure, and
must disclose the resulting acquisition and operating choices. For a physical
check, collect matched non-contrast flow-on/off complex IQ with B-mode,
independent motion and flow-speed measurements, static and controlled moving
runs, PRF, pixel spacing and tube geometry. Predeclare the tube and outside
regions and score tissue suppression, flow-on increment and visible slow flow
separately. Do not claim pixelwise blood preservation from tube dimensions or
pump speed alone.

Reproduce from the [runner](../examples/axis1_stable_tissue_study.py) and
[implementation](../src/btbench/stable_tissue.py) in a *new* output directory:
`OPENBLAS_NUM_THREADS=1 PYTHONPATH=src python examples/axis1_stable_tissue_study.py --split seen_diagnostic --out NEW`, then `--split development --freeze --out NEW`, then
`--split heldout --out NEW`. The committed
[held-out archive](../runs/axis1_stable_tissue_v1/heldout.json) is protected
against scoring the same path again. The [verifier](../scripts/verify_axis1_stable_tissue.py)
checks source leases, split separation, result coverage and numerical replay.
