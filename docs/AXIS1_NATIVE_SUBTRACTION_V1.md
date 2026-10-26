# Axis 1 native-grid tissue subtraction: frozen test

26 October 2026. **Decision: the new candidate fails the predeclared joint
blood-preservation and tissue-suppression criterion.** This is a source-labeled
IQ simulation, not a physical Resolve phantom result. No parameter was changed
after opening the fresh 8301/8302 held-out seeds.

## Why this candidate

The previous [motion-order study](AXIS1_MOTION_ORDER_V1.md) found stronger
tissue suppression after B-mode-derived registration but worse slow-blood IQ.
Its no-filter complex-IQ alignment/restoration round trip changed the moving
signal. Let \(A\) denote the frozen IQ alignment, \(R\) its numerical inverse
resampling, and \(P_T\) the SVD tissue projector fitted to the *observed
mixture*. The previous output was \(R[(I-P_T)A(X)]\). This study changes only
the output construction:

\[
    Y_{\mathrm{native}}=X-R[P_TA(X)].
\]

The estimated tissue is mapped back and subtracted from the **original IQ**.
The blood residual itself does not take an inverse interpolation. This is still
rank-1 SVD with the previous automatic translation fit and global per-view
IQ phase step. It has no tissue-only calibration acquisition, supplied tissue
ROI, extra fitted rank or new output-tuning coefficient. It is not a claim of
novelty for registering before filtering, which has prior ultrasound work.

For the same fitted motion and projector, the difference is exactly
\(Y_{\mathrm{native}}-Y_{\mathrm{return}}=X-RA(X)\). We verified this identity
on mixture and isolated blood to below 10⁻¹⁶ relative norm. It isolates the
round-trip term **without** asserting that the tissue projector itself preserves
blood. The frozen replay remains linear and source-additive.

## Split and decision

The [seen diagnostic](../runs/axis1_native_subtraction_v1/seen_diagnostic.json)
reuses the previously inspected 4101/4102/7201/7202 scenes for mechanism
analysis only. Two new seeds, 5201/5202, form the
[development archive](../runs/axis1_native_subtraction_v1/development.json).
Source/config hashes, those two archives and the decision thresholds were
sealed in the [freeze manifest](../runs/axis1_native_subtraction_v1/freeze.json)
before the new 8301/8302 seeds were scored once in the
[held-out archive](../runs/axis1_native_subtraction_v1/heldout.json).

The fresh bank has four conditions per seed: stationary slow flow, translating
very slow flow, smooth deformation with still slower flow, and an unmodeled
in-plane deformation with added correlated tissue-texture turnover. The
turnover is an explicitly specified complex-speckle surrogate; it is **not** a
3-D/elevational acoustic or RBC model. All arms use the same 64-frame,
three-view 16 × 24 IQ, last 32 frames, and four-pixel common crop.

The frozen gate required, in *every* six moving held-out cases, added-blood IQ
error at most raw SVD plus 0.05, background power remaining at most 0.10, and
zero registration failures; the two stationary cases had the same blood-error
margin. It failed in development and independently failed in held-out.

## Held-out results

Medians over the six moving cases; lower is better in both columns. The
added-blood error compares the independent refit response
\(F(T+N+B)-F(T+N)\) with known blood IQ. Background remaining is the power of
\(F(T+N)\) divided by input tissue power on identical support.

| Method | Added-blood IQ error ↓ | Background power remaining ↓ |
|:--|--:|--:|
| Ordinary raw rank-1 SVD | **0.533** | 0.598 |
| Registered SVD, returned filtered IQ | 1.962 | **0.076** |
| Registered tissue estimate subtracted on native IQ | 1.645 | 0.126 |

Native subtraction lowers added-blood error relative to registered return in
five of six moving cases, but remains much worse than raw SVD on blood
preservation in all six. It also leaves more tissue than registered return;
the turnover cases leave 0.224 and 0.264 of input tissue power. Two moving
cases each had one envelope-correlation registration fallback; the stationary
cases selected identity and showed no meaningful difference between native
subtraction and registered return.

The blood-only IQ round trip has median error 0.183 in moving held-out cases.
The **frozen fitted tissue projector**, returned to the native grid, removes
known blood with median relative norm 0.696. These quantities are vector
errors and must not be added as independent power fractions. The frozen blood
error for native subtraction is 0.696, while its refitted added-blood error is
1.645. This gap shows that source attribution with a fixed fitted operator
does not capture changes in registration and SVD when blood is added to
background. It does not by itself identify which refit change dominates.

## Practical interpretation and measured-data handoff

Avoiding inverse interpolation of the residual addresses one real error term,
but **it does not solve the overlap between slow blood and the learned tissue
subspace**. The present candidate should not be promoted as a robust,
calibration-free clutter filter or retuned on 8301/8302. The next numerical
question is how to estimate tissue without learning away slow blood or making
the refit unstable, while keeping a small, fixed set of operating choices.

To evaluate future candidates physically, acquire **non-contrast** complex IQ
and matching B-mode (or a validated structural envelope) for flow-on and
flow-off repeats, static and controlled in-plane motion, with independent
flow speed, PRF/frequency, pixel spacings and tube geometry. Use flow-off and
independent motion measurements for **evaluation**, never as calibration
inputs to a claimed no-calibration method. Predeclare tube and outside-tube
regions, motion-quality failures, filter settings and comparison metrics.
Report background power in flow-off, the flow-on increment and slow-flow
visibility separately. Flow-on minus flow-off is **not** exact pixelwise blood
IQ truth when speckle and motion differ. The existing Resolve stroke flow
analysis lacks a local raw-IQ replay here and cannot independently establish
this source-preservation result. No physical data are claimed in this study.

Reproduce with
`OPENBLAS_NUM_THREADS=1 PYTHONPATH=src python examples/axis1_native_subtraction_study.py --split seen_diagnostic`,
then `--split development --freeze`, then `--split heldout` on a clean new
output directory. The [implementation](../src/btbench/native_subtraction.py),
[configuration](../configs/axis1_native_subtraction_v1.json) and
[runner](../examples/axis1_native_subtraction_study.py) identify every fitted
operator and source hash. The archived held-out file is not rerunnable in place
because its path is guarded against scoring twice.
