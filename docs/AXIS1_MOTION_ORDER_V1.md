# Axis 1: B-mode motion before or after clutter filtering

26 October 2026. This is an **offline source-separated IQ simulation**, not a
new Resolve stroke flow phantom acquisition or clinical validation. The
development bank was inspected and the identity-versus-translation gate was
corrected there. Source and config hashes were sealed in
[`freeze.json`](../runs/axis1_motion_order_v1/freeze.json) before opening the
independent 7201/7202 held-out seeds. The held-out archive was run once.

## Question and result

Does a simple, data-only envelope motion fit before SVD preserve slow flow
*while* suppressing moving tissue? In this bank, **no single arm does both**.
Registration before SVD markedly reduces simulated background, but also makes
the refitted added-blood IQ less faithful. The automatic gate correctly leaves
both stationary cases at identity. The richer affine and smooth fits are never
selected on held-out scenes. A global per-view phase correction barely changes
the scores. These are paired, small-sample findings, not a universal ranking.

Held-out medians for the **six moving cases** (three imposed 2-D motion classes,
two seeds). Lower is better in both columns. Added-blood IQ error is zero only
for perfect complex-signal preservation; values above one are possible.
Background remaining is output power for independently refitted tissue+noise
divided by the unfiltered tissue power, on the same cropped region.

| Arm | Added-blood IQ error ↓ | Background power remaining ↓ |
|:--|--:|--:|
| Ordinary SVD, raw IQ, rank 1 | **0.584** | 0.500 |
| SVD then geometric registration, rank 1 | 0.686 | 0.351 |
| Geometric registration then SVD, rank 1 | 0.873 | 0.063 |
| Geometry and global IQ phase then SVD, rank 1 | 0.870 | **0.063** |
| Current Axis 1, calibrated, bounded history | 0.811 | 0.212 |
| Fixed overlapping windows, rank 1 | 0.689 | 0.345 |
| Automatic motion-limited windows, rank 1 | 0.738 | 0.309 |
| Ordinary SVD, raw IQ, rank 2 | 0.753 | 0.278 |
| Geometry and phase then SVD, rank 2 | 1.089 | 0.043 |

The registered rank-1 arm improves tissue suppression in all six moving
pairs and worsens added-blood error in all six, relative to raw rank-1 SVD.
Frozen-operator blood IQ error is 0.744 for the registered phase arm versus
0.581 for raw SVD; the paired refit error (0.870 versus 0.584) additionally
captures motion and SVD changing when blood is added. The existing calibrated
Axis 1 has a much smaller *frozen* blood error (~0.249), but its **refitted**
added-blood error is ~0.811. Those metrics answer different questions.

For the two **no-motion cases**, ordinary rank-1 SVD and geometric-then-SVD
have median added-blood errors 0.768 and 0.768, with effectively zero
background remaining. The identity gate selected zero geometric motion in all
63 nonreference frames of each case.

## Exact arms and information available to them

The scene is `shared_scene`, with 64 acquired frames, three coherent IQ views,
16 × 24 spatial pixels, relative slow flow of 1.2 mm/s, tissue motion, and
source-separated blood/tissue/noise for **evaluation only**. Conditions are no
imposed image motion, translation, smooth deformation, and deformation outside
the fitted basis. Only frames 32–63 and a common four-pixel spatial crop are
scored; each arm is reduced to the coherent mean of the three views. There are
two development and two disjoint held-out seeds per condition. These images
are specified statistical IQ fields with imposed spatial maps. They are not
RF/channel acoustics, real red cells, or an out-of-plane motion test.

1. **B-mode-derived geometry.** Form the log-compressed mean of the three
   observed IQ magnitudes. Select the brighter half of its interior support;
   split that support in checkerboard fitting/check pixels. Fit framewise
   subpixel 2-D translation to frame zero, with a maximum of 2.5 pixels.
   Compare identity, translation, optional affine and optional smooth
   sine-basis deformation by held-out image-residual BIC. Affine and smooth
   coefficients are bounded at ±0.5 pixel in normalized coordinates. A check
   normalized correlation below 0.5 triggers identity fallback. Adjacent
   pixels are correlated, so the BIC is a pragmatic automatic gate, not a
   calibrated significance test. No true motion, flow mask, calibration
   acquisition or hand-drawn ROI is used in these new fits.
2. **Geometry applied to IQ.** Bilinearly interpolate real and imaginary IQ
   at the *same* fitted coordinates in every view. Restore the result to the
   original grid using an iterated inverse map with zero extension. Fit motion
   on the observed mixture; replay its frozen map on blood/tissue/noise for
   source attribution. A separate ablation multiplies aligned IQ by the
   conjugate of its bright-tissue complex cross-phase relative to frame zero,
   per view, and undoes that phase at restoration. This is a **single global
   phase per view and frame**, not a validated spatially varying axial carrier
   correction. Refit motion and phase independently on background for the
   added-blood counterfactual.
3. **Filtering order.** `svd_filter` fits an uncentered joint-view Casorati
   operator and removes the first one or two modes. “Filter before” fits on
   native IQ and then performs the same frozen geometry round trip. “Register
   before” fits after geometry (and optional phase), then restores. The
   Baranger arms use the repository's spatial-similarity square-fit rank band
   on the same joint-view matrix, before or after registration. They do not
   use Zhang's Maresca cutoff.
4. **Current Axis 1.** `SharedContext` uses a **separate 32-frame calibration
   sequence and supplied tissue training ROI** to fit phase/view gains and a
   complex-reference smooth motion field. `fit_shared(axis1)` applies rank-1
   spatial moment filtering with 32-frame bounded history, packet size 16 and
   forgetting factor 0.75. It sees more prior information than the new arm;
   its score is a reference comparator, not a no-calibration contest.
5. **Motion-adaptive windows.** A distinct automatic-ROI, Zhang-inspired
   control applies SVD directly on native IQ in overlapping windows with
   an estimated translation path budget of one wavelength, maximum 48 frames,
   minimum 32 and hop eight. Fixed windows use the same overlap/add and maximum.
   Each window uses a frozen rank-1 operator or the existing Maresca curvature
   rule with a **fixed** offset of 26. The actual [Zhang et al. 2022](https://doi.org/10.1088/1361-6560/ac7430)
   method uses motion in a manually selected myocardium ROI, windows adapted
   to its maximal cumulative displacement, and an operator-adjusted offset
   per ensemble. The present control **is not a faithful reproduction** of
   those manual choices or the separate calibrated FMBV attenuation step.

## Ablations, failure checks and limits

- **Registration gate:** the six moving held-out scenes select translation in
  48–51 of 63 nonreference frames, identity in the others. There were zero
  envelope-NCC failure fallbacks. No affine or smooth field was selected; the
  richer models therefore have identical held-out results to translation.
  “Not accepted” frames from existing Axis 1 include identity decisions and
  must not be read as the same failure count.
- **IQ phase:** the maximum paired change in added-blood error from the global
  phase step was below 0.004. This does not rule out a locally phase-consistent
  warp on a higher-resolution IQ acquisition.
- **Rank sensitivity:** rank two lowers remaining tissue relative to rank one
  for ordinary SVD (0.278 versus 0.500) and registered SVD (0.043 versus
  0.063), while increasing added-blood error (0.753 versus 0.584, and 1.089
  versus 0.870). Rank remains a consequential choice.
- **Interpolation:** the no-filter geometry/phase round trip changes complex
  IQ by 0.176–0.201 relative norm in moving held-out cases, on the scored
  region. Some blood penalty belongs to resampling itself. A native-grid
  output and a registration-first output must therefore both report this
  control.
- **Window sensitivity:** on moving cases the adaptive windows all hit their
  32-frame minimum, whereas fixed windows are 32–48 frames. Their modest
  tissue improvement comes with more blood error; the minimum may violate the
  wavelength budget. This bank does not isolate the published unrestricted
  window rule. With offset 26, both Maresca controls essentially zero all
  blood and background; that operator choice is not portable to these short
  simulated ensembles.
- **Baranger:** the automatic band often retains essentially no slow blood in
  the moving conditions (added-blood error near one and frozen-blood error near
  one), despite low tissue output. This is specific to these small simulated
  ensembles. It does not contradict Baranger's clinical utility on other data.
- **Failure coverage:** no measured B-mode/IQ pair, elevational motion,
  decorrelation stress, operator-selected Maresca offset, or independent
  ground-truth source decomposition on the Resolve phantom is in this archive.
  Tube position and pump speed alone cannot give pixelwise blood/tissue IQ
  errors. The new arm uses zero external calibration but still has a motion
  limit, model gate, rank/window settings and phase ablation. It is not
  parameter-free.

The observed trade-off does **not** justify advancing registration-first
rank-1 SVD as the practical Axis 1 solution. A next candidate should reduce
IQ round-trip loss and keep tissue rejection from turning into added-blood
loss. Any revision must use a new development protocol and a fresh held-out
bank; these held-out seeds are now seen.

## Reproduction and archives

`OPENBLAS_NUM_THREADS=1 PYTHONPATH=src python examples/axis1_motion_order_study.py --split development --freeze`
then `--split heldout`, after the sealed sources match. The exact
[configuration](../configs/axis1_motion_order_v1.json),
[source](../src/btbench/motion_first.py),
[runner](../examples/axis1_motion_order_study.py),
[development](../runs/axis1_motion_order_v1/development.json) and
[held-out archive](../runs/axis1_motion_order_v1/heldout.json) are retained.
The archive records per-case metrics, motion-model counts, failure definitions,
phase coherence, window lengths/discards, runtime and the full source lease.
