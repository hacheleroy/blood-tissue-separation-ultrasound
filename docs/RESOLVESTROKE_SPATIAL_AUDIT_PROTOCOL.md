# Resolve Stroke anatomical and temporal robustness audit

26 October 2026. This is an **exploratory extension** of the completed
[single-plane result](RESOLVESTROKE_SINGLE_VIEW_RESULT.md), frozen before
examining new cross-plane filter scores. It reuses the exact 5.25 GB source
HDF5, Zea 0.1.6 complex-IQ delay-and-sum reconstruction with recorded TGC
reversed, 256-frame complex ensembles, and the already fixed unfiltered,
SVD-cut-1 and zero-phase 100 Hz high-pass methods. The
[configuration](../configs/resolvestroke_spatial_audit.json) binds the
[analysis script](../scripts/resolvestroke_spatial_audit.py), beamformer,
predecessor method source and prior 75-window report by SHA-256. The
ground-truth tube mask selects structural planes and supplies evaluation
labels only. No score can adjust a plane, window, setting or variant.

## Fixed extension

Use five x indices **25, 52, 77, 105, 132** on the supplied 0.6 mm grid,
approximately −36, −19.8, −4.8, +12 and +28.2 mm. Geometry-only inspection
shows both labeled tubes in each plane. Keep the earlier 57×77 y/z crop,
outside guard and separate 4 mm / 2 mm labels. For each of the five clips,
score the windows at prior schedule indices **0, 7, 14**, covering early,
middle and late 256-frame intervals. Beamform all five planes from the same
channel blocks, without combining them into virtual transmit views. Average
linear power over the three disjoint windows separately for each plane and
clip. Recalculate AUC, tube-pixel fraction above the depth-matched outside
99th percentile, tube/outside power contrast and outside power relative to
raw IQ; retain both tubes and all five planes without picking a best plane.

On each clip-mean map, additionally replace a tube's positive mask by its
one-pixel cross dilation or one-pixel translation in each of the four grid
directions. Keep the **original** depth-matched, >2-pixel-guarded outside
mask fixed for this sensitivity check. Report all five perturbations and
the original; the lowest AUC is a declared descriptive summary, not a
selection rule. A one-pixel shift is 0.6 mm and especially consequential
for the nine-pixel 2 mm tube. Do not interpret the shifted label as a new
physical ground truth.

Use the prior, SHA-bound 75-window report to plot unfiltered and filtered
linear tube power through all fifteen recorded windows in each clip at
x = −4.8 mm. Keep the unrecorded interclip gaps visible; do not interpolate
or treat adjacent windows as independent acquisitions. Check the designated
−1 s clip as a negative contrast-bolus control, without calling it a
pump-off or independently measured tissue-only acquisition. The five planes,
windows and variants are correlated observations from one bolus exam;
report ranges and failures, not inferential confidence intervals.

All spatial and temporal outcomes measure CEUS signal localization relative
to known tube geometry. They cannot estimate exact blood-IQ preservation,
separate tissue-only leakage, velocity truth, or the three-view phase-10
motion candidate. A low score on an off-center plane is a result, not a
reason to move the plane after inspection.
