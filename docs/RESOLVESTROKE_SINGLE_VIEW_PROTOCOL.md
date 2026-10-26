# Resolve Stroke single-view CEUS tube-mask diagnostic

26 October 2026. This is a **separate exploratory physical-data study**,
fixed before inspecting its filter scores. The source is the CC BY 4.0
[Resolve Stroke flow phantom](https://github.com/open-h/Open-H-RF/blob/main/datasets/resolvestroke/phantom_flow/README.md),
Hugging Face revision `05234aa`, HDF5 SHA-256
`40908ca699f3b88d9cdcdc4cabf9f1e07d481baec901b7044132e42c41f768f5`.
The exact [protocol](../configs/resolvestroke_single_view.json), [complex-IQ
Zea reconstruction](../configs/resolvestroke_complex.yaml), and [analysis
script](../scripts/resolvestroke_single_view.py) are hash-bound. The data are
not copied into this repository. Retain Resolve Stroke/Open-H-RF attribution
when reusing results.

## Relationship to phase 10

The released 20,000 frames are five 1-second intervals from **one**
microbubble-bolus acquisition. Every frame has a single diverging-wave
transmit. There are no independent transmit angles, controlled rigid/nonrigid
tissue motion, matched tissue-only capture or source-resolved blood IQ. The
published `mvi` and radial-velocity maps come from the device processing of
the complete exam; they are not independent flow truth. The `tube_mask`,
however, is known phantom geometry and labels the nominal 4 mm and 2 mm
tubes. The frozen three-view phase-10 mode selector and paired blood/moving
tissue gates therefore do **not** apply. No virtual transmit views are made
from the receive aperture or from repeated frames.

## Fixed analysis

Reverse the recorded axial hardware TGC and run Zea 0.1.6 delay-and-sum
beamforming to **complex IQ**, before envelope detection or log compression.
Use the reference Cartesian coordinate grid's plane nearest x = −5 mm, with
|y| ≤23 mm and 62–96 mm depth. This fixed plane intersects both tubes; the
mask is used to define that anatomical plane and the evaluation ROIs only.
The filter never sees tube labels.

For each of the five labeled 4,000-frame clips, use fifteen disjoint
256-frame ensembles starting at relative frames 80, 336, …, 3664. This covers
3,840 frames per clip and excludes 80 at each edge. Never process across the
unrecorded gaps. Match the input, grid, and slow-time support for:

- Unfiltered complex beamformed IQ;
- Global spatial SVD, discard exactly one singular mode per ensemble;
- Third-order zero-phase 100 Hz Butterworth slow-time high-pass, the conventional wall-filter control.

Integrate linear IQ power separately for each method and ensemble, then
average the fifteen power maps for each clip. For each tube, compare its
ground-truth mask pixels against background pixels at the **same tube depth**
(±2 mm) and at least two 0.6 mm pixels from either tube. Report area under
the ROC curve, tube-pixel detection at the outside-pixel 99th percentile,
mean tube/outside power contrast in dB, and outside residual power relative
to unfiltered IQ. The 4 mm and 2 mm labels are reported separately. Show
all five clips, the baseline-relative tube response, and the per-window
variation without treating correlated windows as independent replicates.

These metrics test *anatomical localization* of CEUS/Doppler signal against
known tube geometry. A bright tube pixel can still contain clutter and a dark
one can reflect lack of bubbles; neither is a direct estimate of blood-signal
retention. The `baseline_minus1s` clip precedes microbubble arrival, but its
name alone does not establish a pump-off or tissue-only control. No threshold,
rank, crop, or preferred method is chosen from the five scored clips.

The result should be published even if filtering fails. It is a single
acquisition diagnostic, not the phase-10 physical validation or a fresh
experimental holdout.
