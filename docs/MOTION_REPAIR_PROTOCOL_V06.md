# Motion compensation repair: protocol v06

Development starts on 23 October 2026 from merged PR #6. Its seen evaluation
seeds are historical diagnostics, not a new held-out bank.

This study repairs registration and adds a data-only rejection/fallback rule,
then repeats the common-view comparison. Development/training seeds are
**601, 787, 919**. The reserved evaluation seeds are **1601, 1907, 2203, 2503,
2801**. None of these evaluation draws may be inspected until the numerical
code, method grid, acceptance criteria and configuration are frozen.

The frozen configuration is `configs/motion_repair_v06.json`: eight original
shared-bank conditions, three training seeds and five new held-out seeds,
three angular views, 8×16 tile, 16-frame packets, 64-frame history, and
score frames [64,128). One matched fixed setting for each of uniform SVD,
weighted SVD, Axis 1, Axis 2, Axis 3 loaded STAP, Baranger SSM, and Lee AER
in temporal and spatial forms is applied to each of none, gated rigid, gated
smooth, and known warp (evaluation-only). Each motion arm also has a
roundtrip control. No per-condition or evaluation-seed tuning is permitted.
The usual global training-only operating-point selector remains available,
but with one setting per family/arm there is no tuning freedom in this grid.

The scene generator remains byte-for-byte governed by the v05 boundary rule.
Filtering alignment uses continuous bilinear zero extension, so a tiny shift
does not erase a whole edge pixel. Registration aligns calibration-derived
complex view templates to each native frame. Per-view complex nuisance gains
are projected out on fitting pixels only; they never correct filtered IQ or
ground truth. A fixed tissue ROI excludes its outer row and column; pixels
where `(z+2*x)%3 == 0` are held out spatially for each frame. The remaining
pixels fit either 2 rigid translation parameters or a four-term smooth basis
per coordinate. Robust least squares uses the exact projected Jacobian,
bounded coefficients and a smooth-shape ridge. A two-start coarse rigid fit
initializes smooth fitting. A smooth fit may replace the rigid one only when
it improves independent spatial validation by at least 0.003 absolute and
10% relative. The exact options, including displacement, conditioning,
support, Jacobian, inverse-consistency and validation gates, are serialized
in the config. Static-consistent, unobservable, unsupported, implausible,
low-improvement or nonconverged fits fall back to rigid or exact zero.
Ungated variants exist only as local diagnostic ablations, not in the common
held-out method grid. All acceptance uses observed IQ and calibration,
never prescribed warp or source truth. A fitted motion trajectory uses no
future frame; the registration reference is from negative-time calibration.

Report every condition and seed, paired actual complex-IQ blood error,
refitted background power, incremental refit responses, motion RMSE against
prescribed coordinates (evaluation-only), gate counts, and roundtrip error.
The software gate for phantom work requires: exact zero on static/mismatched/
unobservable challenge controls; no accepted displacement with invalid map;
lower mean motion RMSE than v05 on both translation and smooth conditions;
and no material degradation (>5% mean paired actual blood error) from the
no-compensation setting for either SVD or Axis 1 on the complete held-out
bank. Failure of any criterion is reported plainly and postpones phantom
validation rather than adjusting the frozen rule to the held-out bank.
This is a statistical-IQ software gate, not an experimental/RF claim.

The user authorized validation and pushes directly to `main` without another
branch, PR or external approval. Historical archives remain unchanged.
