# Axis 1: matched-tissue temporal projection limit

26 October 2026. **The dominant loss on this software bank is intrinsic
slow-time overlap with the tissue direction, not primarily contamination of
the estimated direction by blood.** Longer histories reduce that overlap but
do not meet the strong blood-preservation target for the slow-flow cases.
With correlated tissue-texture turnover, even known imposed geometry and a
true-tissue rank-one or rank-two temporal mode fail the background target.
This is a bound **within the tested fixed-mode native-grid subtraction rays**,
not an impossibility result for all IQ filters or an experimental finding.

## Diagnostic and information boundary

Generate one 256-frame, three-view, 16 × 24 source-separated IQ scene at
1 kHz for each seed, motion condition, and flow speed. The speed is the
**relative** axial flow component added to the tissue motion in the specified
generator. Fixed speeds are 0.4, 0.7, 1.2, and 2.4 mm/s; conditions are no
image motion, translating image motion, and deformation outside a translation
model plus correlated complex tissue-texture turnover. The latter is a
specified surrogate, not an elevational acoustic model. The source generator
has a fixed, synthetic 256-frame motion path for a given scene. History
lengths 64, 128, and 256 take nested intervals **ending at the same time**;
all metrics score their identical last 64 frames and four-pixel crop. Tissue
and noise realizations are also paired across the four speeds. The archive
includes hashes of the common scored blood, tissue, and noise targets.

The evaluator obtains the exact *imposed* displacement field and aligns IQ
with the existing bilinear complex-IQ sampler. Inverse resampling is still
numerical and imperfect at boundaries. It fits one or two right singular
vectors from **isolated true tissue** on each aligned history. This is an
evaluation-only oracle; neither the motion nor true tissue is available to a
practical filter. A mode fitted from the observed mixture using the same
known geometry is a separate diagnostic. Ordinary raw rank-one SVD is a
data-only comparator on the identical history and scored target.

For a fixed tissue-derived projector \(P_k\) and known geometry \(A,R\),
the family is

\[
    F_{k,\alpha}(X)=X-\alpha R[P_kA(X)],\qquad 0\leq\alpha\leq1.
\]

\(\alpha\) is **only a diagnostic sweep of subtraction strength**, not a
new candidate parameter. The tissue constraint is background output power
\(\|F(T+N)\|^2/\|T\|^2\leq0.10\) on the common scored support. For each
fixed projector the evaluator solves the exact quadratic in \(\alpha\) and
records the *smallest* gain satisfying that constraint, if any. Blood IQ
error along this ray is \(\alpha\|RP_kA(B)\|/\|B\|\), so that smallest gain
minimizes blood error among feasible gains on this ray. A separate ambitious
benchmark requires blood IQ error at most 0.25. This declared tolerance is
a research target, not a clinical acceptance threshold; the archive also
records continuous errors and a comparison with raw SVD plus 0.05.

Previously seen seed 9401 was used for numerical diagnosis. Source/config
hashes and the [seen archive](../runs/axis1_oracle_limit_v1/seen_diagnostic.json)
and new [development archive](../runs/axis1_oracle_limit_v1/development.json)
were sealed in the [freeze](../runs/axis1_oracle_limit_v1/freeze.json) before
12701/12702 were scored once in the
[held-out archive](../runs/axis1_oracle_limit_v1/heldout.json). Development
seeds 10601/10602 are disjoint from both. No settings were changed after
development.

## Held-out result

Medians across the eight cases in each motion class (two seeds × four speeds).
All entries are at **full** tissue-mode subtraction \(\alpha=1\). Lower is
better for both columns. Blood error is the norm of complex IQ distortion
relative to true blood, so it is not a percentage of blood power. Tissue
power remaining is evaluated on true tissue plus noise. The rank-one mode
itself comes from true tissue, while its output is returned to the original
IQ grid.

| Motion class | History | Oracle blood IQ error ↓ | Background power remaining ↓ |
|:--|--:|--:|--:|
| Stationary | 64 | 0.853 | 0.001 |
| Stationary | 128 | 0.688 | 0.004 |
| Stationary | 256 | 0.533 | 0.011 |
| Translation | 64 | 0.606 | 0.083 |
| Translation | 128 | 0.493 | 0.070 |
| Translation | 256 | 0.405 | 0.063 |
| Unmodeled motion + turnover | 64 | 0.588 | 0.205 |
| Unmodeled motion + turnover | 128 | 0.481 | 0.222 |
| Unmodeled motion + turnover | 256 | 0.399 | 0.225 |

Across **all 16 moving held-out cases**, median rank-one blood error falls
from 0.595 to 0.488 to 0.400 as history increases, whereas median remaining
background is 0.145, 0.146, and 0.145. At 64 frames, eight translation cases
can attain the tissue-power limit by reducing subtraction strength, but **none
of the 16** can jointly attain it with blood IQ error ≤0.25. At 128 frames the
joint count is also zero. At 256 frames **two of 16** pass both targets, both
the fastest 2.4 mm/s translation cases. None of the 0.4 or 0.7 mm/s moving
cases pass at any tested length. In the turnover class, no rank-one gain in
[0,1] attains the 0.10 tissue limit. Rank two lowers the median turnover
background to 0.157 at 64 frames, but it still misses 0.10 and increases
blood distortion; its medians at 128 and 256 frames are 0.191 and 0.222.

The *known-motion, mixture-fitted* rank-one SVD and the *known-motion,
true-tissue-fitted* rank-one oracle have nearly identical added-blood errors
on this bank: median absolute paired differences in the 16 moving cases are
0.0002, 0.0004, and 0.0004 at 64, 128, and 256 frames (max difference <0.0013).
Thus, the previous idea that blood mainly *contaminated the learned tissue
direction* is not supported here. The tissue direction itself overlaps slow
blood, even when fitted without blood. Separately, the previous
[stable-motion study](AXIS1_STABLE_TISSUE_V1.md) showed sensitivity to
refitting motion on the observed mixture. These are two distinct failure
mechanisms.

## Practical reading

More temporal history helps blood preservation in these paired scenes, but
the 256-frame horizon is only 256 ms and does not test the very long,
motion-stable ensembles required by LE-SVDn. The history comparison changes
the available **past** frames while holding the last 64 scored frames and
physical motion path fixed. It does not claim that extending an acquisition
is free of additional motion or tissue decorrelation. The turnover surrogate
and two held-out seeds per condition limit generalization.

This experiment does **not** rule out a simpler useful filter based on other
observable structure, such as spatial or independently acquired angular
information, and does not establish an information-theoretic limit for all
algorithms. It does show why fitting a cleaner version of the same rank-one
temporal mode is unlikely to satisfy our strong slow-flow goal on these
scenes. The next candidate should justify its *additional discriminating
observable* and number of operating choices, then be tested on a fresh bank
and non-contrast measured flow-on/off IQ. The Resolve bolus dataset cannot
provide pixelwise blood IQ truth for this diagnosis.

The [runner](../examples/axis1_oracle_limit_study.py),
[configuration](../configs/axis1_oracle_limit_v1.json), and
[verifier](../scripts/verify_axis1_oracle_limit.py) reproduce and guard the
numerics. To rerun, use a **new output directory** and execute
`OPENBLAS_NUM_THREADS=1 PYTHONPATH=src python examples/axis1_oracle_limit_study.py --split seen_diagnostic --out NEW`, then `--split development --freeze --out NEW`, then `--split heldout --out NEW`.
