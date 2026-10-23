# Numerical verification

## PDF-audited implementation, v0.3.0 draft

23 October 2026: **134 tests pass** on Python 3.12. Numerical package hash:
`b249dda5158ec26f472e4ba36e09b5c7cadd79fc755751646a7b300ea0af160b`.

New tests cover PDF-derived threshold equations, coherent Song normalization,
Lee temporal/spatial covariance denominators, left/right mixture equivalence
and unequal frozen attribution, refitted injection responses, calibration-only
operating-point selection, Axis 2 ridge/gain controls and batched RF identity.
Experimental source-truth metrics are excluded. These implementation checks
are not validation of experimental blood/tissue separation.

The [PDF audit](PDF_METHOD_AUDIT.md) records unspecified numerical choices and
source limits. The [protocol](FOLLOWUP_PROTOCOL_V04.md) fixes the simulation grids
and public-data calibration/evaluation split. Publication of generated follow-up
figures and result archives remains pending explicit approval; those outputs
are excluded from this draft.

CI checks Python 3.10/3.12. Local package-content, installed-wheel and source
checks are recorded in [the package record](PACKAGING_V03.json). Direct tested
dependency versions are pinned; no full transitive-environment lock is claimed.

## Separation-fidelity milestone, v0.2.0

22 October 2026: **114 tests pass** locally on Python 3.12. Numerical package
fingerprint (sorted module names and bytes):
`b4f4d2323ac7c35a6b9c2001982ac1e5e9bdf56b9039e96cbe1fd847c00bdb8e`.
All five new execution archives contain this source identity. Configurations
were written before execution; generated arrays and detailed reproducible case
traces are excluded from Git. Full portable records are losslessly compressed
as `results.json.gz`.

| Execution | Cases | Method/settings outputs | Returned |
|---|---:|---:|---:|
| [Axes 1–2/source fidelity](../runs/separation_fidelity_v03/SUMMARY.md) | 30 | 1,170 | 1,170 |
| [Literature SVD/EVD](../runs/svd_fidelity_v03/SUMMARY.md) | 24 | 480 | 480 |
| [STAP fidelity](../runs/stap_fidelity_v03/SUMMARY.md) | 42 | 504 | 504 |
| [Waterloo RF pilot](../runs/openhrf_waterloo_pilot_v03/SUMMARY.md) | 1 short packet | 58 | 58 |
| [Default integration](../runs/integration_fidelity_v03/SUMMARY.md) | 1 | 55 | 55 |

Total: **2,267 returned method–case outputs**, including parameter variants and
truth-assisted searches. The 42 LCMV spectra are separate detector outputs.
Returning an output does not establish solver convergence, source
identification, acoustic validity or experimental separation accuracy.

The [archive manifest](../runs/fidelity_v03_manifest.json) records compressed and
decoded SHA256 values, counts, frozen configuration hashes and numerical source
identity. Verify with:

```bash
python scripts/verify_results.py
```

The separation-fidelity archive is stored in six numbered parts to keep upload
requests small. The manifest records every part's size and SHA256. The verifier
checks the parts and their concatenation against the original archive hash;
`--restore` also writes the original `results.json.gz` for plotting or other tools.
All numerical source, configurations, result bytes and original hashes are retained.

New numerical tests cover phase-sensitive blood fidelity, spectral shifts,
empty-output semantics, preservation-constrained/infeasible oracles, frozen
Axes 1–2/ICA linearity through motion/reset/pooled centering, ICA fallback,
transported-point phase for both signs, nonlinear injection interpretation,
embedded views, recorded event timing and RF-to-NPZ-to-STAP integration. A
separate pilot check retains clock origin across calibration/evaluation slices
and removes all truth fields from experimental scoring.

Maximum frozen component-closure errors: **2.47e-13** in the new separation
bank, **8.33e-15** in the literature bank and **8.27e-14** in STAP. These are
implementation identities. The 2-D PSF model adds transport and boundary mixing;
it does not establish realistic RBC/channel or fluid-structure physics.

ICA converges in **19/712** centered 300-iteration fits and **22/712** centered
1500-iteration fits. The default guarded arm explicitly falls back to local SVD
for 693 fits; low-confidence motion fallback is recorded separately. The
uncentered retained-output ablation converges in 63/712 fits. The identifiable
non-Gaussian source-recovery unit test remains successful. Additional iterations
alone do not resolve this short-packet, Gaussian-dominated benchmark.

Axis 2 exhibits large oblique-operator gain in this bank: maxima approximately
66.1 without memory and 64.4 with memory. This is a recorded failure mode;
no general performance improvement is claimed. Three STAP IQ-MSE oracle cases
still select empty output. Component metrics, invalid readouts and candidates
remain in the archives.

The RF pilot decodes 64 Waterloo CCA frames by range reads and benchmarks a
32-frame evaluation packet after disjoint calibration. Both views and their
166.667-microsecond offset are retained. Provider B-mode structural correlation
is 0.700319. There is no experimental source truth or phase-reference calibration.
See [experimental protocol](EXPERIMENTAL_PILOT.md) and [dataset evidence](DATASETS.md).

CI is configured for Python 3.10/3.12 numerical tests, archive integrity, source
and wheel contents, installed-wheel execution outside the source tree and tests
from the source distribution. Local execution uses Python 3.12. The v0.2.0 wheel was installed in a separate
target and imported outside the source tree; its numerical fingerprint matches
the executed studies. All 114 tests also pass from the unpacked source
distribution. [Package verification record](PACKAGING_V02.json). Tested direct numerical/I/O
versions are pinned; no full transitive-environment lock is claimed.


## Historical literature / Axis 3 revision before v0.2.0

Date: 22 October 2026. **85 numerical tests pass.** Numerical implementation
hash (SHA-256 over sorted package module names and bytes):
`f6cfb106623e3631b51c1f98758ec6a13267ac13d17875e1fbbbb40e76060d3c`.
The verified study records and integration run all contain this same hash.
The current wheel builds and installs into a clean target directory. Running
the installed CLI outside the source tree exposes both new study commands;
the installed Butterworth operator and literature/STAP modules import correctly.

| Execution | Cases | Comparison outputs | Result |
|---|---:|---:|---|
| `runs/integration_v02_verified` | 1 | 55 | 55 returned outputs, including the legacy Axes 1/2/AER/ICA/PCP arms |
| `runs/svd_literature_v02_verified` | 24 | 480 | All 20 arms returned on every case |
| `runs/stap_v02_verified` | 42 | 504 | All 12 comparison outputs returned on every case |
| STAP LCMV detector | 42 | 42 spectra | Candidate-grid scores and rejected constraints retained |

The default configuration has 52 data-only entries; the demo adds three
same-case oracle searches. The STAP oracle searches 45 lower/upper bands;
its 42 selected outputs are included in the 504 above. Candidate evaluations
are not counted as independent study cases.

New checks compare the Baranger sum-table objective with an independent,
literal binary-template/corr2 search, including boundaries beyond rank 64.
They also cover complex EVD/direct-SVD equivalence, randomized range exactness,
spatial-subsampling component rejection/preservation, eigenvalue versus singular
value thresholding, Gaussian weighting before covariance estimation, frozen
component linearity, analytic Butterworth cutoff gain for both Doppler signs,
causal initialization, padding limits and acquisition-PRF consistency.

STAP checks cover absolute snapshot timing, complex angular gains, sequential
view offsets, both IQ phase signs, single-column projectors, independent
calibration recovery, guard exclusion, insufficient-training fallback, frozen
view operators, LCMV known-target detection and unidentifiable constraints.
Across the 42-case study, the largest relative frozen-mixture reconstruction
error is **8.27e-14**. The largest accepted LCMV constraint-response error is
**1.47e-12**. These establish numerical identities, not model fidelity.

Complete execution means an output was returned. Invalid velocity/WSS estimates,
empty-output oracle selections and iterative convergence flags remain explicit.
Song's reconstructed tail fits converged in all 576 local fits of the 24-case
bank; this does not establish parity with the unavailable original algorithm.
The ICA integration arm still contains iteration-limit fits on this Gaussian
surrogate (59 of 60 patch/packet fits) and must not be presented as a converged
ICA result. Three STAP IQ-oracle cases select an empty band; their invalid
velocity/WSS outputs remain in the records.

Configs are saved before execution; input and calibration hashes, numerical
environment, fitted cutoffs, fallbacks, oracle objectives and candidate curves
are preserved. Small CPU timings are operational diagnostics. They are not a
controlled hardware throughput comparison. Generated IQ arrays are excluded
from Git. Detailed SVD case traces are reproducible locally; the tracked portable
aggregate preserves all global metrics and fitted-method diagnostics.

See [source fidelity](SVD_LITERATURE_AUDIT.md), [the STAP protocol](STAP_BENCHMARK.md),
and [the first results](RESULTS_FIRST_BENCHMARK.md). Original Song threshold
parity remains unverified. IEEE document 8466597 is identified as Maresca et al.
(2018), DOI 10.1109/TUFFC.2018.2870083; its full-text algorithm check is pending
because publisher and HAL PDF retrieval were blocked.

## Historical initial package (version 0.1.0)

Date: 22 October 2026. CPU, Python 3.12, NumPy/SciPy, one BLAS thread per run.
Pinned tested direct dependencies are in `requirements-tested.txt`.

## Automated numerical tests

54 tests pass. Coverage includes:

- Complex Casorati reshape/restore and uncentered/centered SVD identities.
- Stationary rank-one removal by centering, without conflating it with SVD rank rejection.
- Streamed complex covariance versus direct batch covariance, including global mean.
- Temporal/spatial AER equivalence with matched normalization and regularization.
- Exact whitening, rank-aware inversion, and the regularized eigenvalue formula.
- Local and sliding overlap/add, boundary coverage and identity controls.
- Confidence-mass normalization and complex endpoint centering of lag covariances.
- Complex AJD of known jointly diagonalizable matrices.
- Recovery of identifiable circular non-Gaussian sources by complex ICA (>0.99 matched source correlation).
- Bulk phase demodulation and exact inverse phase restoration.
- Integer rigid alignment and explicit lost-boundary validity.
- Doppler sign/scale for both IQ conventions; known multi-angle vector geometry.
- Fluid-side Poiseuille wall gradients with free and no-slip polynomial fits.
- STAP null/projector identities, distortionless constraints and rejection of overlapping steering constraints.
- Complex PCP reconstruction convergence on its low-rank-plus-sparse test.
- Downstream CCA removal of a known common nuisance with original signal amplitude preserved.
- Frozen SVD/AER/local operators applied to blood/tissue/noise reproduce the mixture output.
- Oracle truth requirements, exhaustive lower-cutoff expansion and candidate records.
- Shared velocity-readout windows with interval-matched dynamic truth.
- Strict JSON, output overwrite protection, safe NPZ round trip, local HDF5 metadata inspection.
- Reference RF DAS Doppler phase on a moving point target.

## Executed software experiments

`runs/verified_seed0`, `verified_seed1`, `verified_seed2` contain the full default
33-entry smoke matrix on distinct deterministic seeds: **99 completed method
runs**. They include centering controls, threshold controls, long/sliding SVD,
full/partial temporal/spatial AER, heuristic data-only AER regularization,
Axis 1 phase/rigid variants, Axis 2 with/without memory, ICA, PCP and three
same-case oracle searches. Each full matrix took approximately 5 seconds in
this environment. These runtimes describe a small CPU software check.

`runs/phase_diagram_validated` contains **16 controlled statistical-IQ cases**,
crossing two flow speeds, two tissue accelerations, two clutter levels, and
two seeds. Five methods per case completed: **80 additional method runs**.
The exact configuration was saved before execution. All raw metric records
are retained; no claim of superiority is inferred from this surrogate bank.

`runs/oracle_objectives_check` additionally verifies three objective-specific
searches: velocity-optimal SVD (64 candidates), WSS-optimal SVD (64 candidates),
and WSS-optimal AER (six candidates). All three searches completed. In total,
the included execution records contain 182 completed method/settings runs.

Important solver diagnostic: on the Gaussian-speckle Poiseuille surrogate,
the circular complex ICA arm reached its iteration limit in 57–59 of 60
patch-packet fits per seed. The convergence flags are retained in results.
This arm must not be interpreted as a converged ICA baseline on those cases.
Its separate identifiable, non-Gaussian source-recovery unit test converges.
Completion of a run means the software returned an output, not that every
underlying scientific assumption or iterative solver succeeded.

## Packaging

The editable install and wheel build succeed. The wheel was installed in a
separate environment and imported outside the source tree; a synthetic-IQ
SVD call succeeded. Source, tests, configurations, examples and selected
verification outputs are included in the delivery archive.

## Not yet validated

No Open H-RF RF acquisition was reconstructed or used in a filter benchmark.
Direct HTTP HDF5 metadata inspection timed out. Candidate screening uses the
consortium's public cards and metadata inventory, with sources recorded in
DATASETS.md. Native-blood/no-injection status still needs explicit confirmation.

The surrogate does not model physical RBC scattering/aggregation, a cell-free
layer, acoustic attenuation/aberration or sub-resolution wall motion. There is
no independent experimental WSS ground truth in this delivery. Registration
and component labeling for Axes 1/2, local ultrasound STAP calibration, the
original-author SSM software parity, and full 3D+t vector WSS remain research
work. All numerical and statistical checks above are implementation evidence.
