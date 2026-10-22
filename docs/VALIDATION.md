# Verification performed for version 0.1.0

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
reference three-region SSM algorithm, and full 3D+t vector WSS remain research
work. All numerical and statistical checks above are implementation evidence.
