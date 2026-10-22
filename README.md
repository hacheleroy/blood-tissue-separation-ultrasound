# Blood/tissue separation benchmark — CPU reference pipeline

Version 0.1.0 · 22 October 2026 · prepared for Henri Leroy

Executable research code implementing the collaboration note and the subsequent
AER/whitening discussion. The default entry point takes **unfiltered complex
beamformed IQ**. It never treats log-compressed B-mode, power Doppler, or a
vendor velocity map as IQ or as independent ground truth.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test,io]'
python -m pytest -q
btbench demo --seed 0 --output runs/my_demo
```

On Windows use WSL, or activate the environment with the corresponding Windows
command. CPU NumPy/SciPy are sufficient. MATLAB and a GPU are not required for
the numerical tests and statistical-IQ demo.

The demo runs 33 method/settings entries, including three truth-assisted
parameter searches. Output: `results.json`, `SUMMARY.md`, and `synthetic_iq.npz`.
The included synthetic data are dense Gaussian-speckle **IQ surrogates**. They
are useful for software verification; they are not a validated acoustic/RBC
simulation or evidence of WSS performance in vivo.

## Implemented scope

| Family | Implementation | Status |
|---|---|---|
| Global SVD | Fixed cutoff, cumulative-energy cutoff, envelope-SSM two-block cutoff; explicit upper bound | Reference algebra; SSM rule is a documented surrogate |
| Local SVD | Overlapping Gaussian spatial patches | Reference implementation of blockwise filtering; no claim of exact Song code reproduction |
| Oracle SVD | Declared lower/upper-cutoff grid; separate IQ, velocity, or WSS objectives | Evaluation-only; same-case truth-assisted bound |
| Long ensemble | Same SVD on a longer contiguous window | Ensemble length explicit; no joining across gaps |
| ZCA-AER | Temporal and spatial Gram/covariance formulations; exact streamed spatial moments | Equation-level implementation of regularized AER |
| AER ablations | Centering on/off; full/partial whitening; relative, absolute, or noise-floor regularization; oracle grid | Parameter rules and gain changes recorded |
| Axis 1 | Bulk phase correction, optional rigid integer registration, Gaussian patches, q-normalized covariance memory, event reset/fallback | Research prototype |
| Axis 2 | Motion demodulation, complex lag covariances, rank-reduced whitening, complex Hermitian-pair AJD, soft reconstruction | Research prototype; component score heuristic |
| Other baselines | Complex kurtosis FastICA, complex PCP/RPCA, temporal Butterworth high-pass | Reference algorithms; assumptions documented |
| Axis 3 | Physical angle-time steering, tissue projector, loaded covariance adaptation, LCMV, guarded training | Separate prototype requiring uncompounded views and verified geometry |
| Downstream fUS | Faure-style whitening/common-mode CCA regression | Separate power-Doppler denoising API; outside IQ benchmark |

Every SVD/AER centering flag means temporal per-pixel mean subtraction **without
mean add-back**. The default is uncentered. Centered and uncentered results are
kept as different experimental conditions.

## Run on your own data

Create an `Acquisition` using `btbench.core` and save with
`btbench.io.save_acquisition`; see [the data contract](docs/DATA_CONTRACT.md).

```bash
btbench write-config configs/my_benchmark.json
btbench run data/unfiltered_iq.npz --config configs/my_benchmark.json --output runs/real_pilot
```

Settings are plain JSON. `configs/benchmark.json` contains the default matrix.
`configs/oracle_examples.json` shows independent objective-specific searches.
Always choose a new output directory; existing benchmark results are protected.

An executable 16-case Cartesian sweep is also provided:

```bash
python examples/run_phase_diagram.py configs/phase_diagram_small.json runs/my_grid
```

It crosses speed, tissue acceleration, clutter level and seed with five methods,
saving the configuration before execution. Extend the generator grid for SNR,
radial amplitude variation and motion-event controls. This remains a statistical
IQ software experiment; realistic acoustic identifiability maps are a later step.

No reference field reaches the data-only filtering API. Oracle flags, selected
parameters, candidate losses, centering, runtime, input hash, and numerical
environment are recorded. An in-vivo reference map remains a comparator.

## Dataset findings

[Open H-RF screening and recommended acquisitions](docs/DATASETS.md).
The most promising reviewed native-blood candidates are Waterloo carotid and
femoral-vein RF sequences. Absence of contrast injection is inferred from their
descriptions and still needs explicit confirmation. Their stored velocity maps
are processed estimates, with no independent near-wall/WSS truth identified.
No Open H-RF acquisition was beamformed or benchmarked in this delivery.

For large raw files, start with metadata and a bounded window:

```bash
btbench inspect-openhrf /path/Acq0.hdf5 --output header.json
```

The inspector also supports HTTP range-readable URLs. Remote inspection was
attempted here but timed out; dataset findings are based on the public data
cards and the consortium's per-file metadata inventory. A small-window
plane-wave RF-to-IQ example is provided in `examples/beamform_openhrf_window.py`.
It requires manually verified time origins and receive-aperture mapping before
execution. Production reconstruction should retain complex IQ from a verified
provider/zea pipeline, before envelope detection and clutter filtering.

## What remains experimental

Axis 1 currently has global phase and integer rigid 2-D registration; affine
and nonrigid tissue warps are not implemented. Axis 2 and ICA have an explicit,
untrained component score; its thresholds must be frozen or calibrated on
training cases. The SSM threshold is not the published three-region sum-table
algorithm. STAP is a far-field angle-time prototype, with full receive-aperture
and near-field manifold validation still needed. WSS is a known-direction,
fluid-side 2-D comparator with supplied wall geometry; the package does not
claim a validated 3D+t vector-WSS reconstruction.

See [methods](docs/METHODS.md), [validation](docs/VALIDATION.md), and
[references](docs/REFERENCES.md).

The repository includes compact verification records in `runs/`. Generated IQ
arrays, acquisition files, and package build products are excluded from Git.
Run the demo commands above to regenerate the synthetic arrays.
