# Blood/tissue separation benchmark — CPU reference pipeline

Development benchmark v0.5.0 · 23 October 2026 · prepared for Henri Leroy

The [motion-repair held-out study](docs/RESULTS_MOTION_REPAIR_V06.md) fixes
registration and gates each displacement on independent spatial pixels. Its
predeclared software gate passes on five reserved seeds, permitting a phantom
experiment while leaving significant reconstruction and blood-response losses.
The [preregistered phantom protocol](docs/PHANTOM_PROTOCOL_V07.md) fixes
physical controls, independent motion references, data eligibility and the
paired flow-response/tissue-leakage gate. It has **no experimental results yet**.
An [eSDIva RF benchmark proposal and Axis 3 literature update](docs/ESDIVA_RF_BENCHMARK_PROPOSAL.md)
lays out an acoustic-simulation tier and its separate source-truth controls;
the [Example 22 reproduction and three-view pilot](docs/RESULTS_ESDIVA_PILOT_V08.md)
record upstream Doppler checks and a short, reduced-density RF-to-IQ bridge.
The [angular calibration and full-density short control](docs/RESULTS_ESDIVA_CALIBRATION_V09.md)
admit the central `0°, -5°, +5°` reference geometry, expose an unsupported
loaded STAP fit, and compare matched paired blood-addition/background outputs
for SVD, SSM, AER and the three axes. They do not establish filter accuracy
under tissue motion, independent scatterer seeds, or in a physical phantom.
The [longer-ensemble eSDIva development controls](docs/RESULTS_ESDIVA_PHASE2_DEVELOPMENT.md)
add rigid/smooth tissue motion and fast/slow/zero flow at 10% particle density.
Guarded loaded STAP now passes nominal support and gain bounds without a
separation benefit over orthogonal STAP; the conservative motion gate falls
back on moving tissue. The [frozen protocol](docs/ESDIVA_PHASE2_PROTOCOL.md)
was applied unchanged to two fresh scatterer seeds. The
[held-out input-gate result](docs/RESULTS_ESDIVA_PHASE2_HELDOUT.md) finds
**0/2 eligible seeds**: cross-angle tissue coherence fails at +5° on one
seed and −5° on the other. No held-out filter ranking is claimed.
The [phase-3 development study](docs/RESULTS_ESDIVA_PHASE3_DEVELOPMENT.md)
tests new scatterer seeds, three-view angular calibration and a data-only
tissue-ROI selector with an RMS-bounded gain. A fixed left ROI fails fresh
development checks. A minimum-support symmetric selector passes new sparse
input checks and two full-density controls; paired source-resolved scenes still
show unrepaired RF motion tracking and no separation benefit for loaded STAP.
The new reserve remains unopened.
The [ULMShare follow-on](docs/ULMSHARE_FOLLOWON.md) is an in vivo microbubble
benchmark design for after the acoustic-simulation study.

The [shared three-axis benchmark](docs/RESULTS_SHARED_V05.md) uses the same
three-view input, disjoint-time calibration, spatial tile and score interval for
all families, including Stephen Lee AER. Its [frozen protocol](docs/SHARED_BENCHMARK_PROTOCOL_V05.md)
separates training-only operating points, matched weighted-history controls,
estimated/known motion, and evaluation-only source-identification diagnostics.

Earlier: [full-PDF method audit](docs/PDF_METHOD_AUDIT.md) and
[frozen follow-up protocol](docs/FOLLOWUP_PROTOCOL_V04.md), including Lee AER,
matched operating points, Axis 2 conditioning and a longer two-subject RF pilot.
The [previous separation-fidelity study](docs/SEPARATION_FIDELITY.md),
[first experimental RF pilot](docs/EXPERIMENTAL_PILOT.md), and
[executed verification records](docs/VALIDATION.md).
The [SVD/EVD source audit](docs/SVD_LITERATURE_AUDIT.md) records paper-level limits.

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

The demo runs 55 method/settings entries, including three truth-assisted
parameter searches. Output: `results.json`, `SUMMARY.md`, and `synthetic_iq.npz`.
The demo uses dense Gaussian-speckle **IQ surrogates**. The dedicated separation
study adds a 2-D Gaussian-PSF model with advected scatterers and moving surrounding
tissue. Neither model is a validated RF/channel or RBC simulation. The common
objective is tissue rejection with blood-information preservation. Doppler, fUSi,
ULM and WSS are downstream endpoints.

The new primary comparison uses actual output-to-blood IQ error, separately
refitted background-only output power, and refitted blood injection responses.
Frozen component diagnostics depend on the fitted operator representation:
equivalent spatial/temporal SVD or AER implementations can produce the SAME
mixture output and DIFFERENT frozen blood scores. They are secondary conditional
diagnostics, not an intrinsic ranking of adaptive filters. No gain or phase is fitted to blood truth.
Preservation-constrained oracles and per-case Pareto sets accompany IQ-MSE oracles;
an empty output is not accepted by a blood-preservation constraint.

## Implemented scope

The three research axes retain the names from the [research note](docs/notes/Research_Directions_Beyond_SVD.pdf).
Short labels such as "Axis 1" in figures and result records refer to these names.
STAP means space-time adaptive processing.

| Family | Implementation | Status |
|---|---|---|
| Global SVD | Fixed cutoff, cumulative-energy cutoff, full-rank Baranger SSM two-boundary square fit; explicit upper bound | SVD algebra and published SSM objective checked independently |
| Local SVD | Generic Gaussian patches; full-PDF Song block thresholds and Eq.1 normalization | Unspecified numerical tolerances logged; older sigmoid reconstruction retained as historical |
| Cardiac sliding SVD | Maresca normalized-spectrum curvature plus fixed operator offset | Original manual calibration is explicit; offset is not universally transferable |
| Randomized SVD | Seeded leading-clutter range finding and optional spatial subsampling | Known-source rejection/preservation checked; sampling variant documented |
| Adaptive EVD | Median-eigenvalue, image-NCC and sigmoid rules; Gaussian covariance localization | Huang equations; paper PD combination is separate from the IQ extension |
| Oracle SVD | Declared bands; IQ, preservation-constrained tissue rejection, velocity or WSS objectives | Evaluation-only; every candidate and infeasible search retained |
| Long ensemble | Same SVD on a longer contiguous window | Ensemble length explicit; no joining across gaps |
| ZCA-AER | Lee centered temporal covariance /P and streamed spatial covariance /(T-1) | Published algorithms checked; absolute and relative regularizers distinguished |
| AER ablations | Centering on/off; full/partial whitening; relative, absolute, or noise-floor regularization; oracle grid | Parameter rules and gain changes recorded |
| [Axis 1: Motion-compensated virtual long ensemble](docs/METHODS.md#axis-1-motion-compensated-virtual-long-ensemble) | Original quality-gated patch memory; shared runner adds finite-history expiration, subpixel smooth registration and exact weighted-SVD controls | Research prototype; matched memory is algebraically equivalent to weighted SVD |
| [Axis 2: Motion-demodulated multi-lag separation](docs/METHODS.md#axis-2-motion-demodulated-multi-lag-separation) | Complex lag AJD, rank/whitening controls, optional spectral-norm cap, lag-signature diagnostics | Conditioning controls do not guarantee source identification |
| Other baselines | Centered complex kurtosis FastICA with explicit convergence fallback, complex PCP/RPCA; Butterworth controls | ICA stopping/conditioning diagnostics; nonlinear RPCA injection responses |
| [Axis 3: Local ultrasound STAP](docs/STAP_BENCHMARK.md) | Physical angle-time steering, tissue projector, loaded covariance adaptation, LCMV, guarded training | 42-case statistical-IQ study; uncompounded views and independent tissue calibration |
| Downstream fUS | Faure-style whitening/common-mode CCA regression | Separate power-Doppler denoising API; outside IQ benchmark |

Every SVD/AER centering flag means temporal per-pixel mean subtraction **without
mean add-back**. The default is uncentered. Centered and uncentered results are
kept as different experimental conditions.

## Run on your own data

Create an `Acquisition` using `btbench.core` and save with
`btbench.io.save_acquisition`; see [the data contract](docs/DATA_CONTRACT.md).

```bash
btbench write-config configs/my_benchmark.json --prf-hz 3000
btbench run data/unfiltered_iq.npz --config configs/my_benchmark.json --output runs/real_pilot
```

Settings are plain JSON. `configs/benchmark.json` contains the default matrix.
`configs/oracle_examples.json` shows independent objective-specific searches.
Always choose a new output directory; existing benchmark results are protected.

An executable 16-case Cartesian sweep is also provided:

```bash
python examples/run_phase_diagram.py configs/phase_diagram_small.json runs/my_grid
```

The new follow-ups have frozen configurations and preserve per-condition results:

```bash
btbench fidelity-followup --config configs/fidelity_followup_v04.json --output runs/my_followup
python examples/prepare_openhrf_followup.py configs/openhrf_followup_v04.json configs/openhrf_waterloo_acq0_geometry.json data/my_followup
btbench openhrf-followup data/my_followup/waterloo_acq0.npz data/my_followup/waterloo_acq12.npz --config configs/openhrf_followup_v04.json --output runs/my_real_followup
```

Run the shared benchmark (40 synthetic cases × 73 settings, with three
counterfactual fits per setting; four CPU workers, one BLAS thread each).
Install the plot extra to regenerate figures:

```bash
python -m pip install -e '.[plot]'
btbench shared-study --config configs/shared_axes_v05.json --output runs/my_shared_study
python scripts/plot_shared_benchmark.py runs/my_shared_study
python scripts/verify_results.py --manifest runs/shared_v05_manifest.json
python scripts/verify_shared_results.py
```

The motion-repair study fixes a new split of three development and five
reserved held-out seeds across the same eight conditions. Its 36 matched
settings cover no compensation, gated rigid, gated smooth and the
evaluation-only prescribed warp. See
[the protocol](docs/MOTION_REPAIR_PROTOCOL_V06.md) for the data-only acceptance
gate and exact release criteria. After the archived run is installed, verify it:

```bash
btbench motion-study --config configs/motion_repair_v06.json --output runs/my_motion_study
python scripts/verify_results.py --manifest runs/motion_v06_manifest.json --require-current-code
python scripts/verify_motion_results.py
```

The score uses all three complex views on frames [64,128), with 16-frame packet
latency and declared 16/64-frame history. A single setting per family, motion arm
and support is selected across all training conditions; held-out seeds never
change that choice. The scene specifies angular IQ and an image warp, with a
shared calibration tissue texture. It does not provide acoustic/experimental
validation. Portable results, selection candidates, CSV tables and four PNG/SVG
figures are included; source and artifact hashes are checked by CI.

Earlier studies, with their historical source fingerprints and interpretation limits:

```bash
btbench separation-study --config configs/separation_fidelity.json --output runs/my_separation
btbench literature-study --config configs/svd_literature.json --output runs/my_literature
btbench stap-study --config configs/stap_first_study.json --output runs/my_stap
```

The separation study has 30 cases and 39 settings, including three evaluation-only
oracles. It crosses statistical IQ and PSF-transport conditions, includes motion
and memory ablations, and records ICA convergence and RPCA counterfactual responses.
The literature study has 24 cases and 20 arms. The STAP study has 42 cases,
12 comparison outputs including the truth-assisted SVD band oracle, and a
separate LCMV velocity-bank detector. All STAP comparators receive the same
angular data. Calibration-demodulated SVD receives the same independent
calibration packet as STAP. These are controlled algorithmic experiments.

The older 16-case phase diagram crosses speed, tissue acceleration, clutter level and seed with five methods,
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
A bounded 64-frame Waterloo CCA RF window has now been reconstructed and passed
through a 58-setting engineering pilot: 32 calibration frames and 32 evaluation
frames. Both views and their offsets are preserved. The short pilot supplies
no independent blood/tissue truth and does not rank experimental separation accuracy.
The [v0.3.0 follow-up](docs/RESULTS_FOLLOWUP_V04.md) now reconstructs 1,024 frames
from each of two subjects, with 256 calibration and 768 evaluation frames per
subject. It retains the same limits on experimental source/phase truth.

For large raw files, start with metadata and a bounded window:

```bash
btbench inspect-openhrf /path/Acq0.hdf5 --output header.json
```

The inspector supports HTTP range-readable URLs. Direct range reads and
metadata inspection now work on the pinned Waterloo file. A bounded
plane-wave RF-to-IQ adapter is provided in `examples/beamform_openhrf_window.py`.
It requires documented transmit time origins and receive-aperture mapping before
execution; the Waterloo contract checks recorded delays and event timing. Production reconstruction should retain complex IQ from a verified
provider/zea pipeline, before envelope detection and clutter filtering.

## What remains experimental

The original Axis 1 packet API has global phase and integer rigid registration.
The shared runner adds subpixel translation and an eight-coefficient smooth
deformation model, with known-warp and round-trip controls. General nonrigid
motion, acoustic realism and transfer to independent tissue/calibration remain
unvalidated. Axis 2 and ICA have an explicit,
untrained component score; its thresholds must be frozen or calibrated on
training cases. Supplied full texts now support the Song/Maresca/Baranger audit.
Song's unspecified detector tolerances and author-software execution parity
remain qualifications; see the [PDF audit](docs/PDF_METHOD_AUDIT.md).
STAP is a far-field angle-time prototype, with full receive-aperture
and near-field manifold validation still needed. WSS is a known-direction,
fluid-side 2-D comparator with supplied wall geometry; the package does not
claim a validated 3D+t vector-WSS reconstruction.

See [methods](docs/METHODS.md), [validation](docs/VALIDATION.md), and
[references](docs/REFERENCES.md).

CI checks numerical tests on Python 3.10/3.12, source/wheel contents and an installed
wheel outside the source tree. Tested direct numerical dependencies are pinned.

The repository includes compact verification records in `runs/`. Generated IQ
arrays, acquisition files, and package build products are excluded from Git.
Run the demo commands above to regenerate the synthetic arrays.


The [v0.3.0 results and figures](docs/RESULTS_FOLLOWUP_V04.md) include 3,850
synthetic and 58 experimental outputs. Lossless XZ archives preserve the exact
frozen JSON bytes in small transfer parts. Verify and restore them with:

```bash
python scripts/verify_results.py --manifest runs/followup_v04_manifest.json --require-current-code --restore
```

Read restored `results.json.xz` with `lzma.open(path, 'rt')`. Runners write
ordinary JSON when rerun. Historical v0.2.0 GZIP records remain available:
`python scripts/verify_results.py --restore` verifies those earlier source hashes
without requiring them to match the current implementation.
