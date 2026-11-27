# Blood/tissue separation for ultrasound

Python benchmarks for suppressing tissue clutter while preserving blood information in unfiltered complex ultrasound IQ.

> **Experimental repository — under active development.** Algorithms, APIs and results may change. The methods are research prototypes and are not validated for clinical use.

## Main findings

- **Motion correction helps in selected simulations**, but slow-flow recovery remains difficult and physical validation is still pending. [RF evaluation](docs/ESDIVA_PHASE11_PROSPECTIVE_RESULT.md)
- **Near-wall residual background is mainly tissue.** Preserving more blood can also preserve clutter; low background power alone does not demonstrate successful separation. [Source analysis](docs/RECEIVE_CORE_SOURCE_PROFILE_V1.md)
- **The tested RPCA variants do not reliably recover blood waveforms.** The Cauchy-core implementation also fails its convergence criterion at the fixed iteration budget. [RPCA comparison](docs/RPCA_COMPARISON_V1.md)

## Methods

Baselines include global/local SVD, adaptive EVD, ZCA-AER, RPCA, ICA and Butterworth filtering. Three experimental research directions are also included:

1. **Motion-compensated virtual long ensemble**
2. **Motion-demodulated multi-lag separation**
3. **Local ultrasound STAP** (space-time adaptive processing)

The benchmarks measure blood-signal fidelity, tissue leakage and downstream Doppler/near-wall velocity errors using statistical IQ models, acoustic RF simulations and limited experimental pilots. They do not establish a validated 3D+t wall-shear-stress pipeline.

## Quick start

Requires Python 3.10 or newer. From the repository directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test,io]'
btbench demo --seed 0 --output runs/my_demo
```

On Windows, use WSL or the corresponding virtual-environment activation command. The demo runs on CPU and writes results, a summary and synthetic IQ. Use a new output directory for each run.

For your own data, start with the [complex-IQ data contract](docs/DATA_CONTRACT.md).

## Documentation

- [Methods](docs/METHODS.md) and [references](docs/REFERENCES.md)
- [Detailed benchmark guide and study history](RESEARCH_ARCHIVE.md)
- [Validation records](docs/VALIDATION.md) and [dataset notes](docs/DATASETS.md)
- [Planned physical flow-on/off experiment](docs/TRANSPORT_LAG_PHYSICAL_PROTOCOL_V1.md)

The research archive preserves the previous README exactly for historical verification. Individual study reports contain their assumptions, failures and reproduction commands.

## License

[MIT](LICENSE) — Copyright (c) 2026 Henri Leroy. Referenced datasets and third-party software retain their own licenses.
