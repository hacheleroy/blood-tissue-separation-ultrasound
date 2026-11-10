# Stable RPCA → Sui Cauchy core → static batch IRLS: IQ benchmark v1

## Scope and fidelity

This is a **source-separated IQ surrogate benchmark**, not an acoustic RF/channel,
RBC, fluid–structure, or in-vivo experiment. The earlier project's RF
checkpoints were not present in the recoverable workspace. We used the existing
`poiseuille` statistical-IQ and `advected_scatterers` 2-D PSF transport
generators and their ten conditions from `configs/separation_fidelity.json`.
Each scene has 96 complex frames of 8×13 pixels at 1000 Hz PRF and 5 MHz
center frequency. The final 16 frames are scored. Three untouched evaluation
seeds (61, 73, 89) yield 30 scenes; the four-condition development bank at
seed 17 selects one fixed parameter set per method. Masks and source truth are
used for tuning and evaluation, never passed into a fit.

The methods, in requested order:

1. **Plain stable RPCA:** convex quadratic-residual nuclear plus complex-ℓ1
   model, with alternating proximal block minimization. Blood candidate is
   its sparse component; residual models noise.
2. **Sui Cauchy core:** reweighted singular-value and entry soft thresholds
   with ADMM on the entire Casorati matrix. This follows the published core
   updates, but omits randomized spatial downsampling and parallel reassembly.
   The Cauchy scale is 0.1 after observed-IQ RMS normalization; it is not the
   paper's unit-dependent 0.01. The original method has no explicit quadratic
   noise residual.
3. **Static batch IRLS reconstruction:** complex Hermitian extension of the
   factorized batch IRLS RPCA updates of Rontogiannis, Giampouras, and
   Koutroumbas (2020). Chen et al. (2025) identify a static IRLS-rPCA arm, but
   their full update rules and initialization were unavailable for verification.
   **These numbers are not a reproduction of Chen's exact static filter**, and
   no online update or spatial-temporal constraint is claimed.

All three use source-blind full-batch inputs. The fit for a counterfactual
`tissue + noise` is rerun independently; the paired blood response is
`F(blood + tissue + noise) − F(tissue + noise)`. The existing full-window
SVD (discard first mode), PCP, and unfiltered input are fixed references.

## Frozen selection

The four candidates per method and objective are in
`configs/rpca_comparison_v1.json`; the complete development observations are
in `results/rpca_comparison_v1/development.json`. The development objective
minimizes the mean of
`log1p(observed near-wall IQ NRMSE) + 0.5 log1p(paired near-wall IQ NRMSE) +
0.2 log1p(100 × tissue-only background power ratio)`.
`selected.json` stores the selected parameters and SHA-256 digests of both
configs and solver source; the evaluation command rejects a mismatch.

| Arm | Frozen parameters |
| --- | --- |
| Stable RPCA | `tau=0.1, sparse_lambda=0.05` |
| Sui Cauchy core | `sparse_lambda=0.03, gamma=0.1` |
| Batch IRLS reconstruction | `low_rank_lambda=0.1, sparse_mu=0.01, rank=12, epsilon=0.02` |

## Independent evaluation

The table reports medians over 30 cases. Lower IQ NRMSE and tissue background
power are better; higher velocity coverage is better. Observed IQ NRMSE uses
the actual output against blood truth. Paired IQ NRMSE measures the
counterfactual blood response and does **not** measure remaining tissue in the
actual output. Tissue background power ratio uses a separate tissue-only
refit in the tissue ROI, normalized by the input tissue power. Velocity
coverage is the fraction of near-wall pixels with output power at least 1% of
the blood-only power and a finite lag-one velocity. All ratios are linear.

| Method | Near-wall observed IQ NRMSE ↓ | Paired near-wall IQ NRMSE ↓ | Tissue-only background power ratio ↓ | Near-wall velocity coverage ↑ | Converged fits / 30 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Stable RPCA | 0.93 | 0.83 | 1.4×10⁻⁶ | 0.91 | 20 |
| Sui Cauchy core | 0.86 | 0.77 | 1.1×10⁻⁵ | 1.00 | **0** |
| Batch IRLS reconstruction | 0.87 | 0.85 | 5.0×10⁻⁶ | 1.00 | 24 |
| SVD, discard first mode | 1.04 | 0.29 | 5.1×10⁻⁵ | 1.00 | closed-form |
| Existing PCP | 0.97 | 0.97 | 3.0×10⁻⁶ | 0.56 | 30 |
| Unfiltered | 8.51 | 0.00 | 1.00 | 1.00 | identity |

Median observed **flow-wide** IQ NRMSE was 0.91 (stable), 0.79 (Cauchy),
and 0.85 (batch IRLS); these are large errors, so none demonstrates reliable
blood waveform recovery. The zero paired error of the unfiltered input is
mathematically expected and underscores why observed error is also reported.

| Condition, median over three evaluation seeds | Stable | Cauchy core | Batch IRLS | SVD |
| --- | ---: | ---: | ---: | ---: |
| Slow blood overlapping tissue frequency | 0.97 | 0.73 | 0.70 | 0.72 |
| PSF scatterers, stationary tissue | 0.37 | 0.23 | 0.40 | 0.17 |
| PSF scatterers, bulk tissue motion | 1.00 | 1.00 | 1.00 | 3.59 |
| PSF scatterers, deforming tissue | 1.00 | 1.01 | 1.00 | 3.14 |
| PSF scatterers, equal static components | 1.00 | 1.00 | 1.00 | 1.00 |

The table's condition values are **observed near-wall IQ NRMSE**; 1.00 is
essentially no recovered blood. Under PSF bulk motion the stable method
covered only 0.19 of near-wall pixels (median); batch IRLS covered 0.50 and
Cauchy 0.69. In equal-static scenes none identifies the blood source.
Median near-wall velocity errors are available per case in `evaluation.json`,
but must be read with coverage and the blood-only PSF Doppler bias: a velocity
estimate from a weak residual can look accurate while the blood IQ is lost.

**Cauchy convergence caveat.** At the frozen 120-iteration budget the published
core updates did not satisfy the 10⁻⁵ residual-and-change rule in any case.
Median relative primal residual is 2.68×10⁻⁴ and median successive sparse
change is 0.142. On the slow-overlap evaluation case with seed 61, increasing
the cap from 120 to 1200 did not resolve the approximately 0.140 change: the
updates settle into an oscillation. Its IQ scores are therefore
iteration-limited diagnostics, not certified converged decompositions.

## Run and inspect

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/benchmark_rpca.py development
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/benchmark_rpca.py evaluation
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p test_rpca_comparators.py -v
```

`results/rpca_comparison_v1/evaluation.json` records all 180 scene-method
rows, method diagnostics, per-case IQ/background/power/velocity measures and
summary medians. Whole-fit timing is recorded as a local implementation
diagnostic, not a paper runtime comparison: the Cauchy arm does not use Sui's
acceleration and the IRLS arm is not Chen's online implementation.

## Literature

- Sui et al., *IEEE TUFFC* 69(8), 2425–2436 (2022), DOI
  [10.1109/TUFFC.2022.3180053](https://doi.org/10.1109/TUFFC.2022.3180053).
- Chen et al., *IEEE TMI* 44(6), 2477–2491 (2025), DOI
  [10.1109/TMI.2025.3535550](https://doi.org/10.1109/TMI.2025.3535550).
- Rontogiannis, Giampouras, and Koutroumbas, *IEEE Signal Processing Letters*
  27, 1340–1344 (2020), DOI
  [10.1109/LSP.2020.3011896](https://doi.org/10.1109/LSP.2020.3011896).
