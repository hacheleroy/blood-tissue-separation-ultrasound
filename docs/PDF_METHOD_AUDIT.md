# Full-PDF method audit, 23 October 2026

The four supplied PDFs were read locally, including the equations and the
published MATLAB appendices. Their byte hashes and the consulted AER source
are in [the provenance record](PDF_SOURCE_PROVENANCE_V04.json). Copyrighted
article PDFs and third-party package source are not redistributed.

| Source | Checked rule | Implementation and remaining qualification |
|---|---|---|
| Song et al., TMI 2017, pp.252-253, Figs.2-3 and Eq.1; DOI 10.1109/TMI.2016.2605819 | Lower cutoff is the larger of the singular-value flattening point and lag-one temporal-singular-vector frequency crossing. Upper cutoff uses a noise-frequency plateau to choose a straight-line fit to the LOG singular-value tail. Normalize each block by the SUM of retained singular values, coherently average overlapping blocks, then calculate PD. | New `song_pdf_local` implements this order. The PDF does not give numerical slope/plateau/deviation detectors, so their tolerances are explicit configuration choices. The older sigmoid/linear-raw-tail `song_local` remains available solely as the historical reconstruction. It is not the full-PDF arm. `normalization=none` is a labeled amplitude-preserving ablation. |
| Maresca et al., TUFFC 2018, pp.2202-2203, Eqs.2-3; DOI 10.1109/TUFFC.2018.2870083 | Normalize singular values by their SUM, convert to dB, find the maximum interior second derivative, add a fixed operator-selected offset k, and reject that many leading modes. Their k=26 was manually adjusted to an initial cutoff of 30 on a 60-frame ensemble; overlap was about 25%. | New `selector=maresca` with a 60-frame/45-frame-hop comparison. The manual offset is exposed; k=26 is not assumed optimal for another acquisition. There is no invented upper noise cutoff. Finite differences, log floor, clipping and one-based-to-zero-based conversion are recorded. Weighted IQ overlap/add is an explicit reconstruction extension beyond the paper's window-wise angiograms. |
| Baranger et al., TMI 2018, Eq.5 and subspace segmentation; DOI 10.1109/TMI.2018.2789499 | Pearson correlation of the envelopes of spatial singular vectors. Fit two adjacent correlated squares for tissue/blood and blood/noise boundaries. | Existing `baranger` implementation matches this objective. It keeps the full available rank and the SSM diagonal. Constant-envelope and uninformative-matrix fallbacks are explicit. This does not reproduce every one of the other threshold estimators surveyed in the article. |
| Baranger, Aguet and Villemain, TUFFC 2023, Eqs.18-21 and Appendix B; DOI 10.1109/TUFFC.2023.3289235 | Sum-table evaluation of the SAME normalized two-square correlation, not a new tissue/blood criterion. | The existing O(r^2) search was checked against the original PDF's sum-table equation and independent literal corr2 templates. Rank boundaries use [k,h): reject the first k modes, retain k through h-1. MATLAB package execution/bitwise parity is not claimed; degenerate cases and deterministic ties are documented. |
| Stephen A. Lee et al., arXiv:2610.23099v1, Eqs.9-15 and Algorithms 1-2 | Center each pixel over time. Temporal covariance is X^H X/P. The streamed spatial algorithm accumulates a GLOBAL mean and uses (sum xx^H - T mean mean^H)/(T-1). Apply regularized ZCA with weights (lambda+epsilon)^(-1/2). | New `lee_aer` uses these two conventions, including the spatial T-1 denominator. Absolute 1/26/1000 regularizers retain the current data's units; they do not reproduce the paper's optimized animal ULM operating points. Relative regularization and unit-operator scaling are separately labeled extensions. Common-denominator controls check algebraic equivalence; no fitted truth-based gain is used. |

The MathWorks [SVD_SSM_AutoTh v1.0.2 description](https://www.mathworks.com/matlabcentral/fileexchange/130139-svd_ssm_autoth)
was inspected. The downloadable package requires sign-in and its public source
viewer returned no code. The supplied 2023 PDF itself contains the author's
optimized MATLAB algorithm in Appendix B, sufficient to audit the objective.
No authentication is needed to continue the independent implementation.

## AER interpretation and acquisition limits

The [Lee preprint](https://arxiv.org/html/2610.23099v1) studies contrast-enhanced
transcranial ULM, long ensembles and localization/tracking. Our native-blood-like
statistical scenes and bounded carotid RF pilot test the transform and transfer;
they do not reproduce those animal experiments or their trajectory counts.

Whitening changes amplitude calibration. Raw IQ errors are retained alongside
normalized spectral/phase diagnostics and incremental source responses. Plot
normalization, where used for visibility, is labeled and does not change scores.
Complex IQ retains a signed Nyquist interval: bins above PRF/2 in an unshifted
FFT correspond to negative Doppler frequencies and are not automatically noise.
No such one-sided noise label is used in this benchmark.

## Correction to the frozen-component interpretation

Equal filtered mixtures can have unequal frozen source attributions. In
particular, temporal and spatial implementations of the same SVD projection,
and common-normalization AER, give equal F(X) while their fitted left/right
operators can give different results on B alone. Axis 1 without memory is also
equivalent to local SVD at matching support. Its previously reported frozen
blood advantage is therefore not, by itself, evidence of a better filter.

Tests explicitly demonstrate these identities and check that refitted injection
responses agree. The [new protocol](FOLLOWUP_PROTOCOL_V04.md) uses actual output
error and separately refitted background response for its primary calibrated
comparison. Conditional frozen attribution remains visible as a secondary
diagnostic, without a representation-independent preservation claim.
