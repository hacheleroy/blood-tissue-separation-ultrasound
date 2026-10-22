# SVD/EVD literature audit

Audit date: 22 October 2026. Matrix convention: **X = U diag(sigma) V^H**,
with complex IQ, spatial rows and slow-time columns. `H` is conjugate transpose.
Every retained band uses zero-based half-open indices `[k,h)`. The default is
uncentered. Centering subtracts each pixel's temporal mean without adding it back.

This is a source and implementation audit of the benchmark's principal method
families. The table records the degree of reproduction; it is not an exhaustive
systematic review of every SVD-related publication since 2015.

| Source / family | Implementation | Evidence and status |
|---|---|---|
| Demene et al., TMI 2015, DOI 10.1109/TMI.2015.2428634 | `svd_filter`, fixed lower/upper bands, global short/long windows, explicit centering | Reference SVD projection. Full retention, known bands, complex reshape and frozen operators tested. |
| Song et al., TMI 2017, DOI 10.1109/TMI.2016.2605819 | `song_local_filter`: overlapping rectangular spatial blocks, independently selected lower and upper cutoffs, separate incoherent power aggregation | **Provisional threshold reconstruction.** Local organization is documented in the original abstract. A slope cutoff and sigmoid-frequency / linear singular-tail fit are implemented. Numerical fit choices are explicit, but original-paper/software parity is not established. |
| Maresca et al., TUFFC 2018, DOI 10.1109/TUFFC.2018.2870083; IEEE 8466597 | Paper-specific implementation **pending full-text verification** | Identity, authors, journal and abstract verified. HAL lists an author-deposited PDF, but HAL downloads returned HTTP 403 and IEEE returned HTTP 418. Existing generic energy/sliding SVD controls are not labeled as reproductions of this paper. |
| Baranger et al., TMI 2018, DOI 10.1109/TMI.2018.2789499 | `spatial_similarity`, `svd_filter(selector="baranger")` | Spatial-envelope Pearson SSM uses the full available rank and retains its diagonal. The original 2018 paper was read from the University of Geneva repository. |
| Baranger, Aguet, Villemain, TUFFC 2023, DOI 10.1109/TUFFC.2023.3289235 | `square_fit` / `baranger_band` | Implements the authors' documented two-square `corr2` objective with an O(r^2) summed-area-table search. Independent literal-template tests verify the objective and both boundaries. This is an independent Python implementation, without a MATLAB author-code parity claim. |
| Temporal singular-vector spectral thresholds, surveyed by Baranger 2018 | `selector="temporal_frequency"` | Threshold on the PSD-weighted **mean absolute frequency**. This is an explicit folded-spectrum convention; positive and negative Doppler do not cancel. The full set of 14 estimators in that paper is not claimed here. |
| Song et al., TUFFC 2017, DOI 10.1109/TUFFC.2017.2665342 | `randomized_svd_filter` | Leading-clutter randomized range finding, oversampling and QR-stabilized power iterations. The uncomputed residual is preserved. Algebraic full-range parity and known-source rejection/preservation are tested. The optional spatial-sampling variant transfers an estimated basis to all pixels; it does not reproduce the authors' parallel partition/reassembly scheme. |
| Huang et al., arXiv:2405.11105 (2024), Eqs. 7-9 | `huang_evd_filter` | Complex temporal-Gram EVD; median-eigenvalue, adjacent-spatial-image NCC and sigmoid weights. Uses **lambda = sigma^2**, not sigma. The median rule is a flow-filter heuristic in this source; it is not the Gavish-Donoho denoising guarantee. |
| Huang et al. 2024, Eqs. 11-13 | `gaussian_local_evd` | Full-FOV Gaussian windows enter **before** covariance estimation: X^H W^2 X. The published normalized power combination is returned separately. The complex-IQ overlap-add output is explicitly named an extension and is evaluated separately from the paper's power output. Gaussian widths/centers are frozen configuration choices. |
| Bjærum et al., TUFFC 2002, DOI 10.1109/58.985705; standard Butterworth wall filter | `temporal_highpass` | Complex, per-pixel slow-time high-pass. Three fixed cutoff arms and one causal arm are in the default matrix. Analytic cutoff gain, DC rejection, Doppler-sign symmetry and frozen linearity are checked. |
| Long ensemble, Lee et al., PNAS 2026, DOI 10.1073/pnas.2509564123 | `sliding_filter(method="svd")` | Contiguous long/sliding ensembles, explicit observation length and hop. The long-ensemble operator is implemented; this small study does not reproduce the paper's 6000-frame contrast-enhanced ULM experiments. |

## Corrections to the initial package

The old `ssm_surrogate` compared two blocks at a single lower cutoff, discarded
the diagonal and searched at most 64 modes. It remains callable for historical
reproduction, but is removed from the default matrix. The new Baranger arm searches
both tissue/blood and blood/noise boundaries over all available modes. The old
Gaussian patch overlap/add arm remains `svd_local`; it is a generic local-SVD
control. It must not be cited as a reproduction of Song's complete adaptive method.

The previous implementation had a Butterworth function but did not include it in
`standard_methods()` or the default configuration. The new matrix includes 10,
30 and 80 Hz cutoffs at a 1000 Hz revisit rate, plus a causal 30 Hz arm. These
values are frozen examples, not parameters optimized on the evaluation bank.

## Baranger objective and indexing

Let C = corr(|U|). For tissue rank k and upper blood boundary h, define

\[
M_{kh}(i,j) = \mathbf{1}_{i,j<k} + \mathbf{1}_{k\le i,j<h}.
\]

The implementation maximizes Pearson correlation of the vectorized matrices:

\[
J(k,h)=\frac{A_{kh}-S_C m_{kh}/r^2}
 {\sqrt{V_C\left(m_{kh}-m_{kh}^2/r^2\right)}},\qquad
m_{kh}=k^2+(h-k)^2,
\]

where A is the sum of C inside both squares, S_C is the total sum of C, and
V_C = sum(C-mean(C))^2. Integral-image rectangle sums compute A in constant
time per candidate. Both high-correlation squares and the remaining noise
region contribute to the objective through its normalization. No singular modes
are truncated to make thresholding cheaper. Spatially constant envelopes and
an uninformative SSM trigger documented, finite fallbacks.

The author's description sometimes labels square size and rank boundary
differently. This package always returns `[k,h)`, with no shared blood/tissue or
blood/noise endpoint. Literal template tests cover the indexing convention.

## Song verification limit

The original full article and original author software could not be retrieved.
The accessible original abstract establishes local blocks and adaptive cutoffs.
Ozgun and Byram (2021), Section III-A, describes Song's noise cutoff procedure:
a sigmoid fit to singular-vector frequencies, then a line fit to the noise
singular-value tail. That sequence is implemented here. The choices of log-slope
normalization, sigmoid bounds, a relative tail-deviation tolerance, fit failure
fallback and coherent IQ reconstruction are fully specified in the code and
results. Those choices are **not verified as the original Song implementation**.
The arm name `song_local_reconstruction_*` and diagnostics preserve that limit.
Do not use its observed performance as a performance claim about Song's paper.

## Maresca 2018: identified; algorithm verification pending

The user supplied the DOI that resolves [IEEE document 8466597](https://ieeexplore.ieee.org/abstract/document/8466597):
David Maresca, Mafalda Correia, Mickael Tanter, Bijan Ghaleh and Mathieu Pernot,
*Adaptive Spatiotemporal Filtering for Coronary Ultrafast Doppler Angiography*,
IEEE TUFFC **65**(11), 2201-2204 (2018),
[DOI 10.1109/TUFFC.2018.2870083](https://doi.org/10.1109/TUFFC.2018.2870083).
The original abstract describes detecting slow flow in rapidly moving myocardium.
PubMed and the authors' laboratory publication list verify the bibliographic record;
the indexed publisher PDF URL contains document number 08466597.

HAL's metadata API identifies deposit `hal-04089080v1` and the public
[Maresca2018.pdf](https://hal.science/hal-04089080/file/Maresca2018.pdf).
Full-text retrieval remained blocked: HAL returned HTTP 403 and IEEE HTTP 418.
The abstract does not specify the numerical threshold rule, normalization or
windowing needed for a faithful implementation. No paper-specific arm is claimed.
A readable copy of the original PDF is needed to finish this algorithm check.

## Further families identified, outside this implementation

Multi-rate HOSVD with optimized tensor-core classification (Kim et al. 2017/2018),
aperture-domain HOSVD (Ozgun and Byram 2021), subspace tracking, and 3-D feature
clustering are recorded as candidates. They require their own classifiers, data
dimensions or validation. Aperture-domain methods cannot be reconstructed from
compounded IQ. Noise equalization and incoherent power processing can improve PD
appearance without defining a phase-preserving blood-IQ operator; WSS cannot be
assigned to such outputs without a separate velocity estimator.

## Sources actually consulted

- [Baranger 2018: author-deposited full article](https://access.archive-ouverte.unige.ch/access/metadata/88011ea7-bad6-4e52-bee4-2fdc8137efc3/download).
- [Baranger's SVD_SSM_AutoTh description and objective](https://www.mathworks.com/matlabcentral/fileexchange/130139-svd_ssm_autoth?s_tid=srchtitle).
- [Song 2017: PubMed / original abstract](https://pubmed.ncbi.nlm.nih.gov/27608455/).
- [Song et al. 2017 randomized SVD and spatial downsampling](https://pubmed.ncbi.nlm.nih.gov/28186887/).
- [Maresca 2018: original abstract and metadata](https://pubmed.ncbi.nlm.nih.gov/30418874/);
  [authors' publication list](https://www.physicsformedicine.espci.fr/all-publications/).
  The HAL PDF was located through repository metadata but could not be downloaded.
- [Huang et al. 2024: full preprint](https://arxiv.org/pdf/2405.11105), especially Eqs. 7-13.
- [Ozgun and Byram 2021: author-deposited full article](https://par.nsf.gov/servlets/purl/10312674), especially Section III-A and Eqs. 2, 12-13.
- Henri Leroy, *Research Directions Beyond SVD for Blood and Tissue Separation*,
  3 October 2026, supplied six-page note. The original attachment is unchanged.

No third-party source code or full article PDF is redistributed in this repository.
