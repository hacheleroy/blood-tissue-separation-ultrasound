# Mathematical conventions and benchmark protocol

## Matrix and centering

The public IQ shape is `(time,z,x)`. Internally X has shape `(pixels,time)`.
H is conjugate transpose. No input samples are conjugated by reshaping.

Uncentered spatial second moment: R0=XX^H/T. Temporal Gram: C0=X^H X/P.
Optional temporal centering uses Xc=X-(X 1/T)1^T. It is a separate operation
that can remove low-frequency blood as well as stationary tissue. Filtering
acts on Xc and the removed mean is not returned. These controls do not equate
centering with rejection of precisely the first SVD mode.

## SVD, locality, and threshold rules

X=U diag(sigma) V^H. Zero-based `discard=k`, `upper=h` retains [k,h):
Y=U[:,k:h] diag(sigma[k:h]) V[:,k:h]^H. An absent upper bound means all
remaining modes, including noise, are retained. Short/global means full spatial
support over a short window; long/global has the same support over a longer
contiguous interval. Local SVD uses positive Gaussian overlap/add weights.
The `sliding` wrapper adds temporal Gaussian overlap/add with an explicit hop;
it accepts SVD or AER so matched sliding-window comparisons can be configured.

Energy threshold: reject the smallest leading block carrying the configured
fraction (default 95%) of total squared singular values; retain at least one
mode. This is data-only and can remove genuine blood.

Historical SSM surrogate (`ssm_surrogate`, removed from the default matrix):
form Pearson correlations between absolute spatial singular
vectors. For each split k, maximize the average within-block off-diagonal
correlation minus the between-block correlation. Integral-image sums accelerate
the search, capped at 64 modes. This implements a transparent two-block
criterion inspired by spatial-similarity filtering. It does not replicate the
published Baranger square-fitting cutoff/sum-table method. Publication comparisons
using this historical arm must retain its surrogate label.

The current `baranger` arm implements the published two-square correlation
objective with a summed-area-table search over both cutoffs. Its SSM includes
the diagonal and all available modes. See the equations and source-fidelity
table in [the literature audit](SVD_LITERATURE_AUDIT.md). New randomized SVD,
complex temporal-Gram EVD, Huang eigenvalue/NCC/sigmoid rules, weighted
Gaussian covariance localization and the provisional Song-style local arm
are also specified there. None accesses the evaluation truth.

## AER/ZCA

For C=V diag(lambda) V^H, W=V diag((lambda+epsilon)^(-alpha)) V^H.
Temporal application: Y=X W. Spatial application: Y=W X, with the spatial
eigenbasis. Full whitening has alpha=1/2; partial whitening has 0<alpha<1/2.
Exact zero regularization uses a pseudoinverse on nonzero eigenmodes, so the
output moment is identity only when full rank; otherwise it is a projector.

With alpha=1/2, output eigenvalues are lambda/(lambda+epsilon). Increasing
epsilon eventually approaches a scalar multiple of the original data. Zero
regularization is the strongest equalization, subject to rank and precision.
Whitening alone makes no tissue/blood label assignment inside a mixed mode.

Regularization choices:

- `epsilon_abs`: numerical units of the specified covariance, for equation-level replication.
- `epsilon_rel`: fraction of the leading eigenvalue, default 1e-3; scale-invariant benchmark convention.
- `epsilon_rule=noise_floor`: median of the lower half of numerical nonzero covariance eigenvalues. This is a heuristic; its assumption can fail for flow-dominated spectra.

The temporal covariance denominator is P, spatial denominator T. Their raw
epsilon values and output gains therefore differ. `normalization=common` uses
T on both sides: for identical support, centering, and regularization, the two
filtered matrices then coincide. This is tested explicitly. Spatial accumulation
uses Chan/Welford merged complex moments and global centering, not separately
centered packet covariances. It needs two passes to estimate then apply the
fixed transform. It does not create additional information by changing sides.

The shipped demo uses 64-frame SVD packets and 256-frame long-SVD/AER windows.
Use additional matched-window arms to isolate filtering effects from observation
time. It is an offline benchmark. Report latency and ensemble support separately
from kernel wall time. A dense P-by-P covariance is capped at 512 MiB before
allocation; full 3-D global spatial whitening is generally impractical.
RSS is recorded as the process's cumulative high-water mark, including the
input/truth arrays and previous methods. It is a resource diagnostic, not an
isolated per-method peak-memory comparison; those require separate processes.

## Motion compensation and covariance memory (Axis 1)

Estimate increments from a separately specified tissue ROI:
delta_phi[t]=arg sum_roi X[t+1] conj(X[t]). Integrate and demodulate each
frame. If no ROI is supplied, the brightest quartile of time-mean power is
used; this assumption is logged through the method description.

Optional integer rigid registration aligns envelope texture to the first frame
using FFT cross-correlation, with a maximum shift. Complex interpolation has
zero boundaries; there is no periodic wrapping. Deformation and subpixel rigid
refinement are future extensions. The reference coordinate system must remain
consistent across all packets; covariance states cannot be pooled across
independently drifting coordinate systems.

For each Gaussian patch, accumulate raw pair moments with
S_tau <- rho S_tau+(1-rho) q Rhat_tau and m <- rho m+(1-rho)q.
Use S_tau/m; for centered operation subtract the matching endpoint means.
The confidence mass prevents artificial power drift when q varies. Low
confidence clears memory and invokes a short local-SVD fallback for that packet.
The fixed leading tissue subspace of the accumulated R0 is projected out.
The prototype uses the minimum within-packet phase coherence as q. Its choice
and reset threshold need validation, particularly for phase wraps and deformation.

After filtering, restore the estimated tissue phase and inverse rigid warp to
lab coordinates. Without restoration, the Doppler estimate would be relative
to the demodulated tissue frame. Boundary validity fraction is reported;
boundary losses remain in evaluation and are not silently masked using truth.

## Multi-lag complex SOBI (Axis 2)

R(tau)=mean_t x(t+tau) x(t)^H for lags 1,2,4,8,16 by default. Whiten a
rank-reduced R0 and jointly diagonalize both Hermitian parts of every lag:
(R+R^H)/2 and (R-R^H)/(2i). The imaginary part carries Doppler direction.
A complex unitary Jacobi algorithm minimizes off-diagonal structure. Records
include diagonalization residual, retained rank, lag coherence, phase, spatial
extent, and soft component weights. Packet memory can be enabled independently
of motion demodulation; the default Axis 2 arm uses current-packet statistics.

The initial blood weight combines lag-one coherence (noise rejection), near-DC
phase/coherence (tissue score), and spatial participation. It is an explicit
heuristic, not a proven classifier. Higher-lag signatures determine the joint
basis; a validated multi-feature component-identification model remains needed.
Blood with the same spatial mixture and lag signature as tissue stays
unidentifiable. Complex kurtosis FastICA uses the same soft classifier as a
separate baseline, allowing the representation and identification steps to be
distinguished. The ICA baseline assumes circular complex, non-Gaussian sources.

## Other baselines and Axis 3

PCP/RPCA solves min ||L||_*+lambda||S||_1, X=L+S with complex soft-thresholding
and inexact augmented-Lagrangian updates. S is candidate blood. Its sparse
pixel-time assumption can fail for dense blood speckle. Convergence is logged.
Butterworth high-pass filtering is in the default matrix, with three fixed
cutoff frequencies and an initialized causal control. The offline
forward/backward arm squares the one-pass amplitude response (one-half
amplitude at its design cutoff); the causal arm has 1/sqrt(2) amplitude there.
Short ensembles are rejected when they cannot support the stated padding.
Frequency-based filters validate their PRF against the acquisition revisit rate.

STAP requires `(time,view,z,x)` data, actual Tx/Rx directions and view offsets.
Its steering phase is 2*pi*f0/c*(n_tx+n_rx).v*t, with explicit IQ phase sign.
Implemented: orthogonal tissue null, loaded oblique covariance-adaptive null,
guard-excluded local training, and LCMV response constraints. The adaptive arm
requires an independently defined tissue training ROI. The same coherent view
mean must form the input to other methods. Full channel/near-field manifolds,
aberrated transmit fields and calibrated view gains need dedicated validation.
The new [STAP benchmark](STAP_BENCHMARK.md) preserves each view, uses absolute
snapshot times and independent tissue calibration, and includes calibration,
training-support and identifiability controls. LCMV detection is separate
from phase-preserving IQ reconstruction and WSS estimation.

Optional downstream CCA works on power-Doppler movies, using separate functional
and nuisance ROIs. It whitens each dataset, SVDs their concatenation, selects
common temporal modes by an explicit correlation threshold, and projects them
out of original amplitudes. It is separate from IQ-level blood/tissue filtering.

## Measurements and oracle rules

Scalar Doppler uses arg mean(y[t+1] conj(y[t])) times c*PRF/(2*pi*f0*q), where
q=(n_tx+n_rx).d for a known flow direction. Default q=2 gives axial velocity.
Multi-view velocity solves the 2-D least-squares geometry explicitly. Revisit
PRF, pulse PRF, sample rate and output-map rate are different quantities.
Aliasing is not automatically corrected. Complex IQ still has a Nyquist
interval; frequencies beyond PRF/2 do not define a noise-only band.

WSS fits a fluid-side polynomial of directional velocity along supplied inward
wall normals. Free-intercept and no-slip options are available in the API.
Viscosity, wall positions, normals and wall velocity are explicit SI inputs.
The comparator is 2-D and directional. It makes no 3D tensor-WSS claim.

Oracle search uses true blood IQ, true velocity, or true WSS on the same case.
It is an optimistic bound over the declared grid, not a deployable method or a
universal bound on every possible SVD estimator. `discard="all"` in an oracle
grid expands all lower cutoffs; upper bounds are optional. The demo searches
only cutoffs 0..8 for speed. Distinct velocity, IQ and WSS oracles can choose
different settings. Real in-vivo processed flow maps must never enter truth slots.

Metrics: raw IQ error, velocity NRMSE and bias, near-wall NRMSE, directional WSS
error and invalid fraction, tissue-ROI residual-power proxy, and timing. For
SVD/AER/local linear maps, frozen-operator component power is also measured:
apply the SAME transform fitted to the mixture to known blood and tissue.
Independently refitting on isolated sources would be a misleading leakage test.
Other methods explicitly report these component diagnostics as unavailable.
AER changes amplitude calibration; no truth-based gain correction is hidden.

The output also includes a shared 64-frame measurement-window series for every
method. Filtering ensemble length and velocity-readout support are independent
controls. Dynamic true velocity fields can be supplied as (T,Z,X); each readout
is compared to the prescribed mean over that same interval. Full-record metrics
are descriptive averages and cannot resolve pulsatility. Temporal variability
is labeled unattributed because actual pulsatility also contributes. Shorter
readout windows need their own noise/gradient observability analysis.

Freeze preprocessing, wall/ROI geometry, velocity estimator, fit support and
search grids. Split by subject/acquisition; pixels from one acquisition cannot
form independent training/test subjects. Retain failed and invalid estimates.
Sweeps should cross relative velocity, tissue acceleration/deformation, clutter
level, SNR, radial amplitude variation, packet size, ensemble length and phase
convention. A clinically useful conclusion requires independent reference
measurements and realistic acoustic scattering beyond the provided smoke tests.
