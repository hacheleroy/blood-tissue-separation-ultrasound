# Additional comparators: slow flow, automation and tissue motion

26 October 2026. The four user-supplied references are assessed against the
current benchmark. The first two are **different low-rank/sparse algorithms**;
the third is partly implemented; the fourth is an unimplemented motion-adaptive
SVD comparator. Insana et al. (2022) is included as direct prior art for the
registration-before-filtering experiment. Adding a citation or a comparator
plan does not retroactively add a method to a frozen numerical result.

| Paper | Actual idea | Current coverage | Fair comparator and limitation |
|---|---|---|---|
| [Bayat & Fatemi 2018](https://doi.org/10.1109/ICASSP.2018.8461638) | Separate strong tissue as low rank and weak native blood as sparse while suppressing noise. The accessible original abstract reports simulations against singular-value thresholding. | Generic complex PCP/RPCA is implemented, but its exact `X=L+S` constraint and parameterization are not verified as the paper's algorithm. | Keep PCP explicitly generic until the full article supplies an auditable objective, noise handling and solver settings. Test dense speckle and slow, tissue-like blood, where sparsity or source identifiability may fail. |
| [Ashikuzzaman et al. 2020 / RAPID](https://doi.org/10.1109/TMI.2019.2941865) | Adapt fast robust matrix completion with an in-face extended Frank–Wolfe solver to extract moving blood without a manually selected clutter-rank cutoff. Original abstract reports simulation, phantom and in-vivo validation; the [author supplement](https://users.encs.concordia.ca/~hrivaz/Ashik_RAPID_TMI_Supp.pdf) compares a 15-frame setting. | **No RAPID implementation.** Generic complex PCP uses another optimization problem and solver. | Verify the full algorithm or author code and input/normalization convention before adding a RAPID arm; freeze any optimization tolerances and resource limits. An automatic rank rule does not imply parameter-free operation or guaranteed preservation of phase-sensitive IQ. |
| [Huang et al. 2024](https://arxiv.org/pdf/2405.11105) | Temporal-Gram EVD with adaptive hard/soft thresholds and full-FOV Gaussian-weighted local covariance. The spatial tissue-versus-blood image NCC minimum is the paper's threshold rule without a predefined *scaling factor*. Other rules and Gaussian geometry involve choices. Its local outputs are combined in power Doppler (PD). | `huang_evd_filter` covers median, image NCC and sigmoid rules; `gaussian_local_evd` computes weighted covariance, returns the published-style PD combination separately and labels its coherent IQ overlap/add as an **extension**. It does not cover every proposed threshold or claim GPU/runtime parity. | Compare NCC with matched full-FOV support and the same fixed Gaussian geometry. Evaluate both paper PD for detection and the IQ extension for blood-phase preservation. Report the user choices (packet, number/width of windows, mode-neighborhood size, and any sigmoid/median scale) and runtime. The paper selects several best CR/CNR values by sweeping a scale; those are not untuned operating points. |
| [Zhang et al. 2022, with Baranger](https://doi.org/10.1088/1361-6560/ac7430) | Estimate tissue motion in a manually selected myocardial ROI and adapt overlapping SVD windows until maximum cumulative displacement is about one wavelength. The cutoff uses Maresca's singular-value acceleration with an operator-adjusted offset **per ensemble**, not Baranger's spatial-similarity threshold. Attenuation compensation for fractional moving blood volume (FMBV) is a separate calibrated measurement step. | Maresca cutoff and fixed/sliding ensembles exist separately. There is **no faithful motion-adaptive-window arm**. Axis 1 registers complex data and accumulates weighted history; it is not the same method. | Reproduce window definition and cutoff with the published manual controls clearly labeled, compare with a fixed window on identical support, and separately test a fully data-derived/no-manual-ROI version as an extension. A CNR improvement does not show that all slow blood is preserved. |
| [Insana et al. 2022](https://doi.org/10.1109/TUFFC.2022.3211469) | Register echo frames **before** PCA clutter filtering for native-blood muscle perfusion, showing narrower clutter spectra and limitations from out-of-plane motion and blood/noise overlap. | Axis 1 and shared SVD experiments use IQ motion correction, but have not established a calibration-free B-mode-estimated registration front end on physical data. | Treat pre-filter registration as established prior art. Compare registered SVD/Baranger with unregistered SVD/Baranger and the Zhang adaptive-window strategy, using an identical post-registration support mask. |

## Proposed benchmark additions (not yet run)

1. **Low-rank/sparse arm:** retain generic complex PCP under its own name. Add
   RAPID and a Bayat-specific arm only after original methods/code can be
   verified. Use the same short/long packet support and include a dense-flow
   counterexample, zero-flow tissue, and motion where tissue and blood overlap.
2. **Motion arms:** estimate smooth frame displacement from B-mode/envelope of
   the *same* acquisition, apply the resulting warp to complex IQ, check
   residual tissue phase, and restore lab-frame coordinates/phase after
   filtering. Compare fixed SVD, Baranger SSM, and (separately) Zhang's
   published motion-adaptive windows; keep manual-ROI and automatic variants
   distinct. Include an out-of-plane or decorrelated-speckle failure control.
3. **Equal evaluation:** report incremental blood-IQ preservation and residual
   tissue power on independently source-labeled *simulations*; noise
   amplification, flow visibility, known flow speed/tube geometry and static
   no-flow controls on a physical phantom; and CNR/PD only as secondary image
   metrics. Known physical pump speed and tube geometry do **not** provide a
   pixelwise blood/tissue IQ decomposition. Record acquisition length, valid
   spatial support, user-set parameters, prior calibration and processing time.
4. **Protocol:** choose methods, operator controls, scales and motion quality
   rules on development data, then freeze them before opening new evaluation
   data. Existing phase-11 results remain unchanged; none of these references
   supplies a new result for that bank.

The Bayat conference paper and RAPID main article were identified from their
original abstracts; the RAPID author supplement and the complete Huang
preprint were readable. Zhang's abstract and author-provided full-text preview
support its window/cutoff summary; Insana's full text is available via its
open author manuscript. Paper-specific numerical fidelity for Bayat and RAPID
still requires their complete algorithm/code. Copyrighted article files are
not redistributed.
