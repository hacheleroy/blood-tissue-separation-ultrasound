# References and source provenance

## SVD and related low-rank methods


- Demene C et al. [*Spatiotemporal clutter filtering of ultrafast ultrasound data highly increases Doppler and fUltrasound sensitivity.*](https://doi.org/10.1109/TMI.2015.2428634) IEEE TMI 34(11),
  2271–2285 (2015).

- Song P et al. [*Ultrasound small vessel imaging with block-wise adaptive local clutter filtering.*](https://doi.org/10.1109/TMI.2016.2605819) IEEE TMI 36(1), 251–262 (2017).
  The local adaptive threshold reconstruction is provisional; see the source audit.

- Song P, Trzasko JD, Manduca A, Qiang B, Kadirvel R, Kallmes DF, Chen S.
  [*Accelerated Singular Value-Based Ultrasound Blood Flow Clutter Filtering With Randomized Singular Value Decomposition and Randomized Spatial Downsampling.*](https://doi.org/10.1109/TUFFC.2017.2665342) IEEE TUFFC 64(4), 706-716 (2017).

- Baranger J et al. [*Adaptive spatiotemporal SVD clutter filtering for ultrafast Doppler imaging using similarity of spatial singular vectors.*](https://doi.org/10.1109/TMI.2018.2789499) IEEE TMI
  37(7), 1574–1586 (2018).
  Current Baranger arm uses the full SSM and published square-fitting objective;
  the historical two-block surrogate is preserved under its explicit name.

- Maresca D, Correia M, Tanter M, Ghaleh B, Pernot M. [*Adaptive Spatiotemporal Filtering for Coronary Ultrafast Doppler Angiography.*](https://doi.org/10.1109/TUFFC.2018.2870083) IEEE TUFFC 65(11),
  2201-2204 (2018). IEEE document 8466597.
  Identity/abstract verified; full-text algorithm verification remains pending.

- Ozgun KA, Byram B. [*Multidimensional Clutter Filtering of Aperture Domain Data for Improved Blood Flow Sensitivity.*](https://doi.org/10.1109/TUFFC.2021.3073292) IEEE TUFFC 68(8), 2645-2656 (2021).
  Section III-A documents the Song noise-cutoff procedure used for reconstruction.

- Baranger J, Aguet J, Villemain O. [*Fast Thresholding of SVD Clutter Filter Using the Spatial Similarity Matrix and a Sum-Table Algorithm.*](https://doi.org/10.1109/TUFFC.2023.3289235) IEEE TUFFC
  70(8), 821-830 (2023).

- Riemer K et al. [*On the use of singular value decomposition as a clutter filter for ultrasound flow imaging.*](https://arxiv.org/abs/2304.12783)

- Zhang N, Nguyen MB, Mertens L, Barron DJ, Villemain O, Baranger J.
  [*Improving coronary ultrafast Doppler angiography using fractional moving blood volume and motion-adaptive ensemble length.*](https://doi.org/10.1088/1361-6560/ac7430)
  Phys Med Biol 67(12), 125021 (2022). Motion-dependent sliding SVD
  windows; a manually selected myocardium region and per-window cutoff adjustment
  distinguish this published comparator from an automatic registration-first arm.

- Lee SA et al. [*Assessing cerebral capillary function and stalling using single capillary reporters in ultrasound localization microscopy.*](https://doi.org/10.1073/pnas.2509564123) PNAS
  123(2), e2509564123 (2026). Long-ensemble motivation, as cited in the research note.

- Candes EJ et al. [*Robust principal component analysis?*](https://arxiv.org/abs/0912.3599) JACM 58(3), 11 (2011).
  PCP objective; complex extension in this code.

- Bayat M, Fatemi M. [*Concurrent Clutter and Noise Suppression via Low Rank Plus Sparse Optimization for Non-Contrast Ultrasound Flow Doppler Processing in Microvasculature.*](https://doi.org/10.1109/ICASSP.2018.8461638)
  ICASSP, 1080–1084 (2018). Native-blood low-rank/sparse comparator;
  the current generic complex PCP is not a verified reproduction.

- Ashikuzzaman M, Belasso C, Kibria MG, Bergdahl A, Gauthier CJ, Rivaz H.
  [*Low Rank and Sparse Decomposition of Ultrasound Color Flow Images for Suppressing Clutter in Real-Time.*](https://doi.org/10.1109/TMI.2019.2941865)
  IEEE TMI 39(4), 1073–1084 (2020; online 2019). RAPID uses fast robust
  matrix completion/in-face Frank–Wolfe; not implemented by the PCP arm.

- Insana MF, Dai B, Babaei S, Abbey CK.
  [*Combining Spatial Registration With Clutter Filtering for Power-Doppler Imaging in Peripheral Perfusion Applications.*](https://doi.org/10.1109/TUFFC.2022.3211469)
  IEEE TUFFC 69(12), 3243–3254 (2022). Direct precedent for registering
  frames before PCA filtering in contrast-free perfusion imaging.

## CCA, ICA and eigen-based separation

Temporal Gram EVD and SVD describe equivalent temporal subspaces. Papers that
present eigenvalue processing are grouped here by their implementation.


- Lee SA, Jelenkovic I, Kokoeva MV, Provost J. [*Spatiotemporal Adaptive Eigenspectrum Remapping in Transcranial Ultrasound Localization Microscopy.*](https://arxiv.org/html/2610.23099v1) arXiv:2610.23099v1 (2026).
  Implements the regularized spectral transform; absolute parameter sweeps
  from the paper are not automatically portable across covariance normalization.

- Huang C, Lok U-W, Zhang J, Liu H, Chen S. [*Real-time Adaptive and Localized Spatiotemporal Clutter Filtering for Ultrasound Small Vessel Imaging.*](https://arxiv.org/abs/2405.11105) arXiv:2405.11105 (2024).
  EVD rules and Gaussian-weighted statistics; paper PD and IQ extension distinguished.
  See [the additional comparator audit](ADDITIONAL_COMPARATORS.md) for its
  implemented subset, tuning requirements and the three other supplied papers.

- Belouchrani A, Abed-Meraim K, Cardoso JF, Moulines E. [*A blind source separation technique using second-order statistics.*](https://doi.org/10.1109/78.554307) IEEE TSP 45(2), 434–444 (1997).
  Complex lag-based joint diagonalization principle.

- Faure F et al. [*Functional ultrasound (fUS) detects mild cerebral alterations using canonical correlation analysis denoising and dynamic functional connectivity analysis.*](https://pmc.ncbi.nlm.nih.gov/articles/PMC12406050/) Imaging Neuroscience 3, IMAG.a.128 (2025).

## STAP and angle-time processing


- Melvin WL. [*A STAP overview.*](https://doi.org/10.1109/MAES.2004.1263229) IEEE AES Systems Magazine 19(1), 19–35 (2004).
  Physical angle-time clutter manifold motivation; ultrasound implementation
  here remains a prototype with acquisition-specific steering assumptions.

## Temporal wall filtering


- Bjaerum S, Torp H, Kristoffersen K. [*Clutter filter design for ultrasound color flow imaging.*](https://doi.org/10.1109/58.985705) IEEE TUFFC 49(2), 204-216 (2002).
  Temporal wall-filter background.

## WSS, blood scattering and contrast-agent effects


- Riemer K et al. [*Contrast Agent-Free Assessment of Blood Flow and Wall Shear Stress in the Rabbit Aorta using Ultrasound Image Velocimetry.*](https://doi.org/10.1016/j.ultrasmedbio.2021.10.010) Ultrasound
  Med Biol 48(3), 437–449 (2022).

- Pourjafar-Chelikdani M et al. [*Small-amplitude forced oscillation of encapsulated gas bubbles in elastic microchannels: effect of blood shear-thinning.*](https://doi.org/10.1016/j.jnnfm.2025.105514) J Non-Newtonian Fluid Mech (2025).

## Datasets and acquisition resources


- Stevens TSW et al. [*OpenH-RF: A Large-Scale Ultrasound Raw Channel Dataset Initiative.*](https://arxiv.org/abs/2610.13616)

- Open H-RF consortium. [*Data cards and file inventories.*](DATASETS.md)
  Dataset license CC BY 4.0. Cite the parent initiative and the selected
  acquisition's contributor references when using experimental data.

## Collaboration note and proposed research axes


- Leroy H. [*Research Directions Beyond SVD for Blood and Tissue Separation.*](notes/Research_Directions_Beyond_SVD.pdf)
  Collaboration note, 3 October 2026. Original supplied PDF.

The three research axes retain the note's names:

- [Axis 1: Motion-compensated virtual long ensemble](METHODS.md#axis-1-motion-compensated-virtual-long-ensemble)
- [Axis 2: Motion-demodulated multi-lag separation](METHODS.md#axis-2-motion-demodulated-multi-lag-separation)
- [Axis 3: Local ultrasound STAP](STAP_BENCHMARK.md)

Their equations and proposed controls follow this collaboration note.
They are implemented as testable prototypes; the note did not supply a validated
classifier, motion-registration implementation or ultrasound STAP calibration.
