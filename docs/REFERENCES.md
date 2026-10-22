# References and source provenance

## Primary methods

- Lee SA, Jelenkovic I, Kokoeva MV, Provost J. *Spatiotemporal Adaptive
  Eigenspectrum Remapping in Transcranial Ultrasound Localization Microscopy.*
  arXiv:2610.23099v1 (2026). https://arxiv.org/html/2610.23099v1
  Implements the regularized spectral transform; absolute parameter sweeps
  from the paper are not automatically portable across covariance normalization.
- Lee SA et al. *Assessing cerebral capillary function and stalling using
  single capillary reporters in ultrasound localization microscopy.* PNAS
  123(2), e2509564123 (2026). Long-ensemble motivation, as cited in the research note.
- Demene C et al. *Spatiotemporal clutter filtering of ultrafast ultrasound data
  highly increases Doppler and fUltrasound sensitivity.* IEEE TMI 34(11),
  2271–2285 (2015). DOI:10.1109/TMI.2015.2428634.
- Baranger J et al. *Adaptive spatiotemporal SVD clutter filtering for ultrafast
  Doppler imaging using similarity of spatial singular vectors.* IEEE TMI
  37(7), 1574–1586 (2018). DOI:10.1109/TMI.2018.2789499.
  Current Baranger arm uses the full SSM and published square-fitting objective;
  the historical two-block surrogate is preserved under its explicit name.
- Song P et al. *Ultrasound small vessel imaging with block-wise adaptive local
  clutter filtering.* IEEE TMI 36(1), 251–262 (2017).
  DOI:10.1109/TMI.2016.2605819.
  The local adaptive threshold reconstruction is provisional; see the source audit.
- Maresca D, Correia M, Tanter M, Ghaleh B, Pernot M. *Adaptive Spatiotemporal
  Filtering for Coronary Ultrafast Doppler Angiography.* IEEE TUFFC 65(11),
  2201-2204 (2018). DOI:10.1109/TUFFC.2018.2870083. IEEE document 8466597.
  https://pubmed.ncbi.nlm.nih.gov/30418874/
  https://hal.science/hal-04089080/file/Maresca2018.pdf
  Identity/abstract verified; full-text algorithm verification remains pending.
- Baranger J, Aguet J, Villemain O. *Fast Thresholding of SVD Clutter Filter
  Using the Spatial Similarity Matrix and a Sum-Table Algorithm.* IEEE TUFFC
  70(8), 821-830 (2023). DOI:10.1109/TUFFC.2023.3289235.
  https://www.mathworks.com/matlabcentral/fileexchange/130139-svd_ssm_autoth
- Song P, Trzasko JD, Manduca A, Qiang B, Kadirvel R, Kallmes DF, Chen S.
  *Accelerated Singular Value-Based Ultrasound Blood Flow Clutter Filtering
  With Randomized Singular Value Decomposition and Randomized Spatial
  Downsampling.* IEEE TUFFC 64(4), 706-716 (2017).
  DOI:10.1109/TUFFC.2017.2665342. https://pubmed.ncbi.nlm.nih.gov/28186887/
- Huang C, Lok U-W, Zhang J, Liu H, Chen S. *Real-time Adaptive and Localized
  Spatiotemporal Clutter Filtering for Ultrasound Small Vessel Imaging.*
  arXiv:2405.11105 (2024). https://arxiv.org/abs/2405.11105
  EVD rules and Gaussian-weighted statistics; paper PD and IQ extension distinguished.
- Ozgun KA, Byram B. *Multidimensional Clutter Filtering of Aperture Domain
  Data for Improved Blood Flow Sensitivity.* IEEE TUFFC 68(8), 2645-2656 (2021).
  DOI:10.1109/TUFFC.2021.3073292. https://par.nsf.gov/servlets/purl/10312674
  Section III-A documents the Song noise-cutoff procedure used for reconstruction.
- Bjaerum S, Torp H, Kristoffersen K. *Clutter filter design for ultrasound
  color flow imaging.* IEEE TUFFC 49(2), 204-216 (2002).
  DOI:10.1109/58.985705. Temporal wall-filter background.
- Belouchrani A, Abed-Meraim K, Cardoso JF, Moulines E. *A blind source separation
  technique using second-order statistics.* IEEE TSP 45(2), 434–444 (1997).
  DOI:10.1109/78.554307. Complex lag-based joint diagonalization principle.
- Candes EJ et al. *Robust principal component analysis?* JACM 58(3), 11 (2011).
  https://arxiv.org/abs/0912.3599. PCP objective; complex extension in this code.
- Faure F et al. *Functional ultrasound (fUS) detects mild cerebral alterations
  using canonical correlation analysis denoising and dynamic functional
  connectivity analysis.* Imaging Neuroscience 3, IMAG.a.128 (2025).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC12406050/
- Melvin WL. *A STAP overview.* IEEE AES Systems Magazine 19(1), 19–35 (2004).
  Physical angle-time clutter manifold motivation; ultrasound implementation
  here remains a prototype with acquisition-specific steering assumptions.
- Riemer K et al. *On the use of singular value decomposition as a clutter
  filter for ultrasound flow imaging.* https://arxiv.org/abs/2304.12783

## Scattering, WSS and datasets

- Riemer K et al. *Contrast Agent-Free Assessment of Blood Flow and Wall Shear
  Stress in the Rabbit Aorta using Ultrasound Image Velocimetry.* Ultrasound
  Med Biol 48(3), 437–449 (2022). DOI:10.1016/j.ultrasmedbio.2021.10.010.
- Pourjafar-Chelikdani M et al. *Small-amplitude forced oscillation of encapsulated
  gas bubbles in elastic microchannels: effect of blood shear-thinning.*
  J Non-Newtonian Fluid Mech (2025). DOI:10.1016/j.jnnfm.2025.105514.
- Stevens TSW et al. *OpenH-RF: A Large-Scale Ultrasound Raw Channel Dataset
  Initiative.* https://arxiv.org/abs/2610.13616
- Open H-RF consortium data cards and file inventories, linked in DATASETS.md.
  Dataset license CC BY 4.0. Cite the parent initiative and the selected
  acquisition's contributor references when using experimental data.

The Axis 1/2/3 equations and proposed controls follow Henri Leroy's collaboration
note *Research Directions Beyond SVD for Blood and Tissue Separation*, recovered
from its 3 October 2026 source package. They are implemented as testable
prototypes; the note did not supply a validated classifier, motion-registration
implementation or ultrasound STAP calibration.
