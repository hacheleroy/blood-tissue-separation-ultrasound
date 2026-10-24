# ULMShare as the next benchmark tier

Status: design only. Complete and freeze the eSDIva acoustic-simulation study
before selecting ULMShare acquisitions or tuning another algorithm on them.

[ULMShare](https://github.com/provostultrasoundlab/ulmshare) describes 99 in vivo
transcranial acquisitions from 61 mice. Its repository holds metadata summaries,
MATLAB reconstruction examples and selected density maps and trajectories;
the [raw RF and metadata archive](https://doi.org/10.20383/103.01550) is hosted
separately at FRDR. These are contrast-enhanced microbubble ULM acquisitions,
which are relevant to the experimental setting of Stephen Lee's AER paper.
They do **not** give isolated native red-cell blood IQ, tissue IQ or known
displacement. Consequently ULMShare can test downstream microbubble recovery,
reproducibility and processing robustness, but cannot directly furnish the
source-resolved blood preservation and tissue leakage scores from eSDIva.

## Entry conditions after eSDIva

1. Freeze the eSDIva reconstruction, matched-frame protocol, gates and selected
   SVD, SSM, AER and Axis 1/2/3 settings. Publish both passed and failed
   motion/slow-flow strata, without tuning on an evaluation seed.
2. Read ULMShare's summary and a small acquisition's metadata first. Verify
   RF dimensions, acquisition timestamps, angles, transducer, frame grouping,
   PRF, receive-channel geometry, gain and probe-motion information. Choose
   distinct mice and sessions for development and evaluation before viewing
   their algorithm outputs. Record download and derived-data checksums.
3. Establish the supplied reconstruction as a reproducibility control on a
   small subset. Retain one identical beamformed IQ stream, temporal support,
   emission budget and calibration for all methods. If only one compounded
   stream is reconstructible, compare single-stream SVD and AER there and mark
   three-view Axis 3 inapplicable; never synthesize missing angles.
4. Compare the three research axes under their names: Axis 1:
   Motion-compensated virtual long ensemble; Axis 2: Motion-demodulated
   multi-lag separation; Axis 3: Local ultrasound STAP. Add global/local SVD,
   SSM and Stephen Lee AER on matched input. Report a calibrated false-detection
   control, the retained candidate distribution, localization/track consistency,
   vessel coverage and repeatability between independent buffers and mice.
   A supplied ULM map or tracks are a reference pipeline output, not truth.

Do not treat a better-looking density image, an AER microbubble result, or
agreement with the supplied tracks as proof of native red-cell blood fidelity.
Keep the physical, independently measured phantom protocol as a separate tier.
