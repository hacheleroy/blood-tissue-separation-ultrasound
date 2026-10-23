# Open H-RF acquisition screening

Reviewed on 22 October 2026: all 46 dataset/subdataset cards in the official
GitHub catalog, followed by the acquisition inventories of the flow candidates.
GitHub snapshot: `5eae31257f55391129c5aadb178659948448da5d`.
Direct HDF5 reads below use Hugging Face revision
`05234aa084b27c8e4d0096c3a9adc33fbbc123d8`.

The target is tissue rejection with preservation of the blood signal. In-vitro
flow with an independent reference is preferred, followed by in-vivo recordings
with tissue motion; absence of bubbles is preferred. A geometry mask, prescribed
pump voltage or an ultrasound-derived velocity map supplies a different level
of evidence from an independent flow measurement or source-separated IQ truth.

## In-vitro candidates first

| Candidate | Acquisition and contrast | Reference status | Use |
|---|---|---|---|
| [Twente vortex-flow phantom](https://github.com/open-h/Open-H-RF/blob/main/datasets/twente-vortexflow/README.md) | Six acquisitions; short-pulse and chirp tracks; 750 frames/track at about 500 Hz. Explicit microbubbles and optical glass beads. | The card describes synchronized camera recordings. The pinned folder contains six HDF5 files and illustration assets; the representative HDF5 contains RF and ultrasound images, with no raw camera sequence or velocity field found. Optical registration/PIV are not provided. | Separate contrast arm only. Raw optical-reference access is still a gap. The curved array emits a diverging wave and needs its own verified beamforming model. |
| [ResolveStroke flow phantom](https://github.com/open-h/Open-H-RF/blob/main/datasets/resolvestroke/phantom_flow/README.md) | Explicit microbubble bolus; five 1-s clips at 4 kHz, including pre-bolus baseline. Matrix-array DDC channel IQ. | Tube masks provide geometry. Device-computed velocity/MVI maps are processed references, not independent physical truth. The released clips do not reproduce the full-bolus proprietary processing. | Separate contrast/flow-on/off study. Pre-bolus frames may provide a negative control; native-blood or blood-mimicking scattering there is unconfirmed. |
| [Oslo SWE/ARFI phantoms](https://github.com/open-h/Open-H-RF/blob/main/datasets/oslo/F_motion/README.md) | Four 200-frame high-rate tissue push/tracking acquisitions. | Moving tissue with no blood-flow target; no independent blood component. | Tissue-motion negative control. Respect push/track boundaries. A future semi-synthetic blood-injection test must be labeled separately from physical flow truth. |

No documented bubble-free in-vitro Doppler acquisition with an independent flow
reference was identified. The Twente card's optical-reference description should
not be treated as an already usable released velocity ground truth. A display
GIF is insufficient for calibrated velocimetry.

## In-vivo shortlist

| Candidate | Acquisition | Contrast and reference status | Priority |
|---|---|---|---|
| [Waterloo carotid](https://github.com/open-h/Open-H-RF/blob/main/datasets/waterloo-carotid/README.md) | 93 acquisitions, eight subjects; 30,000–36,000 frames each. Longitudinal CCA and cross-sectional ICA, baseline and maneuvers. | Native-blood imaging is inferred; injection status is not explicitly stated. Supplied VFI maps use a wall-filtered pipeline. | First pilot: nonmultiplexed baseline CCA, then maneuvers. |
| [Waterloo femoral vein](https://github.com/open-h/Open-H-RF/blob/main/datasets/waterloo-femoralvein/README.md) | 15 participants/acquisitions, 12,000 frames each, two transmits and 192 elements; muscular contraction with 40-degree head-up tilt. | Injection status unspecified. Derived/dealiased VFI maps can be noise-sensitive. | Slow-flow and tissue-motion challenge after CCA. |
| [TU/e carotid](https://github.com/open-h/Open-H-RF/blob/main/datasets/tue-carotid/README.md) | 150 frames at 30 Hz; 128 focused events interleaved with 21 plane waves. | No independent source/flow reference identified. | Secondary motion/imaging source. Pulse rate must not be confused with same-view revisit rate. |

Waterloo femoral Acq0 has about 2,850-Hz revisit timing and 4.2 seconds of RF;
other files differ. Some Waterloo US4R-lite acquisitions use receive
multiplexing; avoid them until their event mapping is verified. Split by subject
across Waterloo contributions because participants can overlap.

[ULMShare](https://github.com/open-h/Open-H-RF/blob/main/datasets/ulmshare/README.md)
is an in-vivo contrast arm for later ULM transfer. Twente microbubble simulations
and passive cavitation recordings do not provide a native-blood Doppler core.
Single/few-frame B-mode phantoms, CT geometry references and static Field II
examples do not establish a Doppler flow reference.

## Directly verified Waterloo Acq0 and completed pilot

[Pinned acquisition](https://huggingface.co/datasets/nvidia/Open-H-RF/blob/05234aa084b27c8e4d0096c3a9adc33fbbc123d8/waterloo-carotid/data/Acq0.hdf5).
The actual HDF5 header confirms baseline longitudinal CCA, US4R, 30,000 frames,
two plane-wave transmits at -10/+10 degrees and 128 receive channels. The file
uses 5-MHz center frequency, 32.5-MHz sampling and 100-ns A/D start time. Event
intervals are 166.667 microseconds; same-view revisit rate is approximately
3 kHz. Stored element delays fit the supplied plane-wave time origin with a
maximum residual of 7.76e-6 samples.

Only frames 0:64 were decoded: 128 MiB of real RF, followed by complex DAS on a
96-by-128 grid. Both views, their 166.667-microsecond offset and coherent-mean
baseline are embedded in NPZ. The [experimental pilot](EXPERIMENTAL_PILOT.md)
uses 32 calibration frames and 32 disjoint evaluation frames. It establishes
executable processing, not native-blood confirmation or separation accuracy.
No supplied velocity/PD field is promoted to truth.

The measured B-mode log-amplitude spatial correlation with the provider image
is about 0.700, after interpolation to the exported grid. Receive steering and
apodization differ. Structural agreement does not establish experimental Doppler
phase accuracy. The known moving-point numerical phase test is a separate check.

## Next experimental requirements

1. Obtain a bubble-free blood/blood-mimicking flow phantom with independently
   controlled flow and repeatable tissue-motion/background conditions.
2. Verify receive event mapping, timing, beamforming phase and amplitude from
   calibration acquisitions before quantitative experimental comparisons.
3. Freeze acquisition/subject splits and operating points; include slow flow,
   motion and weak-signal regions. Keep processed flow maps as comparators.
4. Retain disjoint in-vivo acquisitions for robustness and downstream Doppler,
   fUSi, ULM or WSS tests appropriate to the available data.

Independent velocity truth supports a downstream test. It does not by itself
prove preservation of every blood-IQ component. Controlled source/injection
experiments and complementary spatial/spectral diagnostics remain useful.

Inventory sources: [Waterloo carotid](https://github.com/open-h/Open-H-RF/blob/main/stats/data/files/waterloo-carotid.jsonl),
[Waterloo femoral vein](https://github.com/open-h/Open-H-RF/blob/main/stats/data/files/waterloo-femoralvein.jsonl),
[Twente](https://github.com/open-h/Open-H-RF/blob/main/stats/data/files/twente-vortexflow.jsonl),
[ResolveStroke](https://github.com/open-h/Open-H-RF/blob/main/stats/data/files/resolvestroke.jsonl).
The upstream data are CC BY 4.0; retain contributor attribution for derived figures.

## Interpretation of the two contrast references

Riemer et al. (2022) study native-blood versus contrast-enhanced velocimetry and
investigate scatter density, radial intensity and aggregation effects. It
supports controlling the scatterer regime. Their results also show useful
agreement between contrast-free and enhanced measurements, so contrast does not
imply a universal large WSS bias.

Pourjafar-Chelikdani et al. (2025) numerically analyze acoustically forced
encapsulated bubbles in elastic channels filled with shear-thinning fluid.
Bubble-induced stresses are a distinct physical mechanism. Their magnitude in
a particular Doppler acquisition depends on pressure/frequency, bubble size,
wall distance and pulse sequence. The non-contrast exclusion is a sensible
benchmark design choice; the paper alone does not quantify bias for every ULM
or large-artery WSS protocol.

Sources: [Riemer](https://pmc.ncbi.nlm.nih.gov/articles/PMC8843088/),
[bubble dynamics paper](https://doi.org/10.1016/j.jnnfm.2025.105514).
