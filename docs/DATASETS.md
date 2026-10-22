# Open H-RF suitability screening

Screened on 22 October 2026. Evidence: consortium data cards plus the public
per-file metadata inventories; no experimental benchmark has yet been run.
GitHub tree inspected: `5eae31257f55391129c5aadb178659948448da5d`.
Catalog snapshot reports revision `88730aa47007553da09ffea61b499002e8b552b1`.
The live HDF5 range-header request timed out; inventory-derived timing still
needs confirmation against the actual file before quantitative reconstruction.

## Recommended starting point

**Waterloo carotid, `waterloo-carotid/data/Acq0.hdf5`.** The card describes human
carotid plane-wave RF with processed vector velocities. There are 93 acquisitions
from eight subjects. Acq0 is a longitudinal common-carotid baseline recording:
30,000 frames, two transmit angles spanning 20 degrees, 128 receive elements,
5 MHz Tx center frequency, 32.5 MHz sampling. The inventory reports 3 kHz frame
rate and 6 kHz pulse rate, giving about 10 seconds. Its decoded RF is about
62.9 GB; use windowed access. Two angles support an initial angle-aware study.
Known issues include receive multiplexing in some US4R-lite files, noise and
occasional aliasing. The card describes no contrast injection; native-blood
status is a strong inference that needs explicit confirmation. Stored velocities
follow a wall-filtered processing chain and cannot be used as oracle truth.

Sources: [data card](https://huggingface.co/datasets/nvidia/Open-H-RF/blob/main/waterloo-carotid/README.md),
[consortium file inventory](https://github.com/open-h/Open-H-RF/blob/main/stats/data/files/waterloo-carotid.jsonl).

## Secondary and excluded candidates

| Dataset | Evidence and role | Decision |
|---|---|---|
| Waterloo femoral vein | 15 acquisitions/subjects; 12,000 frames each; two plane waves, 192 elements; slow venous flow under contraction/tilt. Native-blood imaging is inferred from the description; explicit injection status unconfirmed. | Strong slow-flow/motion stress-test candidate; changing geometry complicates WSS. |
| TU/e carotid 2023 | 150 frames with 128 focused events interleaved with 21 plane waves. | A secondary imaging/motion source. Do not assume the plane-wave events are a uniform high-PRF Doppler ensemble. |
| Oslo carotid | Focused imaging with single/few frames. | Insufficient temporal ensembles for this comparison. |
| Resolve Stroke flow phantom | Matrix-array flow-on/off acquisition explicitly using microbubbles. | Exclude from the native-blood benchmark. |
| Twente vortex-flow phantom | Optical reference is attractive, but the card explicitly states acoustical scattering is enhanced with microbubbles. | Exclude from the non-contrast core; possible separate contrast/control experiment. |
| ULMShare | Microbubble transcranial ULM. | Exclude from the non-contrast core. |

Femoral Acq0 inventory: 12,000 frames at approximately 2,850 Hz, about 4.2 s;
other acquisitions have different frame rates/sample counts. Revisit timing must
be obtained per file. Stored dealiased vector fields are algorithm-derived and
can be noise-sensitive. Split all Waterloo data by subject across contributions
because participants can overlap.

Sources: [femoral card](https://github.com/open-h/Open-H-RF/blob/main/datasets/waterloo-femoralvein/README.md),
[femoral inventory](https://github.com/open-h/Open-H-RF/blob/main/stats/data/files/waterloo-femoralvein.jsonl),
[TU/e](https://huggingface.co/datasets/nvidia/Open-H-RF/blob/main/tue-carotid/README.md),
[Oslo](https://github.com/open-h/Open-H-RF/blob/main/datasets/oslo/B_carotid/README.md),
[Resolve Stroke](https://github.com/open-h/Open-H-RF/blob/main/datasets/resolvestroke/phantom_flow/README.md),
[Twente](https://github.com/open-h/Open-H-RF/blob/main/datasets/twente-vortexflow/README.md),
[ULMShare](https://github.com/open-h/Open-H-RF/blob/main/datasets/ulmshare/README.md).

## What these data can establish

They can support native-blood feasibility, motion robustness, measurement
repeatability, near-wall signal visibility and consistency across methods.
The inspected cards do not establish independent ground-truth blood/tissue IQ,
wall velocity, wall position, viscosity or WSS. Agreement with supplied VFI
maps is a processed-reference comparison, and depends partly on their upstream
wall filter. It cannot certify the recovered near-wall gradient.

I did not identify a non-contrast Poiseuille phantom with independent near-wall
velocity/WSS truth among the screened catalog entries. This is a screening
result, not proof that no such data exist elsewhere or in undocumented files.

Recommended study order:

1. Statistical IQ algebra/implementation tests (provided), then realistic
   non-contrast RBC/scatterer simulations or a controlled blood-mimicking phantom.
2. Freeze the method/measurement configurations and independently verify flow,
   wall geometry, viscosity, phase sign and sample timing in that controlled test.
3. Beamform a small unfiltered Waterloo CCA window with preserved view data;
   compare B-mode and phase against the provider's reconstruction.
4. Run the shared filter grid on held-out subjects. Use in-vivo observations
   as feasibility/consistency evidence, with oracle arms disabled.

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
