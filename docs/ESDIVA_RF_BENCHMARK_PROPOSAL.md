# eSDIva RF benchmark and Axis 3 literature update

24 October 2026. **Design assessment; no eSDIva benchmark has been run.** This
proposed acoustic-simulation tier sits between the existing statistical-IQ
studies and the independently acquired physical-phantom protocol
([v07](PHANTOM_PROTOCOL_V07.md)). It does not change v07's frozen acquisition
or its experimental validation gate.

## What Example 22 provides

Use [eSDIva Example 22](https://estebanrivera08.github.io/eSDIva/examples/example22_flow_doppler/)
at upstream commit
[`ab7eaef6d23428762ab5b9079ab03a5aa99af085`](https://github.com/EstebanRivera08/eSDIva/tree/ab7eaef6d23428762ab5b9079ab03a5aa99af085/examples/example22_flow_doppler)
(BSD-3-Clause). Its [phantom definition](https://github.com/EstebanRivera08/eSDIva/blob/ab7eaef6d23428762ab5b9079ab03a5aa99af085/examples/example22_flow_doppler/step1_define_flow_phantom.py)
returns scatterer positions for every emission, amplitudes, tissue/blood row
counts and blood velocities. Its [acquisition step](https://github.com/EstebanRivera08/eSDIva/blob/ab7eaef6d23428762ab5b9079ab03a5aa99af085/examples/example22_flow_doppler/step2_acquire_RF.py)
calls `sequence_rf` with those moving scatterers and checkpointed per-event
channel RF. The [processing step](https://github.com/EstebanRivera08/eSDIva/blob/ab7eaef6d23428762ab5b9079ab03a5aa99af085/examples/example22_flow_doppler/step3_doppler_processing.py)
beamforms, demodulates, filters and saves complex IQ. This is **a generator and
an example workflow**, not a published, ready-to-download bank of the separate
source IQ acquisitions required for this benchmark.

The default comparison reproduces the plane-wave-compounding question of
[Bercoff et al. (2011)](https://pubmed.ncbi.nlm.nih.gov/21244981/):

| Arm | Emissions | Saved slow-time frames | Use here |
|:--|--:|--:|:--|
| A: 72 focused lines × 8 | 576 | 8 per line | Check focused/plane-wave acquisition differences; not a matched input for the three-view runner. |
| B: 9 angles × 8 | 72 | 8 | Short stationary-flow smoke check only. |
| C: 9 angles × 32 | 288 | 32 | Single-stream SVD/wall-filter sanity check after compounding. |

The original phantom has a 3 MHz probe, nine angles from -20° to 20°, a
3,700 Hz **same-angle revisit rate**, a 33,300 Hz emission rate and a 0.45 m/s
peak flow **along a vessel** at 60° to the beam. The stationary tissue and nominal
0.2 blood-scatterer amplitude relative to tissue (about -14 dB in amplitude)
make its default setting much easier than a moving-tissue, weak/slow-flow test.
The number is a scatterer amplitude setting, **not** a measured per-voxel
blood-to-tissue power ratio. The existing Step 3 compounds all nine angles
within each frame and writes `B_iq` and `C_iq` as time × x × z; it does not
retain a view axis or blood/tissue-only IQ. Eight or 32 frames cannot supply
the shared runner's 64-frame warm-up/history, let alone disjoint calibration.

## Reproducible extension before comparing filters

1. First reproduce upstream Example 22 unchanged, logging its commit, package
   environment, per-event RF coordinates, seed and upstream output checksums.
   Treat this as simulator/beamformer acceptance, **not** an Axis 3 result.
2. Create a separate, declared acquisition with three selected plane-wave
   angles (for example -15°, 0°, +15°), one IQ volume **per emission and view**,
   at least 64 disjoint calibration frames and 256 subsequent frames (first 64
   for history warm-up, then 192 scored). Preserve
   the actual per-angle transmit order, `t0_per_event`, TX/RX geometry, emission
   timestamps and same-angle PRF. Beamform each emission separately with its
   matching event, use one calibration-derived demodulation frequency and
   fixed IQ phase convention, then crop a common 8 × 16 scoring tile. Verify
   view gain/phase and forward-model steering on independent tissue calibration.
   This is a **new acquisition**, not the nine-angle compounded upstream arm.
3. Repeat the *same positions and events* with tissue amplitudes set to zero,
   and with blood amplitudes set to zero. In eSDIva's
   [`sequence_rf`](https://github.com/EstebanRivera08/eSDIva/blob/ab7eaef6d23428762ab5b9079ab03a5aa99af085/src/esdiva/reception/base.py),
   RF is linear in scatterer amplitudes, and the original positions can be
   retained even for zero-amplitude scatterers to hold the time grid fixed.
   Check numerically that combined RF equals blood RF + tissue RF before
   injecting one fixed receiver-noise realization. Apply identical beamforming
   and fixed demodulation to all three. Save the raw components and a manifest,
   but expose only mixtures, geometry and independently obtained calibration
   to the filter API; components and true motion/velocity are scoring-only.
4. Declare stationary, rigid and spatially varying moving-tissue scenes;
   flow-off and tissue-only controls; strong-clutter and slower-flow settings;
   and several development seeds. The upstream trajectory is stationary-tissue
   fast flow and **does not itself test motion compensation**. Revisit rate,
   duration and frequency bands must be chosen for each velocity regime before
   a new held-out seed bank. At 3 MHz and 3,700 Hz, 3 mm/s axial motion gives
   only about 0.0032 Doppler cycles per frame (about half this at 60°). A
   16-frame snapshot alone cannot separate it from a nearly coincident tissue
   manifold. Report these near-overlap cases as failures/inconclusive when the
   signal is not identifiable; use a slower revisit rate and/or longer support
   as an explicitly different acquisition, subject to aliasing at fast flow.
5. Run the existing three-axis runner and global/local SVD, Baranger SSM and
   Stephen Lee AER baselines on **identical uncompounded input, calibration,
   training split and score frames**. Retune only on declared development
   scenes, freeze geometry, reconstruction, algorithm settings and score rules,
   then evaluate independent held-out seeds. Score complex blood-IQ error,
   isolated blood response, tissue-only leakage, output/noise amplification,
   blood preservation versus leakage by motion/velocity/clutter stratum,
   calibration error, displacement error and computation. A refitted
   mixture-only operator and an operator frozen on a mixture answer different
   questions: retain both diagnostics and never rank nonlinear/adaptive methods
   solely by a frozen component decomposition. Velocity/PD are secondary.

Before full simulation, run the upstream runtime estimator on the local machine.
Full-length, separately simulated components multiply RF work and disk use.
The extension should check RF additivity, phase continuity, zero-flow false
positives, interpolation/roundtrip controls and component-to-mixture processing
parity on a tiny pilot **before** spending the full acquisition budget.

## Literature map for Axis 3: Local ultrasound STAP

The current [Axis 3 prototype](STAP_BENCHMARK.md) combines *transmit-angle*
and slow-time beamformed IQ, uses a tissue steering null and optionally loaded
local covariance. The LCMV velocity bank is a **detector**, not a reconstructed
blood-IQ filter. Its prototype does not use a full receive-channel aperture or
provide a distortionless response for every possible blood velocity. Classical
[radar STAP](https://doi.org/10.1109/MAES.2004.1263229) suggests the joint
space-time covariance/steering formalism; its array geometry and clutter
training assumptions cannot simply be transplanted to plane-wave ultrasound.

| Primary source | Evidence relevant to our experiment | Action / boundary |
|:--|:--|:--|
| [Kargel et al. 2003](https://pmc.ncbi.nlm.nih.gov/articles/PMC2715560/) | Adaptive ultrasound eigenfilter for spatially/temporally varying tissue in stationary and moving flow phantoms. | Compare tissue-only guarded local training and calibration mismatch; do not infer independent snapshots from overlapping pixels. |
| [Bjærum et al. 2002](https://pubmed.ncbi.nlm.nih.gov/12075963/); [Tierney et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5605932/) | Motion-adapted phase increments and adaptive amplitude/frequency demodulation reduce tissue clutter spread. | Compare motion-demodulated controls with the **same** calibration support as STAP before attributing gains to covariance adaptation. |
| [Song et al. 2017](https://pubmed.ncbi.nlm.nih.gov/27608455/); [Huang et al. 2024](https://arxiv.org/abs/2405.11105) | Overlapping local SVD and Gaussian-windowed local EVD handle heterogeneous tissue/noise; locality and threshold selection matter. | Match spatial support; log covariance sample count/effective rank, loading, guard width, leakage and runtime. These methods are ultrasound comparators, not identical STAP implementations. |
| [Huang et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38579486/) | Angular/temporal/spatial HOSVD for coronary flow under tissue motion; fewer angles can trade angular diversity for more time samples. | Add an angular HOSVD comparator if its tensor and classifier are reproduced; preserve views and emission budget. |
| [Ozgun and Byram 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8345228/) | Delayed receive-channel aperture × space × time HOSVD, including separately simulated blood/tissue and moving tissue at -40 dB in Field II. | eSDIva channel RF enables a **separate aperture-domain track** after delay calibration. Do not describe transmit angles as receive elements or claim current Axis 3 reproduces this method. |
| [Lee et al. 2026](https://arxiv.org/html/2610.23099v1) | AER continuously remaps covariance eigenvalues rather than hard truncating; published outcome is transcranial *microbubble* ULM. | Keep temporal/spatial AER as explicit same-input comparator. Its microbubble improvement is no prediction for native RBC flow or absolute blood-IQ fidelity; inspect regularization and noise gain. |

**First STAP upgrade to test:** derive tissue angle-time steering from a
disjoint tissue calibration with measured view phases and nonconstant tissue
motion; compare the orthogonal null to loaded/local covariance on the same
samples and to the calibration-demodulated SVD control. Reject poorly supported
or ill-conditioned fits, and report blood attenuation even when tissue leakage
falls. Scan the fixed blood-velocity LCMV detector separately with explicit
constraint overlap and white-noise gain. Near-identical tissue and blood
steering cannot be both nulled and preserved, regardless of loading.

The intended result is an **acoustic-simulation validation of separation**, not
an experimental phantom result or a claim that the upstream Doppler/velocity
plots already validate our filters. The physical v07 protocol remains pending
until independent data exist.
