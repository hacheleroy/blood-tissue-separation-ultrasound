# Paired physical transport-lag information test

27 October 2026. **Preacquisition protocol. No physical flow-on/off IQ has
been scored.** This study follows the [source-separated numerical information
test](TRANSPORT_LAG_INFORMATION_V1.md). The frozen
[configuration](../configs/transport_lag_physical_v1.json),
[scorer and input check](../scripts/transport_lag_physical.py), and
[unfilled inventory](../runs/transport_lag_physical_v1/manifest_template.json)
specify the measurement. The broader
[phase-10 phantom bridge](PHANTOM_PHASE10_BRIDGE.md) addresses filter transfer;
its 180 **unfilled** captures do not supply this timing comparison.

## Physical question

Does a longer observation reveal a reproducible blood-associated change in
the *actual* unfiltered IQ mixture after tissue-frame registration, without
making stationary-tissue fluctuations look like flow? Compare two acquisition
schedules with **80 cycles × 3 actual views = 240 transmissions per capture**:

| Schedule | Same-view cycle times | Span of first to last emission |
|:--|:--|--:|
| Dense | 80 consecutive cycles at 3700 cycles/s | about 21.5 ms |
| Burst | 10 groups of 8 cycles at the same within-group rate, groups 22 ms apart | about 200.1 ms |

Transmit waveform, view order (0°, −5°, +5°), receive aperture, gain, carrier,
beamformer, probe position, and output grid are fixed within a session. The
clock file records **all 240 actual emission times**, and the scorer rejects
a schedule more than 10 μs away from the specified relative timing. The
longer schedule changes elapsed time, not transmission count. A scanner that
cannot produce the three genuine views and both schedules needs an amended
protocol frozen *before* capture; one single diverging transmit must not be
called three independent angles.

Use a bubble-free blood-mimicking suspension in a structurally measured flow
tube and tissue-mimicking surround. Maintain stable temperature and record
pump/flowmeter data. For each motion condition (static, rigid, nonrigid) and
flow target (mean 3 and 20 mm/s), acquire **five matched four-capture blocks**:
dense pump-off, dense pump-on, burst pump-off, burst pump-on. A seeded order is
stored in the inventory before capture; allow washout/settling between pump
states. Record independent stage encoder/optical motion on the IQ grid and its
uncertainty at each cycle, including out-of-plane motion if measurable. Repeat
the 30 blocks in a **fresh-remount evaluation session**. Total: 240 captures
and 57,600 emissions; the experimental unit is a matched four-capture block,
not a correlated near-wall pixel.

An independently dated acquisition card must identify scanner/probe,
waveform/receive aperture, measured angles, clocks, TGC and beamformer, tube
structure, optical-to-IQ mapping and synchronization, remounts and flowmeter
calibration before the first capture. The inventory's claimed timestamps and
hashes can be checked for consistency; the software cannot authenticate lab
provenance. Keep native unfiltered complex IQ per view in
`iq: (80,3,z,x)` and emission times in `event_times_s: (80,3)`, in a
`.npz` with no object arrays. Also retain the original raw channel data
where available and its hash outside Git. The scorer expects each independent
optical/stage file to provide `time_s: (80,)`,
`displacement_mm: (80,z,x,2)` and equally shaped nonnegative
`uncertainty_mm`. Hash all capture and structure files in the inventory.

## Fixed regions and transform

Before any pump-on capture, fix the structural grid and Boolean masks from a
separate structural image: `lumen`, an **inner near-wall** subset of the
lumen, and a tissue-control region with at least a one-pixel guard *outside*
the lumen. Store increasing `z_mm` and `x_mm` axes; require at least eight
pixels in each scored region. Use the innermost-wall band that the measured
resolution actually supports, and document its width and partial-volume
uncertainty. The geometric band is a proxy for weak slow flow; it cannot be
renamed a measured low-velocity subset until independent flow reference
confirms that interpretation. Apply the same masks to all four captures.
No source-resolved power, bolus intensity, pump-on contrast or evaluation
session output may select pixels.

For each capture separately, compound the three native complex IQ views,
estimate a bounded B-mode translation to its first image on the whole field,
then estimate a single global complex-IQ phase. The independent optical
displacement is **recorded for evaluation only**, never supplied to this
transform. With registered IQ `I_t(p)`, compute normalized covariance
`C_l(p)` as defined in the [numerical test](TRANSPORT_LAG_INFORMATION_V1.md),
then the fixed, dimensionless feature

\[
S(p)=|C_1(p)|-\operatorname{Re} C_{48}(p).
\]

Short-lag pairs spanning burst gaps are excluded. The lag-48 delay is about
13 ms in dense mode and 132 ms in burst mode; comparing those distinct
physical intervals is the proposed information test, not an equal-lag
software comparison. No threshold, classifier, directional speed search,
source truth or calibration phase is fitted.

For each four-capture block and schedule, compute the **median pixelwise
pump-on minus pump-off feature** in the fixed near-wall band and separately
in tissue control. Subtract tissue-control response from near-wall response
to describe local specificity. The burst advantage is this value in burst
minus the same value in dense. A simple frame-by-frame subtraction between
different pump acquisitions would be invalid because their speckle is not
phase matched. Report every block, motion amplitude/uncertainty and
registration shift as well as the cell summaries.

Use a paired bootstrap over the five whole blocks (2000 resamples, seed 1907)
for each of the six motion × speed cells. The **preregistered information
gate** requires the untouched evaluation session to have a strictly
positive lower 95% bootstrap bound for both burst specificity and burst
advantage in **every** cell. The development session may reveal failures
and support future revisions, but a revision requires a new acquisition and
freeze before testing it in a fresh session. A positive gate would justify
developing a source-blind filter; it would **not** establish blood IQ
preservation, tissue-power suppression, post-filter blood/background ratio
or an accurate velocity profile. These require a later fixed filter evaluated
against independent flow measurements.

## Input and independence gates

The [scorer](../scripts/transport_lag_physical.py) rejects missing pump-off
captures, reused IQ filenames within a block, timing/count mismatches,
non-complex or post-filtered declared IQ, geometry overlap, missing motion
fields or uncertainty, flowmeter values outside ±20% of the stated on
target, pump-off above 0.6 mm/s, missing hashes, and repeated remount IDs.
An evaluation run additionally requires a lock recording the unchanged
analysis/configuration/feature hashes, structural geometry and beamformer,
and a completed development report hash. It cannot prove that a field
claimed to be an optical measurement really came from a camera.

The [Resolve Stroke public flow phantom](https://github.com/open-h/Open-H-RF/blob/main/datasets/resolvestroke/phantom_flow/README.md)
can supply a **separate CEUS retrospective check**, not this primary
measurement. Its five 4000-frame clips comprise a pre-microbubble baseline
and later wash-in/plateau intervals, use a single diverging-wave transmit,
and contain microbubble contrast. The tube mask is known geometry; the
provided radial-velocity map is a device-derived estimate from the full
bolus passage, rather than an independent local velocity measurement from
the released clips. It has no synchronized imposed tissue-motion record for
this protocol. Calling baseline versus wash-in a bubble-free matched
pump-off versus pump-on pair would answer a different question.

## Reproducible handoff

```bash
python scripts/transport_lag_physical.py scaffold \
  --out runs/transport_lag_physical_v1/manifest_template.json
python scripts/transport_lag_physical.py preflight \
  --manifest runs/transport_lag_physical_v1/manifest_template.json \
  --data-root /path/to/physical_capture --session development
```

The committed template intentionally fails preflight. Once all development
captures and metadata exist, score the filled inventory with
`score --manifest ... --data-root ... --session development --out ...`;
run `lock --manifest ... --development-report ... --report-relative-to-data-root
... --out ...` before seeing the evaluation set. The analogous evaluation
`score` requires `--lock`. Store raw measurements, filled inventory,
development lock and physical results with the acquisition laboratory; publish
only a consented compact result archive, not unreviewed raw captures.
