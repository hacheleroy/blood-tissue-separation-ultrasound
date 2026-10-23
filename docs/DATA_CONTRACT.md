# Input and acquisition contract

## Unfiltered IQ

```python
from btbench import Acquisition
from btbench.io import save_acquisition

a = Acquisition(
    iq=complex_iq,              # (T,Z,X), BEFORE clutter filtering/envelope/log
    prf_hz=3000.0,              # revisit rate of this same IQ view
    f0_hz=5e6,
    sound_speed_mps=1540.0,
    phase_sign=-1,             # verify on a known moving reflector
    timestamps_s=frame_times,  # seconds; optional but recommended
    z_m=z_grid, x_m=x_grid,
    tissue_roi=tissue_mask,    # anatomical/calibration ROI, never velocity truth
    flow_mask=lumen_mask,      # used for evaluation, not data-only source labels
    metadata={"source": "acquisition identifier", "contrast_status": "confirmed_absent"}
)
save_acquisition("unfiltered_iq.npz", a)
```

NPZ storage uses arrays and JSON metadata; pickle loading is disabled.
Uniform timestamps must agree with `prf_hz`. Gaps/nonuniformity are rejected;
split into contiguous segments without joining independent ensembles.

`phase_sign=+1` means a positive true projected velocity increases IQ phase;
`-1` means it decreases phase. The RF reference DAS convention normally needs
`-1`, tested on a moving point reflector. Never infer sign from a plausible
looking vascular image. Store actual transmit and receive directions for VFI.

Optional arrays `blood_true`, `tissue_true`, `velocity_true_mps` are reserved for
independent simulation/physical truth. Vendor velocity maps belong in separate
comparison data. They are not true blood signal or valid oracle targets.

For WSS set `wall` with physical positions `(z,x)`, inward unit normals, strictly
fluid-side distances, viscosity, optional wall velocity, and optional analytic
WSS. Supply frame-specific geometry separately for moving walls; the simple
Acquisition currently stores a static geometry for one measurement interval.
Use the vector-Doppler API for known, observable multi-angle acquisition
geometries. The default scalar estimator supplies axial velocity, which by
itself generally cannot determine vessel-tangential WSS.

## Raw Open H-RF data

Waterloo raw data are real RF channel arrays. They must be beamformed back to
complex IQ. The stored `image/values` and `power_doppler/values` fields are in dB
and cannot substitute for IQ. `vector_velocity_*` fields are processed outputs.

Track roots vary with format version, e.g. `tracks/track_0`. Inspect the file
inventory first. Most Waterloo RF chunks contain one full frame with all
transmits/channels; select a short temporal window before doing any reconstruction.
Install `.[io]` to register HDF5 Blosc filters (`hdf5plugin`) for actual RF reads.

Prefer a provider/zea reconstruction that preserves IQ before envelope/log
processing and before the provider's existing wall filter. Their B-mode example
alone does not establish quantitative phase accuracy.

The CPU reference adapter `export_openhrf_rf_window` takes an explicit geometry JSON with:

- `frames_start_stop`: bounded contiguous range, e.g. [0,64].
- `transmit_indices`: selected stored plane-wave events.
- `root`: HDF5 track path, e.g. "tracks/track_0".
- `z_m`, `x_m`, `tx_angles_rad`, `fs_hz`, `f0_hz`, `sound_speed_mps`.
- `initial_times_s`: A/D start time for each event.
- `tx_reference_times_s`, optional `tx_origins_zx_m`: the verified far-field Tx
  wavefront time origin satisfying the equation in `beamform.py`.
- Optional `revisit_prf_hz`, `frame_times_s`: checked against recorded event intervals.
  The adapter derives same-view revisit timing and constant view offsets from
  `scan/time_to_next_transmit`; it rejects gaps/varying schedules.
- `phase_sign`, optional `f_number`, `max_raw_bytes`.
- `time_reference_verified=true`, `full_receive_aperture_verified=true`.

The last flags are gates, not automatic validation. Only set them after checking
the actual sequence and comparing with the provider reconstruction. Some US4R-lite
files multiplex two receive subapertures across separate events. Their assembly
and motion-dependent phase cannot be guessed from the element count. Avoid those
files for the first pilot. Compounding different transmit times also requires
an explicit motion/phase model. Preserve uncompounded views for Axis 3: Local ultrasound STAP.

## STAP API

```python
from btbench.stap import local_stap

result = local_stap(
    iq_views, prf_hz, f0_hz, tx_rad, rx_rad,
    tissue_velocities_zx_mps=estimated_tissue_velocity_grid,
    view_offsets_s=transmit_offsets,
    phase_sign=-1, taps=4,
    adaptive=True, training_roi=tissue_mask,
    guard=1, training_radius=3,
)
```

`iq_views` is `(T,A,Z,X)`. The optional runner arm `method="stap"` uses these
same arguments under `params` and requires `Acquisition.iq_views` plus a matching
`Acquisition.iq = iq_views.mean(axis=1)` baseline. The tissue velocity grid is
an independently measured/estimated nuisance parameter; truth-assisted steering
must be labeled explicitly as an oracle study. Significant tissue deformation,
strong aberration and unknown channel gains violate this initial model.


## Embedded multi-view acquisition

RF export embeds all complex views, with `iq = iq_views.mean(axis=1)`, and stores
`tx_angles_rad`, `rx_angles_rad`, `view_offsets_s` and `timestamps_s`. The loader
restores these arrays directly; a sidecar is no longer needed. The standard STAP
runner derives physical frequencies, geometry and absolute packet origin from
these fields and rejects mismatched declared settings.

Independent simulated view components can be supplied as `blood_views_true` and
`tissue_views_true`, together with their coherent means in `blood_true` and
`tissue_true`. Their dimensions and mean consistency are checked. Those fields
are evaluation-only. A scalar/compounded component cannot substitute for its
unknown angular distribution in STAP diagnostics.

The `write-config` CLI accepts `--prf-hz`, `--packet` and `--long-window`:

```bash
btbench write-config configs/my_3khz.json --prf-hz 3000 --packet 32 --long-window 96
```

The [bounded experimental pilot](EXPERIMENTAL_PILOT.md) preserves disjoint
calibration/evaluation packets and their clock origin. It removes every source
truth field, including any that may have accompanied an input simulation.
