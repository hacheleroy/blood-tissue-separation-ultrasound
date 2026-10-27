"""Preacquisition plan, strict input check, and paired physical information scorer.

No physical acquisition is bundled with this source. Source-free preprocessing
reuses the frozen numerical feature; optical motion and pump state are never
passed to registration or feature computation.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_dilation

from transport_lag_information import align, covariance

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/transport_lag_physical_v1.json"
FEATURE = ROOT / "scripts/transport_lag_information.py"
ROLES = ("dense_off", "dense_on", "burst_off", "burst_on")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sources():
    return {"config_sha256": sha(CONFIG), "analysis_sha256": sha(__file__),
            "feature_sha256": sha(FEATURE)}


def scaffold():
    c = json.loads(CONFIG.read_text())
    rng = np.random.default_rng(1907)
    sessions = []
    for session in c["sessions"]:
        blocks = []
        for motion in c["motion_conditions"]:
            for speed in c["speed_conditions"]:
                for repeat in range(1, c["repeats_per_cell"] + 1):
                    order = list(rng.permutation(ROLES))
                    blocks.append({"motion": motion, "speed": speed,
                                   "repeat": repeat, "capture_order": order,
                                   "captures": {
                                       role: {"iq_path": None, "iq_sha256": None,
                                              "motion_path": None, "motion_sha256": None,
                                              "flowmeter_mm_s": None,
                                              "capture_utc": None}
                                       for role in ROLES}})
        sessions.append({"id": session, "remount_id": None, "blocks": blocks})
    return {"schema_version": 1, "status": "INCOMPLETE: no physical measurements",
            **sources(), "scanner_id": None, "probe_id": None,
            "beamformer_sha256": None, "native_iq_unfiltered": None,
            "measured_tx_angles_deg": None, "receive_aperture_id": None,
            "pulse_id": None, "carrier_hz": None, "sound_speed_m_s": None,
            "contrast_agent": None, "flow_off_is_pump_off": None,
            "flow_reference": None, "acquisition_card_utc": None,
            "first_capture_utc": None, "geometry_from_independent_structural_scan": None,
            "geometry_path": None, "geometry_sha256": None,
            "optical_to_iq_transform_sha256": None,
            "clock_synchronization_record_sha256": None,
            "sessions": sessions}


def expected_times(schedule, c):
    cycle = np.arange(c["n_cycles"]) / c["cycle_rate_hz"]
    if schedule == "burst":
        cycle = (np.arange(c["n_cycles"]) // c["cycles_per_burst"]
                 * c["burst_start_spacing_s"]
                 + np.arange(c["n_cycles"]) % c["cycles_per_burst"]
                 / c["cycle_rate_hz"])
    return cycle[:, None] + np.arange(c["n_views"])[None, :] / (
        c["cycle_rate_hz"] * c["n_views"])


def within(root, relative):
    if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
        raise ValueError("Path must be a nonempty relative path within the data root")
    target = (root / relative).resolve()
    if not target.is_relative_to(root.resolve()):
        raise ValueError("Path escapes the data root")
    return target


def checked_file(root, relative, digest):
    path = within(root, relative)
    if not path.is_file() or not isinstance(digest, str) or len(digest) != 64:
        raise ValueError(f"Missing file or SHA-256: {relative}")
    if sha(path) != digest:
        raise ValueError(f"SHA-256 mismatch: {relative}")
    return path


def geometry(path):
    with np.load(path, allow_pickle=False) as d:
        if not {"near_wall", "tissue_control", "lumen", "z_mm", "x_mm"} <= set(d.files):
            raise ValueError("Missing predeclared structural mask or axes")
        wall, tissue, lumen = (np.asarray(d[k]) for k in
                               ("near_wall", "tissue_control", "lumen"))
        z, x = d["z_mm"], d["x_mm"]
    if any(m.ndim != 2 or m.dtype != bool for m in (wall, tissue, lumen)):
        raise ValueError("Structural masks must be 2D Boolean arrays")
    if wall.shape != tissue.shape or tissue.shape != lumen.shape:
        raise ValueError("Structural mask shapes differ")
    if z.shape != (wall.shape[0],) or x.shape != (wall.shape[1],):
        raise ValueError("Structural grid axes do not match masks")
    if min(wall.sum(), tissue.sum()) < 8 or not np.all(wall <= lumen):
        raise ValueError("At least eight inner-wall and tissue-control pixels required")
    if np.any(binary_dilation(lumen, iterations=1) & tissue):
        raise ValueError("Tissue control needs a one-pixel guard outside the lumen")
    if not (np.all(np.diff(z) > 0) and np.all(np.diff(x) > 0)):
        raise ValueError("Grid axes must increase")
    return wall, tissue, lumen, z, x


def checked_capture(item, role, root, mask_shape, c, target_speed):
    iq_path = checked_file(root, item.get("iq_path"), item.get("iq_sha256"))
    motion_path = checked_file(root, item.get("motion_path"), item.get("motion_sha256"))
    with np.load(iq_path, allow_pickle=False) as d:
        if set(d.files) != {"iq", "event_times_s"}:
            raise ValueError("IQ file must contain only iq and event_times_s")
        iq = d["iq"]
        times = np.asarray(d["event_times_s"], float)
    if (iq.shape != (c["n_cycles"], c["n_views"], *mask_shape)
            or not np.issubdtype(iq.dtype, np.complexfloating)
            or not np.all(np.isfinite(iq))):
        raise ValueError("Native IQ must be finite complex cycle × view × z × x")
    if times.shape != (c["n_cycles"], c["n_views"]) or not np.all(np.isfinite(times)):
        raise ValueError("Missing per-emission clocks")
    expected = expected_times(role.split("_")[0], c)
    if np.max(abs((times - times[0, 0]) - expected)) > c["timing_tolerance_s"]:
        raise ValueError("Emission clocks do not match declared fixed-count schedule")
    if not np.all(np.diff(times.ravel()) > 0):
        raise ValueError("Emissions are not strictly ordered")
    with np.load(motion_path, allow_pickle=False) as d:
        if not {"time_s", "displacement_mm", "uncertainty_mm"} <= set(d.files):
            raise ValueError("Independent motion record lacks clock/displacement/uncertainty")
        mt = np.asarray(d["time_s"], float)
        field = np.asarray(d["displacement_mm"], float)
        uncertainty = np.asarray(d["uncertainty_mm"], float)
    if (mt.shape != (c["n_cycles"],) or
            field.shape != (c["n_cycles"], *mask_shape, 2) or
            uncertainty.shape != field.shape or
            not np.all(np.isfinite(field)) or
            not np.all(np.isfinite(uncertainty)) or
            np.any(uncertainty < 0) or
            not np.all(np.isfinite(mt)) or
            np.max(abs(mt - times[:, 0])) > c["timing_tolerance_s"]):
        raise ValueError("Independent motion field/uncertainty is incomplete or unsynchronized")
    speed = item.get("flowmeter_mm_s")
    if not isinstance(speed, (int, float)) or not np.isfinite(speed):
        raise ValueError("Missing independent mean flowmeter speed")
    if role.endswith("_off"):
        if abs(speed) > 0.2 * c["speed_conditions"]["slow"]:
            raise ValueError("Pump-off flowmeter is not near zero")
    elif abs(speed - target_speed) > c["flow_speed_tolerance_fraction"] * target_speed:
        raise ValueError("Flowmeter outside predeclared target tolerance")
    if not isinstance(item.get("capture_utc"), str) or not item["capture_utc"]:
        raise ValueError("Missing capture clock metadata")
    return iq, times, field, uncertainty


def validate(manifest, root, session, c, load_data=False):
    errors = []
    if manifest.get("schema_version") != 1:
        errors.append("Unsupported manifest schema")
    for key, expected in sources().items():
        if manifest.get(key) != expected:
            errors.append(f"Frozen source mismatch: {key}")
    required = ["scanner_id", "probe_id", "beamformer_sha256", "flow_reference",
                "receive_aperture_id", "pulse_id", "carrier_hz", "sound_speed_m_s",
                "acquisition_card_utc", "first_capture_utc",
                "geometry_sha256", "optical_to_iq_transform_sha256",
                "clock_synchronization_record_sha256"]
    for key in required:
        if not manifest.get(key):
            errors.append(f"Missing independent acquisition field: {key}")
    if manifest.get("native_iq_unfiltered") is not True:
        errors.append("Unfiltered native complex IQ is required")
    angles = np.asarray(manifest.get("measured_tx_angles_deg", []), float)
    if (angles.shape != (c["n_views"],) or
            not np.all(np.isfinite(angles)) or
            np.max(abs(angles - c["tx_angles_deg"])) > c["tx_angle_tolerance_deg"]):
        errors.append("Three measured true views must match the predeclared angles")
    if manifest.get("contrast_agent") != "none":
        errors.append("Primary study requires bubble-free fluid; CEUS belongs in a separate analysis")
    if manifest.get("flow_off_is_pump_off") is not True:
        errors.append("Flow-off must be a measured pump-off capture")
    if manifest.get("geometry_from_independent_structural_scan") is not True:
        errors.append("Structure/masks must precede flow captures")
    if manifest.get("acquisition_card_utc") and manifest.get("first_capture_utc"):
        if manifest["acquisition_card_utc"] >= manifest["first_capture_utc"]:
            errors.append("Card is not dated before first capture")
    masks = None
    try:
        masks = geometry(checked_file(root, manifest.get("geometry_path"),
                                      manifest.get("geometry_sha256")))
    except (ValueError, OSError, KeyError) as e:
        errors.append(f"Geometry: {e}")
    sessions = manifest.get("sessions", [])
    found = next((s for s in sessions if s.get("id") == session), None)
    if found is None:
        errors.append(f"Missing session: {session}")
    elif not found.get("remount_id"):
        errors.append(f"Missing session remount ID: {session}")
    if len(sessions) != 2 or {s.get("id") for s in sessions} != set(c["sessions"]):
        errors.append("Development and fresh-remount evaluation sessions required")
    elif all(s.get("remount_id") for s in sessions):
        if sessions[0]["remount_id"] == sessions[1]["remount_id"]:
            errors.append("Evaluation must have a distinct remount ID")
    recorded_iq = [item.get("iq_sha256") for s in sessions
                   for b in s.get("blocks", [])
                   for item in b.get("captures", {}).values()
                   if isinstance(item, dict) and item.get("iq_sha256")]
    if len(recorded_iq) != len(set(recorded_iq)):
        errors.append("A native IQ file/hash is reused across captures or sessions")
    data = []
    if found:
        blocks = found.get("blocks", [])
        expected_blocks = {(m, s, r) for m in c["motion_conditions"]
                           for s in c["speed_conditions"]
                           for r in range(1, c["repeats_per_cell"] + 1)}
        keys = [(b.get("motion"), b.get("speed"), b.get("repeat")) for b in blocks]
        if len(keys) != len(expected_blocks) or set(keys) != expected_blocks:
            errors.append("Incomplete or duplicated motion/speed/repeat blocks")
        for b in blocks:
            if set(b.get("capture_order", [])) != set(ROLES) or len(b.get("capture_order", [])) != 4:
                errors.append(f"Invalid balanced capture order: {keys[blocks.index(b)]}")
            if set(b.get("captures", {})) != set(ROLES):
                errors.append(f"Incomplete four-capture block: {keys[blocks.index(b)]}")
                continue
            all_paths = [item.get("iq_path") for item in b["captures"].values()]
            if len(set(all_paths)) != 4:
                errors.append(f"Four roles need distinct native IQ files: {keys[blocks.index(b)]}")
            for role, item in b["captures"].items():
                if masks is None:
                    continue
                try:
                    payload = checked_capture(item, role, root, masks[0].shape,
                                              c, c["speed_conditions"][b["speed"]])
                    if load_data:
                        data.append((b["motion"], b["speed"], b["repeat"], role, payload))
                except (ValueError, OSError, KeyError) as e:
                    errors.append(f"{session}/{b.get('motion')}/{b.get('speed')}/"
                                  f"{b.get('repeat')}/{role}: {e}")
    return errors, masks, data


def feature(iq, times, c):
    mix = np.asarray(iq).mean(axis=1)
    registered, shifts, phases = align(mix, {"mixture": mix})
    x = registered["mixture"]
    cycle_times = times[:, 0]
    short, long = c["feature_short_lag"], c["feature_long_lag"]
    close = np.diff(cycle_times) < 2 / c["cycle_rate_hz"]
    short_c = covariance(x, short, pair_mask=close)
    long_c = covariance(x, long)
    return (abs(short_c) - long_c.real,
            {"max_estimated_translation_pixels": float(np.max(np.linalg.norm(shifts, axis=1))),
             "phase_span_rad": float(np.ptp(np.unwrap(np.angle(phases))))})


def paired_bootstrap(values, seed, n):
    values = np.asarray(values, float)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(values), (n, len(values)))
    estimates = np.mean(values[draws], axis=1)
    return {"mean": float(np.mean(values)),
            "lower_95": float(np.quantile(estimates, .025)),
            "upper_95": float(np.quantile(estimates, .975)),
            "positive_fraction": float(np.mean(values > 0))}


def score(manifest, root, session, c, lock=None):
    errors, masks, _ = validate(manifest, root, session, c)
    if errors:
        raise ValueError("Physical preflight failed: " + "; ".join(errors[:8]))
    if session == "evaluation":
        if lock is None:
            raise ValueError("Independent evaluation requires a development lock")
        for key, digest in sources().items():
            if lock.get(key) != digest:
                raise ValueError(f"Evaluation lock/source mismatch: {key}")
        for key in ("geometry_sha256", "beamformer_sha256"):
            if lock.get(key) != manifest[key]:
                raise ValueError(f"Evaluation lock/measurement mismatch: {key}")
        checked_file(root, lock.get("development_report_path"),
                     lock.get("development_report_sha256"))
    wall, tissue = masks[:2]
    captures = []
    blocks = []
    physical_session = next(s for s in manifest["sessions"] if s["id"] == session)
    for block in physical_session["blocks"]:
        motion, speed, repeat = (block[k] for k in ("motion", "speed", "repeat"))
        maps = {}
        for role in ROLES:
            iq, times, field, uncertainty = checked_capture(
                block["captures"][role], role, root, wall.shape, c, c["speed_conditions"][speed])
            f, registration = feature(iq, times, c)
            maps[role] = f
            captures.append({"motion": motion, "speed": speed, "repeat": repeat,
                             "role": role, "near_wall_median": float(np.median(f[wall])),
                             "tissue_control_median": float(np.median(f[tissue])),
                             "independent_motion_rms_mm": float(np.sqrt(np.mean(field**2))),
                             "independent_motion_uncertainty_rms_mm": float(np.sqrt(np.mean(uncertainty**2))),
                             "registration": registration})
            del iq, field, uncertainty
        row = {"motion": motion, "speed": speed, "repeat": repeat}
        for schedule in ("dense", "burst"):
            on, off = maps[schedule + "_on"], maps[schedule + "_off"]
            row[schedule + "_wall_on_off"] = float(np.median((on - off)[wall]))
            row[schedule + "_tissue_on_off"] = float(np.median((on - off)[tissue]))
            row[schedule + "_specific"] = (row[schedule + "_wall_on_off"]
                                           - row[schedule + "_tissue_on_off"])
        row["burst_advantage"] = row["burst_specific"] - row["dense_specific"]
        blocks.append(row)
    cells = {}
    for motion in c["motion_conditions"]:
        for speed in c["speed_conditions"]:
            key = motion + "/" + speed
            subset = [b for b in blocks if b["motion"] == motion and b["speed"] == speed]
            cells[key] = {
                metric: paired_bootstrap([b[metric] for b in subset],
                                         c["bootstrap_seed"], c["bootstrap_resamples"])
                for metric in ("dense_specific", "burst_specific", "burst_advantage",
                               "burst_tissue_on_off")}
    gate = all(cell["burst_specific"]["lower_95"] > 0 and
               cell["burst_advantage"]["lower_95"] > 0 for cell in cells.values())
    return {"status": "physical_information_scored",
            "interpretation": "Feature separability only; no filtered blood preservation, tissue leakage or velocity-profile claim",
            "session": session, **sources(),
            "geometry_sha256": manifest["geometry_sha256"],
            "beamformer_sha256": manifest["beamformer_sha256"],
            "n_captures": len(captures), "transmissions_per_capture": 240,
            "captures": captures, "paired_blocks": blocks, "cells": cells,
            "information_gate_passed": gate if session == "evaluation" else None}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="action", required=True)
    sc = sub.add_parser("scaffold")
    sc.add_argument("--out", type=Path, required=True)
    for action in ("preflight", "score"):
        x = sub.add_parser(action)
        x.add_argument("--manifest", type=Path, required=True)
        x.add_argument("--data-root", type=Path, required=True)
        x.add_argument("--session", choices=("development", "evaluation"), required=True)
        if action == "score":
            x.add_argument("--out", type=Path, required=True)
            x.add_argument("--lock", type=Path)
    lk = sub.add_parser("lock")
    lk.add_argument("--manifest", type=Path, required=True)
    lk.add_argument("--development-report", type=Path, required=True)
    lk.add_argument("--report-relative-to-data-root", required=True)
    lk.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    c = json.loads(CONFIG.read_text())
    if args.action == "scaffold":
        payload = scaffold()
    else:
        manifest = json.loads(args.manifest.read_text())
        if args.action == "preflight":
            errors, _, _ = validate(manifest, args.data_root, args.session, c)
            payload = {"physical_eligible": not errors, "session": args.session,
                       "issues": errors[:24], "total_issues": len(errors)}
        elif args.action == "score":
            lock = json.loads(args.lock.read_text()) if args.lock else None
            payload = score(manifest, args.data_root, args.session, c, lock)
        else:
            dev = json.loads(args.development_report.read_text())
            if dev.get("status") != "physical_information_scored" or dev.get("session") != "development":
                raise ValueError("Lock requires a completed physical development report")
            for k, digest in sources().items():
                if dev.get(k) != digest or manifest.get(k) != digest:
                    raise ValueError(f"Development source mismatch: {k}")
            payload = {**sources(), "geometry_sha256": manifest["geometry_sha256"],
                       "beamformer_sha256": manifest["beamformer_sha256"],
                       "development_report_path": args.report_relative_to_data_root,
                       "development_report_sha256": sha(args.development_report),
                       "status": "development locked before evaluation"}
    if args.action == "preflight":
        print(json.dumps(payload, indent=2))
    else:
        if args.out.exists():
            raise FileExistsError(args.out)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
        print(args.out)


if __name__ == "__main__":
    main()
