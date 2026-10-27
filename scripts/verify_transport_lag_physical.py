"""Verify the preacquisition lease and arithmetic with synthetic IQ only."""
import copy
import importlib.util
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "transport_lag_physical", ROOT / "scripts/transport_lag_physical.py")
import sys
sys.path.insert(0, str(ROOT / "scripts"))
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)


def check():
    c = json.loads(study.CONFIG.read_text())
    template = json.loads((ROOT / "runs/transport_lag_physical_v1/manifest_template.json").read_text())
    freeze = json.loads((ROOT / "runs/transport_lag_physical_v1/freeze.json").read_text())
    assert template == study.scaffold()
    assert template["status"].startswith("INCOMPLETE")
    for name, digest in study.sources().items():
        assert freeze[name] == digest and template[name] == digest
    assert freeze["template_sha256"] == study.sha(
        ROOT / "runs/transport_lag_physical_v1/manifest_template.json")
    assert freeze["protocol_sha256"] == study.sha(
        ROOT / "docs/TRANSPORT_LAG_PHYSICAL_PROTOCOL_V1.md")
    assert c["n_cycles"] * c["n_views"] == c["transmissions_per_capture"] == 240
    assert len(template["sessions"]) == 2
    assert all(len(s["blocks"]) == 30 for s in template["sessions"])
    assert all(len(b["captures"]) == 4 for s in template["sessions"] for b in s["blocks"])
    with TemporaryDirectory() as temp:
        root = Path(temp)
        errors, _, _ = study.validate(template, root, "development", c)
        assert errors and any("Geometry" in e for e in errors)
        wall = np.zeros((12, 12), bool)
        wall[5, 2:10] = True
        lumen = wall.copy()
        lumen[6, 2:10] = True
        tissue = np.zeros_like(wall)
        tissue[0, 2:10] = True
        geometry = root / "geometry.npz"
        np.savez(geometry, near_wall=wall, lumen=lumen, tissue_control=tissue,
                 z_mm=np.arange(12, dtype=float), x_mm=np.arange(12, dtype=float))
        case = copy.deepcopy(template)
        case.update(status="SYNTHETIC FIXTURE; NOT PHYSICAL", scanner_id="fixture",
                    probe_id="fixture", beamformer_sha256="f" * 64,
                    native_iq_unfiltered=True, measured_tx_angles_deg=c["tx_angles_deg"],
                    receive_aperture_id="fixture", pulse_id="fixture",
                    carrier_hz=2.5e6, sound_speed_m_s=1540,
                    contrast_agent="none", flow_off_is_pump_off=True,
                    flow_reference="synthetic", acquisition_card_utc="2026-01-01T00:00:00Z",
                    first_capture_utc="2026-01-02T00:00:00Z",
                    geometry_from_independent_structural_scan=True,
                    geometry_path="geometry.npz", geometry_sha256=study.sha(geometry),
                    optical_to_iq_transform_sha256="e" * 64,
                    clock_synchronization_record_sha256="d" * 64)
        short = copy.deepcopy(c)
        short["motion_conditions"] = ["static"]
        short["speed_conditions"] = {"slow": 3.0}
        short["repeats_per_cell"] = 1
        for s in case["sessions"]:
            s["remount_id"] = s["id"] + "_sim"
            s["blocks"] = [s["blocks"][0]]
        rng = np.random.default_rng(30)
        base = rng.normal(size=(12, 12)) + 1j * rng.normal(size=(12, 12))
        for role in study.ROLES:
            schedule = role.split("_")[0]
            times = study.expected_times(schedule, short)
            iq = np.broadcast_to(base, (80, 3, 12, 12)).astype(np.complex64).copy()
            if role.endswith("_on"):
                iq[:, :, wall] += 1.2 * np.exp(1j * 12 * times[:, 0])[:, None, None]
            noise = (rng.normal(size=iq.shape) + 1j * rng.normal(size=iq.shape)) * 0.01
            iq += noise.astype(np.complex64)
            iq_path = root / (role + ".npz")
            np.savez(iq_path, iq=iq, event_times_s=times)
            motion_path = root / (role + "_motion.npz")
            np.savez(motion_path, time_s=times[:, 0],
                     displacement_mm=np.zeros((80, 12, 12, 2)),
                     uncertainty_mm=np.full((80, 12, 12, 2), .01))
            item = case["sessions"][0]["blocks"][0]["captures"][role]
            item.update(iq_path=iq_path.name, iq_sha256=study.sha(iq_path),
                        motion_path=motion_path.name, motion_sha256=study.sha(motion_path),
                        flowmeter_mm_s=3.0 if role.endswith("_on") else 0.,
                        capture_utc="2026-01-02T00:00:01Z")
        errors, _, _ = study.validate(case, root, "development", short)
        assert not errors, errors
        result = study.score(case, root, "development", short)
        assert result["n_captures"] == 4 and result["paired_blocks"][0]["burst_advantage"] > 0.1
        repeated = copy.deepcopy(case)
        repeated["sessions"][1]["blocks"][0]["captures"] = copy.deepcopy(
            case["sessions"][0]["blocks"][0]["captures"])
        errors, _, _ = study.validate(repeated, root, "evaluation", short)
        assert any("reused" in e for e in errors)
        for role in study.ROLES:
            source = root / (role + ".npz")
            with np.load(source) as d:
                eval_iq = d["iq"].copy() + (1e-4 + 1e-4j)
                eval_times = d["event_times_s"]
            target = root / ("evaluation_" + role + ".npz")
            np.savez(target, iq=eval_iq, event_times_s=eval_times)
            repeated["sessions"][1]["blocks"][0]["captures"][role].update(
                iq_path=target.name, iq_sha256=study.sha(target))
        errors, _, _ = study.validate(repeated, root, "evaluation", short)
        assert not errors, errors
        try:
            study.score(repeated, root, "evaluation", short)
        except ValueError as e:
            assert "lock" in str(e)
        else:
            raise AssertionError("Evaluation without a lock succeeded")
        dev_report = root / "dev.json"
        dev_report.write_text(json.dumps(result))
        lock = {**study.sources(), "geometry_sha256": repeated["geometry_sha256"],
                "beamformer_sha256": repeated["beamformer_sha256"],
                "development_report_path": "dev.json",
                "development_report_sha256": study.sha(dev_report)}
        held = study.score(repeated, root, "evaluation", short, lock)
        assert held["n_captures"] == 4 and isinstance(held["information_gate_passed"], bool)
        lock["feature_sha256"] = "0" * 64
        try:
            study.score(repeated, root, "evaluation", short, lock)
        except ValueError as e:
            assert "lock/source" in str(e)
        else:
            raise AssertionError("Modified frozen source passed evaluation")
        dense_item = case["sessions"][0]["blocks"][0]["captures"]["dense_on"]
        bad = dict(dense_item, iq_sha256="0" * 64)
        try:
            study.checked_capture(bad, "dense_on", root, wall.shape, short, 3.0)
        except ValueError as e:
            assert "SHA-256" in str(e)
        else:
            raise AssertionError("Changed capture bytes passed")
        wrong_times = study.expected_times("burst", short).copy()
        wrong_times[9:] += 0.001
        wrong_path = root / "wrong_times.npz"
        np.savez(wrong_path, iq=iq, event_times_s=wrong_times)
        bad = dict(dense_item, iq_path=wrong_path.name, iq_sha256=study.sha(wrong_path))
        try:
            study.checked_capture(bad, "dense_on", root, wall.shape, short, 3.0)
        except ValueError as e:
            assert "clocks" in str(e)
        else:
            raise AssertionError("Wrong schedule passed")
    print("Transport-lag physical protocol: frozen and synthetic arithmetic/negative controls verified")


if __name__ == "__main__":
    check()
