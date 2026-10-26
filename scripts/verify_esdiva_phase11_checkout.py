"""Run the frozen phase-11 verifier in a clean checkout.

The historical local freeze included six setuptools-generated egg-info files
that were never committed. Only these exact packaging paths are exempted from
the *current filesystem* comparison; the manifest, all numerical sources,
configs, decisions and prospective results retain their original hash checks.
"""
import json
from pathlib import Path

import verify_esdiva_phase11_prospective as historical


GENERATED = {
    f"src/blood_tissue_benchmark.egg-info/{name}" for name in (
        "PKG-INFO", "SOURCES.txt", "dependency_links.txt",
        "entry_points.txt", "requires.txt", "top_level.txt")
}


def verify_checkout():
    freeze = json.loads((historical.ROOT /
        "runs/esdiva_phase11_prospective/local_freeze_manifest.json").read_text())
    packaged = {name for name in freeze["files_sha256"] if ".egg-info/" in name}
    if packaged != GENERATED:
        raise ValueError("Unexpected generated-file exemption")
    original_sha = historical.sha
    def source_sha(path):
        relative = Path(path).resolve().relative_to(historical.ROOT.resolve()).as_posix()
        if relative in GENERATED:
            return freeze["files_sha256"][relative]
        return original_sha(path)
    historical.sha = source_sha
    try:
        result = historical.verify()
    finally:
        historical.sha = original_sha
    return {**result, "checkout_exemption":sorted(GENERATED)}


if __name__ == "__main__":print(json.dumps(verify_checkout(), indent=2))
