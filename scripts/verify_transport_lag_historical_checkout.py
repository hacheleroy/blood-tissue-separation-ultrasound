"""Run the original transport-lag archive verifier against its frozen docs.

README and CI gained *only* the later physical-study additions. Recover their
original exact bytes for the historical hash check, restore even on failure.
The old archive, manifest, numerical source and verifier are never rewritten.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
CI = ROOT / ".github/workflows/ci.yml"
MANIFEST = ROOT / "runs/transport_lag_information_v1/manifest.json"

README_ADD = (
    "The [paired physical transport-lag protocol](docs/TRANSPORT_LAG_PHYSICAL_PROTOCOL_V1.md)\n"
    "freezes a bubble-free pump-off/pump-on test of those two timings, with native IQ,\n"
    "independent tissue-motion records and predeclared near-wall regions. Its\n"
    "capture inventory is unfilled; there is no physical outcome yet.\n"
)
CI_OLD = "      - run: python scripts/verify_transport_lag_information.py\n"
CI_RPCA_SOURCE_OLD = "--exclude-current-source esdiva_aperture.py\n"
CI_RPCA_SOURCE_NEW = "--exclude-current-source esdiva_aperture.py --exclude-current-source rpca_comparators.py\n"
CI_RPCA_V07_OLD = "                     Path('src/btbench/esdiva_aperture.py')]\n"
CI_RPCA_V07_NEW = ("                     Path('src/btbench/esdiva_aperture.py'),\n"
                   "                     Path('src/btbench/rpca_comparators.py')]\n")
CI_TEMP_OLD = "with TemporaryDirectory() as temp:"
CI_TEMP_NEW = "with TemporaryDirectory(dir=Path('src/btbench')) as temp:"
CI_NEW = (
    "      - run: python scripts/verify_transport_lag_historical_checkout.py\n"
    "      - run: python scripts/verify_transport_lag_physical.py\n"
)


def original_bytes():
    readme, ci = README.read_bytes(), CI.read_bytes()
    addition = README_ADD.encode()
    current_ci = CI_NEW.encode()
    rpca_source = CI_RPCA_SOURCE_NEW.encode()
    rpca_v07 = CI_RPCA_V07_NEW.encode()
    temp_change = CI_TEMP_NEW.encode()
    if (readme.count(addition) != 1 or ci.count(current_ci) != 1 or
            ci.count(rpca_source) != 1 or ci.count(rpca_v07) != 1 or
            ci.count(temp_change) != 1):
        raise ValueError("Unexpected README/CI edits; historical checkout cannot be reconstructed")
    original_readme = readme.replace(addition, b"", 1)
    original_ci = ci.replace(current_ci, CI_OLD.encode(), 1)
    original_ci = original_ci.replace(rpca_source, CI_RPCA_SOURCE_OLD.encode(), 1)
    original_ci = original_ci.replace(rpca_v07, CI_RPCA_V07_OLD.encode(), 1)
    original_ci = original_ci.replace(temp_change, CI_TEMP_OLD.encode(), 1)
    expected = json.loads(MANIFEST.read_text())["sha256"]
    for path, data in (("README.md", original_readme),
                       (".github/workflows/ci.yml", original_ci)):
        if hashlib.sha256(data).hexdigest() != expected[path]:
            raise ValueError(f"Unexpected changes to frozen historical {path}")
    return readme, ci, original_readme, original_ci


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--check-only", action="store_true")
    args = p.parse_args()
    current_readme, current_ci, old_readme, old_ci = original_bytes()
    if args.check_only:
        print("Historical README/CI reconstruction matches frozen SHA-256")
        return
    try:
        README.write_bytes(old_readme)
        CI.write_bytes(old_ci)
        subprocess.run([sys.executable, str(ROOT / "scripts/verify_transport_lag_information.py")],
                       cwd=ROOT, check=True)
    finally:
        README.write_bytes(current_readme)
        CI.write_bytes(current_ci)


if __name__ == "__main__":
    main()
