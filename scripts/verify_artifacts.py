"""Check portable source/wheel contents and write their SHA256 manifest."""
import argparse
import hashlib
import json
from pathlib import Path,PurePosixPath
import tarfile
import zipfile


def verify(directory):
    directory=Path(directory);manifest=[]
    wheels=list(directory.glob("*.whl"));sources=list(directory.glob("*.tar.gz"))
    if len(wheels)!=1 or len(sources)!=1:raise ValueError("Expected exactly one wheel and one source distribution")
    for path in wheels+sources:
        if path.suffix==".whl":
            with zipfile.ZipFile(path) as f:names=f.namelist()
            required=["btbench/measurements.py","btbench/axes.py","btbench/experimental.py","btbench/separation_experiment.py"]
        else:
            with tarfile.open(path) as f:
                members=f.getmembers()
                if any(m.issym() or m.islnk() for m in members):raise ValueError("Archive links are not permitted")
                names=[m.name for m in members]
            root=names[0].split('/')[0]+"/"
            required=[root+p for p in ["src/btbench/core.py","tests/test_frozen.py","tests/test_multiview_io.py",
                "tests/test_separation_study.py","configs/separation_fidelity.json","scripts/verify_artifacts.py"]]
        for name in names:
            p=PurePosixPath(name)
            if p.is_absolute() or ".." in p.parts:raise ValueError("Unsafe archive path: "+name)
            if p.suffix.lower() in (".npz",".npy",".h5",".hdf5",".mat"):
                raise ValueError("Acquisition/generated arrays leaked into package: "+name)
        missing=set(required)-set(names)
        if missing:raise ValueError("Missing package entries: "+str(sorted(missing)))
        manifest.append({"file":path.name,"bytes":path.stat().st_size,
            "sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"entries":len(names)})
    (directory/"ARTIFACT_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
    return manifest


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("directory");args=p.parse_args()
    print(json.dumps(verify(args.directory),indent=2))
