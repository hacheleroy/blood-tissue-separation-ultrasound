"""Verify archived milestone bytes, frozen configurations and source identity."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def archive_bytes(root,item):
    if not item.get("parts"):
        return (root/item["path"]).read_bytes()
    parts=[]
    for part in item["parts"]:
        data=(root/part["path"]).read_bytes()
        if len(data)!=part["bytes"] or hashlib.sha256(data).hexdigest()!=part["sha256"]:
            raise ValueError("Archive part mismatch: "+part["path"])
        parts.append(data)
    return b"".join(parts)


def verify(root,require_current_code=False,restore=False):
    root=Path(root)
    manifest=json.loads((root/"runs/fidelity_v03_manifest.json").read_text())
    for item in manifest["archives"]:
        data=archive_bytes(root,item)
        if hashlib.sha256(data).hexdigest()!=item["sha256"]:raise ValueError("Archive hash mismatch: "+item["path"])
        raw=gzip.decompress(data)
        if hashlib.sha256(raw).hexdigest()!=item["decoded_sha256"]:raise ValueError("Decoded hash mismatch")
        report=json.loads(raw)
        records=report.get("records")
        rows=report["results"] if records is None else [r for c in records for r in c["results"]]
        if len(rows)!=item["outputs"]:raise ValueError("Output count mismatch")
        hashes={report["code_sha256"]} if "code_sha256" in report else {c["code_sha256"] for c in records}
        if hashes!={manifest["code_sha256"]}:raise ValueError("Mixed source hashes")
        if item.get("config_path"):
            cfg=(root/item["config_path"]).read_bytes()
            if hashlib.sha256(cfg).hexdigest()!=item["config_sha256"]:raise ValueError("Frozen config hash mismatch")
        for row in rows:
            m=row.get("metrics",{})
            closure=m.get("frozen_component_closure_nrmse",m.get("frozen_mixture_linearity_error"))
            if closure is not None and closure>1e-9:raise ValueError("Component closure failed")
        if restore and item.get("parts"):
            path=root/item["path"]
            if path.exists() and path.read_bytes()!=data:
                raise ValueError("Existing archive differs: "+item["path"])
            if not path.exists():path.write_bytes(data)
    if require_current_code:
        h=hashlib.sha256()
        for p in sorted((root/"src/btbench").glob("*.py")):h.update(p.name.encode());h.update(p.read_bytes())
        if h.hexdigest()!=manifest["code_sha256"]:raise ValueError("Current numerical source differs from archived milestone")
    return {"archives":len(manifest["archives"]),"outputs":sum(x["outputs"] for x in manifest["archives"]),
            "code_sha256":manifest["code_sha256"]}


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("root",nargs="?",default=".")
    p.add_argument("--require-current-code",action="store_true")
    p.add_argument("--restore",action="store_true",help="Restore split archives after checking their hashes")
    args=p.parse_args()
    print(json.dumps(verify(args.root,args.require_current_code,args.restore),indent=2))
