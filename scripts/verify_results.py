"""Verify archived milestone bytes, frozen configurations and source identity."""
import argparse
import gzip
import hashlib
import json
import lzma
import math
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


def verify(root,require_current_code=False,restore=False,manifest_path="runs/fidelity_v03_manifest.json",
           exclude_current_source=()):
    root=Path(root)
    manifest=json.loads((root/manifest_path).read_text())
    for item in manifest["archives"]:
        data=archive_bytes(root,item)
        if hashlib.sha256(data).hexdigest()!=item["sha256"]:raise ValueError("Archive hash mismatch: "+item["path"])
        codec=item.get("codec","gzip")
        if codec not in ("gzip","xz"):raise ValueError("Unsupported archive codec: "+codec)
        raw=(lzma.decompress if codec=="xz" else gzip.decompress)(data)
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
            for key in item.get("required_metrics",[]):
                if row.get("status")!="ok" or m.get(key) is None or not math.isfinite(m[key]):
                    raise ValueError("Missing/invalid primary metric: "+key)
            if item.get("independent_source_truth_available") is False:
                if m.get("frozen_component_metrics_available") or "blood_iq_nrmse" in m:
                    raise ValueError("Experimental source truth was incorrectly promoted")
            closure=m.get("frozen_component_closure_nrmse",m.get("frozen_mixture_linearity_error"))
            if closure is not None and closure>1e-9:raise ValueError("Component closure failed")
        if restore and item.get("parts"):
            path=root/item["path"]
            if path.exists() and path.read_bytes()!=data:
                raise ValueError("Existing archive differs: "+item["path"])
            if not path.exists():path.write_bytes(data)
    if require_current_code:
        h=hashlib.sha256()
        excluded=set(exclude_current_source)
        current={p.name:p for p in (root/"src/btbench").glob("*.py")}
        if excluded-current.keys():raise ValueError("Excluded source does not exist")
        for name,p in sorted(current.items()):
            if name not in excluded:h.update(name.encode());h.update(p.read_bytes())
        if h.hexdigest()!=manifest["code_sha256"]:raise ValueError("Current numerical source differs from archived milestone")
    for item in manifest.get('artifacts',[]):
        data=(root/item['path']).read_bytes()
        if len(data)!=item['bytes'] or hashlib.sha256(data).hexdigest()!=item['sha256']:
            raise ValueError('Derived artifact mismatch: '+item['path'])
    return {"archives":len(manifest["archives"]),"outputs":sum(x["outputs"] for x in manifest["archives"]),
            "code_sha256":manifest["code_sha256"]}


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("root",nargs="?",default=".")
    p.add_argument("--require-current-code",action="store_true")
    p.add_argument("--exclude-current-source",action="append",default=[],
                   help="New module absent from an old source hash; verify it under its own frozen archive")
    p.add_argument("--restore",action="store_true",help="Restore split archives after checking their hashes")
    p.add_argument("--manifest",default="runs/fidelity_v03_manifest.json")
    args=p.parse_args()
    print(json.dumps(verify(args.root,args.require_current_code,args.restore,args.manifest,
                            args.exclude_current_source),indent=2))
