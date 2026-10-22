import json
from pathlib import Path
import numpy as np
from .core import Acquisition


def save_acquisition(path,acquisition):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    a=acquisition
    arrays={k:getattr(a,k) for k in ("iq","iq_views","timestamps_s","z_m","x_m","tissue_roi","flow_mask",
        "near_wall_mask","blood_true","tissue_true","velocity_true_mps") if getattr(a,k) is not None}
    meta={k:getattr(a,k) for k in ("prf_hz","f0_hz","sound_speed_mps","phase_sign","wall","metadata")}
    np.savez_compressed(path,**arrays,metadata_json=np.array(json.dumps(meta)))


def load_acquisition(path):
    with np.load(path,allow_pickle=False) as f:
        meta=json.loads(str(f["metadata_json"]))
        return Acquisition(**{k:f[k] for k in f.files if k!="metadata_json"},**meta)


def _open_h5(path):
    import h5py
    if str(path).startswith(("https://","http://")):
        import fsspec
        # HTTPFile uses range reads; it raises if the server cannot seek.
        stream=fsspec.open(path,"rb",block_size=256*1024,cache_type="blockcache",
                          client_kwargs={"trust_env":True}).open()
        return h5py.File(stream,"r")
    return h5py.File(path,"r")


def inspect_openhrf(path):
    """Read dataset metadata/inventory only. No channel-array reads."""
    import h5py
    with _open_h5(path) as f:
        result={"source":str(path),"datasets":{},"scalars":{}}
        def visit(name,obj):
            if not isinstance(obj,h5py.Dataset):
                return
            result["datasets"][name]={"shape":list(obj.shape),"dtype":str(obj.dtype),
                "chunks":list(obj.chunks) if obj.chunks else None,
                "decoded_bytes":int(obj.size*obj.dtype.itemsize)}
            if obj.size<=32 and ("scan/" in name or "probe/" in name):
                a=obj[()]
                if isinstance(a,bytes): a=a.decode(errors="replace")
                elif isinstance(a,np.ndarray): a=a.tolist()
                elif isinstance(a,np.generic): a=a.item()
                result["scalars"][name]=a
        f.visititems(visit)
        return result


def export_h5_iq(source,iq_key,output,metadata,time_axis=0,start=0,stop=None):
    """Only explicitly identified unfiltered COMPLEX IQ, never raw_data or dB maps.

    Caller must supply correct compounded/revisit PRF and phase sign. Stored
    velocity/PD maps are not promoted to truth. Multi-view export is separate.
    """
    if any(a in iq_key.lower() for a in ("raw_data","power_doppler","color_doppler","vector_velocity")):
        raise ValueError("This key is not an unfiltered complex beamformed IQ field")
    with _open_h5(source) as f:
        ds=f[iq_key]
        if ds.dtype.kind!="c" or ds.ndim!=3:
            raise ValueError("Need actual complex IQ (3 axes), not RF/envelope/intensity")
        sl=[slice(None)]*3;sl[time_axis]=slice(start,stop)
        iq=np.moveaxis(ds[tuple(sl)],time_axis,0)
    a=Acquisition(iq=iq,metadata={"source":str(source),"iq_key":iq_key,
                  "independent_truth_available":False},**metadata)
    save_acquisition(output,a)
    return a
