import json
from pathlib import Path
import numpy as np
from .core import Acquisition


def save_acquisition(path,acquisition):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    a=acquisition
    arrays={k:getattr(a,k) for k in ("iq","iq_views","timestamps_s","z_m","x_m","tissue_roi","flow_mask",
        "near_wall_mask","blood_true","tissue_true","velocity_true_mps","tx_angles_rad","rx_angles_rad",
        "view_offsets_s","blood_views_true","tissue_views_true") if getattr(a,k) is not None}
    meta={k:getattr(a,k) for k in ("prf_hz","f0_hz","sound_speed_mps","phase_sign","wall","metadata")}
    np.savez_compressed(path,**arrays,metadata_json=np.array(json.dumps(meta)))


def load_acquisition(path):
    with np.load(path,allow_pickle=False) as f:
        meta=json.loads(str(f["metadata_json"]))
        return Acquisition(**{k:f[k] for k in f.files if k!="metadata_json"},**meta)


def _open_h5(path,block_size=256*1024):
    import h5py
    try:
        import hdf5plugin  # registers provider Blosc filters, including metadata arrays
    except ImportError:
        pass  # uncompressed files and metadata-only inspection still work
    if str(path).startswith(("https://","http://")):
        import fsspec
        # HTTPFile uses range reads; it raises if the server cannot seek.
        stream=fsspec.open(path,"rb",block_size=block_size,cache_type="blockcache",
                          cache_options={"maxblocks":16},
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


def selected_view_timing(intervals_s,transmit_indices,start=0,stop=None):
    """Recover absolute event times from recorded sequential TX intervals.

    Every selected view must have a constant offset and uniform revisit rate.
    Bursts/gaps and varying schedules must be split before this adapter is used.
    """
    dt=np.asarray(intervals_s,float)
    tx=np.asarray(transmit_indices,int)
    stop=len(dt) if stop is None else stop
    if (dt.ndim!=2 or not 0<=start<stop<=len(dt) or stop-start<3 or
            tx.ndim!=1 or not len(tx) or len(set(tx))!=len(tx) or np.any(tx<0) or np.any(tx>=dt.shape[1])):
        raise ValueError("Invalid frame range or transmit selection")
    flat=dt[:stop].ravel()
    if not np.isfinite(flat).all() or np.any(flat[:-1]<=0):
        raise ValueError("Transmit intervals must be finite and positive")
    events=np.r_[0.,np.cumsum(flat[:-1])].reshape(stop,dt.shape[1])[start:stop,tx]
    times=events[:,0]
    step=np.mean(np.diff(times))
    offsets=events-times[:,None]
    if (not np.allclose(np.diff(events,axis=0),step,rtol=1e-4,atol=1e-9) or
            not np.allclose(offsets,offsets[:1],rtol=1e-4,atol=1e-9)):
        raise ValueError("Gapped/nonuniform view timing: split contiguous acquisition segments")
    return times,offsets[0],float(1/step)


def export_openhrf_rf_window(source,geometry,output):
    """Bounded real-RF DAS export with all views embedded in the acquisition.

    Geometry is an explicit verified far-field plane-wave contract. Hardware
    gain and complex phases are preserved; no wall filtering, normalization,
    truth-based calibration or de-aliasing is applied.
    """
    from .beamform import plane_wave_das
    g=dict(geometry)
    if not g.get("time_reference_verified") or not g.get("full_receive_aperture_verified"):
        raise ValueError("Verify transmit time reference and simultaneous receive aperture")
    start,stop=g["frames_start_stop"]
    if int(start)!=start or int(stop)!=stop or not 0<=start<stop or not 3<=stop-start<=2048:
        raise ValueError("RF adapter requires 3..2048 contiguous frames")
    tx=np.asarray(g["transmit_indices"],int)
    if tx.ndim!=1 or not len(tx) or len(set(tx))!=len(tx):raise ValueError("Select distinct transmit indices")
    root=g.get("root","").strip("/")
    def key(s):return root+"/"+s if root else s
    with _open_h5(source,block_size=g.get("http_block_bytes",4*1024**2)) as f:
        ds=f[key("data/raw_data")]
        if ds.ndim!=5 or ds.shape[-1]!=1 or ds.dtype.kind not in "fiu":
            raise ValueError("Need real RF (time,transmit,sample,element,1); DDC IQ needs a separate phase model")
        if stop>len(ds) or np.any(tx<0) or np.any(tx>=ds.shape[1]):raise ValueError("RF selection exceeds stored dimensions")
        # Waterloo chunks contain all transmits; read/decompress each frame once.
        raw_bytes=(stop-start)*np.prod(ds.shape[1:])*ds.dtype.itemsize
        batch_size=g.get("frame_batch_size",stop-start)
        if int(batch_size)!=batch_size or batch_size<1:
            raise ValueError("frame_batch_size must be a positive integer")
        batch_size=min(int(batch_size),stop-start)
        batch_bytes=batch_size*np.prod(ds.shape[1:])*ds.dtype.itemsize
        if batch_bytes>g.get("max_raw_bytes",512*1024**2):
            raise MemoryError("Decoded RF batch exceeds memory budget")
        dt=f[key("scan/time_to_next_transmit")][:stop]
        times,offsets,prf=selected_view_timing(dt,tx,start,stop)
        if "revisit_prf_hz" in g and not np.isclose(prf,g["revisit_prf_hz"],rtol=1e-4):
            raise ValueError("Declared PRF disagrees with recorded event timing")
        if "frame_times_s" in g and not np.allclose(times,g["frame_times_s"],rtol=1e-6,atol=1e-9):
            raise ValueError("Declared frame times disagree with recorded event timing")
        for field,h5name in (("fs_hz","sampling_frequency"),("f0_hz","center_frequency"),("sound_speed_mps","sound_speed")):
            value=float(f[key("scan/"+h5name)][()])
            if field in g and not np.isclose(g[field],value,rtol=1e-5):
                raise ValueError(f"Declared {field} disagrees with stored scan metadata")
            g[field]=value
        angles=np.asarray(f[key("scan/polar_angles")][()])[tx]
        if not np.allclose(angles,g["tx_angles_rad"],atol=1e-6):
            raise ValueError("Declared transmit angles disagree with stored scan metadata")
        initial=np.asarray(f[key("scan/initial_times")][()])[tx]
        if not np.allclose(initial,g["initial_times_s"],rtol=1e-5,atol=.1/g["fs_hz"]):
            raise ValueError("Declared A/D start time disagrees with stored metadata")
        elements=np.asarray(f["probe/probe_geometry"][()])[:,[2,0]]
        if elements.shape!=(ds.shape[3],2):raise ValueError("Receive geometry does not match channels")
        attrs={k:str(v) for k,v in f.attrs.items()}
        z,x=np.asarray(g["z_m"]),np.asarray(g["x_m"])
        views=np.empty((stop-start,len(tx),len(z),len(x)),complex)
        for first in range(start,stop,batch_size):
            last=min(first+batch_size,stop)
            raw=np.asarray(ds[first:last])[...,0][:,tx]
            views[first-start:last-start]=plane_wave_das(raw,elements,z,x,angles,g["fs_hz"],g["f0_hz"],initial,
                g["tx_reference_times_s"],g.get("tx_origins_zx_m"),g["sound_speed_mps"],g.get("f_number",1.5))
            if g.get("report_progress",False):
                print(json.dumps({"beamformed_frames":last-start,"total_frames":stop-start}),flush=True)
    acquisition=Acquisition(iq=views.mean(axis=1),iq_views=views,prf_hz=prf,f0_hz=g["f0_hz"],
        sound_speed_mps=g["sound_speed_mps"],phase_sign=g["phase_sign"],timestamps_s=times,
        tx_angles_rad=angles,rx_angles_rad=np.zeros(len(tx)),view_offsets_s=offsets,z_m=z,x_m=x,
        metadata={"source":str(source),"source_attributes":attrs,"unfiltered_iq":True,
            "subject_id":g.get("source_metadata",{}).get("subject_id"),
            "verified_source_metadata":g.get("source_metadata"),
            "frames_start_stop":[start,stop],"transmit_indices":tx.tolist(),"decoded_raw_bytes":int(raw_bytes),
            "frame_batch_size":batch_size,"maximum_decoded_batch_bytes":int(batch_bytes),
            "contrast_status":g.get("contrast_status","unconfirmed"),"independent_truth_available":False,
            "time_reference_evidence":g.get("time_reference_evidence"),
            "receive_aperture_evidence":g.get("receive_aperture_evidence"),
            "amplitude_processing":"hardware gain retained; no frame normalization",
            "baseline":"coherent mean of embedded uncompounded views",
            "beamformer":"CPU far-field plane-wave DAS, normal receive; experimental phase validation pending"})
    save_acquisition(output,acquisition)
    return acquisition


def verify_waterloo_geometry(source,template,frames_start_stop=(0,1024),grid_stride=2,frame_batch_size=16):
    """Verify each selected nonmultiplexed Waterloo file before RF export.

    Fits the stored per-element transmit-delay reference; does not calibrate
    experimental Doppler phase. The template supplies the declared grid and
    IQ sign convention already checked with a moving-point simulation.
    """
    g=dict(template);root=g.get('root','tracks/track_0').strip('/')
    with _open_h5(source) as f:
        description=str(f.attrs.get('description','')).lower()
        if 'us4r' not in description or 'lite' in description:
            raise ValueError('This adapter requires explicitly identified non-lite us4r acquisition')
        scan=f[root+'/scan'];raw=f[root+'/data/raw_data']
        elements=np.asarray(f['probe/probe_geometry'][()])[:,[2,0]]
        tx=np.asarray(g['transmit_indices'],int)
        angles=np.asarray(scan['polar_angles'][()])[tx]
        delays=np.asarray(scan['t0_delays'][()])[tx]
        apod=np.asarray(scan['tx_apodizations'][()])[tx]
        origins=np.asarray(scan['transmit_origins'][()])[tx][:,[2,0]]
        fs=float(scan['sampling_frequency'][()]);c=float(scan['sound_speed'][()])
        if raw.shape[3]!=len(elements) or delays.shape!=(len(tx),len(elements)) or not np.all(apod!=0):
            raise ValueError('Unverified element mapping, transmit delay shape or aperture')
        direction=np.column_stack((np.cos(angles),np.sin(angles)))
        relative=np.einsum('aep,ap->ae',elements[None,:,:]-origins[:,None,:],direction)/c
        reference=np.mean(delays-relative,axis=1)
        residual=float(np.max(abs(delays-relative-reference[:,None]))*fs)
        if residual>.1:raise ValueError('Stored delays do not fit the declared plane-wave time reference')
        def text(key):
            v=f['metadata/'+key][()]
            return v.decode() if isinstance(v,bytes) else str(v)
        metadata={'subject_id':text('subject/id'),'anatomy':text('annotations/anatomy'),
                  'view':text('annotations/view'),'condition':text('annotations/label'),
                  'raw_shape':list(raw.shape),'source_description':description,
                  'delay_fit_max_residual_adc_samples':residual}
        g.update(fs_hz=fs,f0_hz=float(scan['center_frequency'][()]),sound_speed_mps=c,
                 tx_angles_rad=angles.tolist(),initial_times_s=np.asarray(scan['initial_times'][()])[tx].tolist(),
                 tx_reference_times_s=reference.tolist(),tx_origins_zx_m=origins.tolist())
    if int(grid_stride)!=grid_stride or grid_stride<1:raise ValueError('grid_stride must be a positive integer')
    g.update(frames_start_stop=list(frames_start_stop),frame_batch_size=frame_batch_size,
             z_m=template['z_m'][::grid_stride],x_m=template['x_m'][::grid_stride],
             max_raw_bytes=64*1024**2,time_reference_verified=True,full_receive_aperture_verified=True,
             source_metadata=metadata,report_progress=True,
             time_reference_evidence=f'Per-file stored delay fit: maximum residual {residual:.6g} ADC samples; experimental phase calibration pending',
             receive_aperture_evidence='Per-file non-lite us4r description, 128 channels matched to probe, full transmit apodization')
    return g
