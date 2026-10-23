"""Verify and beamform the declared bounded Waterloo windows in small RF batches."""
import argparse
import hashlib
import json
from pathlib import Path
from btbench.io import verify_waterloo_geometry,export_openhrf_rf_window


def main():
    p=argparse.ArgumentParser();p.add_argument('config');p.add_argument('template_geometry');p.add_argument('output')
    args=p.parse_args();config=json.loads(Path(args.config).read_text());template=json.loads(Path(args.template_geometry).read_text())
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    for spec in config['acquisitions']:
        target=out/(spec['name']+'.npz')
        if target.exists():raise FileExistsError('Refusing to overwrite an existing reconstructed acquisition: '+str(target))
        url=('https://huggingface.co/datasets/nvidia/Open-H-RF/resolve/'+config['dataset_revision']+
             '/waterloo-carotid/data/'+spec['file'])
        geometry=verify_waterloo_geometry(url,template,config['frames_start_stop'],
            config['reconstruction_grid_stride'],config['frame_batch_size'])
        if geometry['source_metadata']['subject_id']!=spec['subject_id']:
            raise ValueError('Unexpected subject in declared acquisition')
        encoded=json.dumps(geometry,sort_keys=True,indent=2)+'\n'
        (out/(spec['name']+'_geometry.json')).write_text(encoded)
        print(json.dumps({'reconstructing':spec['name'],'subject':spec['subject_id'],
                          'geometry_sha256':hashlib.sha256(encoded.encode()).hexdigest()}),flush=True)
        a=export_openhrf_rf_window(url,geometry,target)
        print(json.dumps({'reconstructed':spec['name'],'iq_shape':list(a.iq.shape),
                          'decoded_raw_bytes':a.metadata['decoded_raw_bytes'],
                          'maximum_decoded_batch_bytes':a.metadata['maximum_decoded_batch_bytes']}),flush=True)


if __name__=='__main__':main()
