"""Render the fixed Resolve Stroke tube-mask comparison from archived power maps."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


METHODS = ('raw', 'svd_cut1', 'highpass_100hz')
LABELS = ('Unfiltered IQ', 'SVD cut 1', '100 Hz high-pass')
CLIP_LABELS = ('−1 s', '+5 s', '+10 s', '+15 s', '+20 s')


def plot(report_path, maps_path, output):
    report = json.loads(Path(report_path).read_text())
    with np.load(maps_path) as data:
        power, grid, mask = data['power'], data['coordinates'], data['tube_mask']
    if (power.shape != (5, 3, *mask.shape) or
            report['methods'] != list(METHODS) or
            len(report['clip_scores']) != 5):
        raise ValueError('Incomplete or changed Resolve Stroke result')
    y, z = grid[0, :, 1] * 1e3, grid[:, 0, 2] * 1e3
    fig, axes = plt.subplots(3, 5, figsize=(17, 10), sharex=True, sharey=True,
                             constrained_layout=True)
    for row, label in enumerate(LABELS):
        row_max = max(float(np.max(power[:, row])), 1e-30)
        for col, clip in enumerate(report['clip_scores']):
            db = 10 * np.log10(np.maximum(power[col, row], 1e-30) / row_max)
            ax = axes[row, col]
            image = ax.pcolormesh(y, z, db, shading='auto', cmap='magma', vmin=-40, vmax=0)
            ax.contour(y, z, mask == 1, levels=[.5], colors=['#3fe0dc'], linewidths=1.2)
            ax.contour(y, z, mask == 2, levels=[.5], colors=['#8df350'], linewidths=1.2)
            ax.set_aspect('equal')
            if row == 0:
                ax.set_title(CLIP_LABELS[col])
            if col == 0:
                ax.set_ylabel(label + '\nDepth (mm)')
            if row == 2:
                ax.set_xlabel('Lateral y (mm)')
        fig.colorbar(image, ax=axes[row, :].tolist(), shrink=.55, pad=.01,
                     label='dB relative to row maximum')
    axes[0, 0].invert_yaxis()
    fig.suptitle('Resolve Stroke CEUS phantom · x = −4.8 mm cross-section · geometric tube masks', fontsize=15)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, metadata={'Creator': 'btbench Resolve Stroke fixed single-view diagnostic'})
    plt.close(fig)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report', type=Path, required=True)
    p.add_argument('--maps', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    plot(a.report, a.maps, a.output)
