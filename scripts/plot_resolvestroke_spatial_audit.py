"""Plot the released multi-plane Resolve Stroke tube-mask audit."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


METHODS = ('svd_cut1', 'highpass_100hz')
METHOD_LABELS = ('SVD cut 1', '100 Hz high-pass')
CLIPS = ('−1', '+5', '+10', '+15', '+20')


def plot(report_path, maps_path, prefix):
    report = json.loads(Path(report_path).read_text())
    with np.load(maps_path) as d:
        power, grid, masks = d['power'], d['coordinates'], d['tube_mask']
    if power.shape != (5, 5, 3, 57, 77) or len(report['central_plane_temporal_power']) != 75:
        raise ValueError('Incomplete multi-plane archive')
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    xlabels = [f'{v*1e3:+.1f}' for v in report['grid']['x_m']]

    fig, axes = plt.subplots(2, 4, figsize=(16, 8), constrained_layout=True)
    for col, (method, label) in enumerate((method, label) for method in METHODS for label in ('1', '2')):
        for row, key in enumerate(('auc', 'tpr_at_outside_99th_percentile')):
            values = np.array([[clip['planes'][i]['scores'][method][label]['original'][key]
                                for i in range(5)] for clip in report['spatial_scores']])
            ax = axes[row, col]
            im = ax.imshow(values, cmap='viridis', vmin=0, vmax=1, aspect='auto')
            for y in range(5):
                for x in range(5):
                    ax.text(x, y, f'{values[y, x]:.2f}', ha='center', va='center',
                            fontsize=8, color='black' if values[y, x] > .65 else 'white')
            ax.set_xticks(range(5), xlabels)
            ax.set_yticks(range(5), CLIPS)
            ax.set_title(f'{METHOD_LABELS[col//2]} · {4 if label == "1" else 2} mm')
            if row == 1:
                ax.set_xlabel('Plane x (mm)')
            if col == 0:
                ax.set_ylabel(('Tube/outside AUC' if row == 0 else 'Tube fraction at 1% outside')
                              + '\nClip start (s)')
    fig.colorbar(im, ax=axes[:, :].ravel().tolist(), label='Score', shrink=.7)
    fig.savefig(prefix.with_name(prefix.name + '_scores.png'), dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True, constrained_layout=True)
    colors = ('#575757', '#0a7285', '#c05422')
    for tube_index, ax in enumerate(axes):
        label = str(tube_index + 1)
        for method_index, method in enumerate(report['methods']):
            baseline = np.mean([row['tube_power'][method][label]
                                for row in report['central_plane_temporal_power'][:15]])
            for clip_index, start in enumerate((-1, 5, 10, 15, 20)):
                rows = report['central_plane_temporal_power'][clip_index*15:(clip_index+1)*15]
                t = [start + (r['absolute_frames'][0] - clip_index*4000 + 128)/4000 for r in rows]
                db = [10*np.log10(r['tube_power'][method][label]/baseline) for r in rows]
                ax.plot(t, db, '.-', lw=1.3, ms=3, color=colors[method_index],
                        label=('Unfiltered IQ', 'SVD cut 1', '100 Hz high-pass')[method_index]
                        if clip_index == 0 else None)
        ax.axhline(0, color='gray', lw=.6)
        ax.set_ylabel(f'{4 if label == "1" else 2} mm tube\nPower vs −1 s baseline (dB)')
        ax.grid(alpha=.2)
    axes[0].legend(ncol=3, loc='upper right')
    axes[-1].set_xlabel('Time relative to injection (s); gaps are unrecorded')
    fig.savefig(prefix.with_name(prefix.name + '_temporal.png'), dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(3, 5, figsize=(17, 10), sharex=True, sharey=True,
                             constrained_layout=True)
    for j, name in enumerate(('Unfiltered IQ', 'SVD cut 1', '100 Hz high-pass')):
        row_max = max(float(np.max(power[1, :, j])), 1e-30)
        for i in range(5):
            ax = axes[j, i]
            db = 10*np.log10(np.maximum(power[1, i, j], 1e-30)/row_max)
            y, z = grid[i, 0, :, 1]*1e3, grid[i, :, 0, 2]*1e3
            im = ax.pcolormesh(y, z, db, cmap='magma', shading='auto', vmin=-40, vmax=0)
            ax.contour(y, z, masks[i] == 1, levels=[.5], colors=['#3fe0dc'], linewidths=1.2)
            ax.contour(y, z, masks[i] == 2, levels=[.5], colors=['#8df350'], linewidths=1.2)
            ax.set_aspect('equal')
            if j == 0:
                ax.set_title(f'x = {xlabels[i]} mm')
            if i == 0:
                ax.set_ylabel(name + '\nDepth (mm)')
            if j == 2:
                ax.set_xlabel('Lateral y (mm)')
        fig.colorbar(im, ax=axes[j, :].tolist(), shrink=.55, pad=.01,
                     label='dB relative to row maximum')
    axes[0, 0].invert_yaxis()
    fig.suptitle('Resolve Stroke +5 s · five structural planes · geometric 4/2 mm tube masks')
    fig.savefig(prefix.with_name(prefix.name + '_maps.png'), dpi=120)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--maps', type=Path, required=True)
    parser.add_argument('--prefix', type=Path, required=True)
    args = parser.parse_args()
    plot(args.report, args.maps, args.prefix)
