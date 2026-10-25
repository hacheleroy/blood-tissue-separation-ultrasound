"""Render the frozen phase-8 paired gates from its complete decision."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt

plt.rcParams['svg.fonttype'] = 'none'


def plot(decision_path, output_path):
    decision = json.loads(Path(decision_path).read_text())
    if decision['missing'] and not decision['failures']:
        raise ValueError('No eligible failure explains the incomplete phase-8 comparison')
    seeds = (2309, 2311, 2503, 2521, 2707, 2729)
    scenes = ('rigid_fast', 'smooth_slow')
    colors = dict(zip(seeds, ('#4c78a8', '#f58518', '#54a24b', '#b279a2',
                             '#e45756', '#72b7b2')))
    methods = {'svd_cut1': 'o', 'axis1_cut1': 's'}
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.8), sharey=True)
    fig.subplots_adjust(left=.09, right=.98, top=.86, bottom=.29, wspace=.14)
    for ax, scene in zip(axes, scenes):
        for row in decision['rows']:
            if row['scene'] != scene:
                continue
            ax.scatter(row['tissue_leakage_ratio'], row['blood_nrmse_increase'],
                       c=colors[row['seed']], marker=methods[row['method']],
                       edgecolors='black', linewidths=.45, s=58, zorder=3)
        ax.axvline(.9, color='#555555', linestyle='--', linewidth=1)
        ax.axhline(.02, color='#555555', linestyle='--', linewidth=1)
        ax.set_title('Rigid fast flow' if scene == 'rigid_fast' else 'Smooth tissue, slow flow')
        ax.set_xlabel('Tissue leakage ratio to no correction (≤ 0.900)')
        ax.grid(alpha=.2, zorder=0)
        ax.margins(.12)
    axes[0].set_ylabel('Paired blood NRMSE increase (≤ 0.020)')
    from matplotlib.lines import Line2D
    scored_seeds = {row['seed'] for row in decision['rows']}
    legend = [Line2D([0], [0], marker='o', color='none', markerfacecolor=colors[s],
                     markeredgecolor='black', markersize=7, label=str(s))
              for s in seeds if s in scored_seeds]
    legend += [Line2D([0], [0], marker=m, color='none', markerfacecolor='white',
                      markeredgecolor='black', markersize=7, label=n)
               for n, m in [('SVD cutoff 1', 'o'), ('Axis 1 cutoff 1', 's')]]
    fig.legend(handles=legend, ncol=4, loc='lower center',
               bbox_to_anchor=(.5, .095), frameon=False)
    fig.suptitle('Phase 8 · guarded calibration and measured phase selection', fontsize=13)
    omitted = sorted(set(seeds) - scored_seeds)
    if omitted:
        fig.text(.5, .02, 'Unscored: ' +
                 ', '.join(decision['missing']) + '; ' +
                 ', '.join(decision['failures']) +
                 '. No prospective release.', ha='center', fontsize=8)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches='tight', metadata={'Creator':'btbench phase-8 decision'})
    plt.close(fig)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--decision', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    plot(a.decision, a.output)
