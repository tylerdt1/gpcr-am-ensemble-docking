#!/usr/bin/env python3
"""
ROC-AUC and logAUC calculator for virtual screening benchmarking.

Usage:
    python auc_calculator.py \
        --csv_files M2R_PDB_Glide.csv M2R_Ensemble_Glide_BEmin.csv M2R_Ensemble_Glide_BEavg.csv \
        --method_names PDB BEmin BEavg \
        --score_columns glide_gscore BE_min BE_avg \
        --program glide \
        --target_name M2R \
        --save_plot M2R_Glide_ROC.png
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

LOGAUC_MIN = 0.001
LOGAUC_MAX = 1.0
RANDOM_LOGAUC = (LOGAUC_MAX - LOGAUC_MIN) / np.log(10) / np.log10(LOGAUC_MAX / LOGAUC_MIN)

COLORS = {
    'glide':  '#c0392b',
    'vina':   '#27ae60',
    'dock':   '#2980b9',
    'boltz2': '#e67e22',
}


def roc_curve(y_true, y_scores):
    idx = np.argsort(-y_scores)
    y = np.array(y_true)[idx]
    n_act, n_dec = y.sum(), len(y) - y.sum()
    pct_a = np.concatenate([[0], np.cumsum(y)]) / n_act * 100
    pct_d = np.concatenate([[0], np.cumsum(1 - y)]) / n_dec * 100
    return pct_d, pct_a


def interpolate(points):
    i = next((i for i, p in enumerate(points) if p[0] >= 0.1), None)
    if not i:
        return points
    slope = (points[i][1] - points[i-1][1]) / (points[i][0] - points[i-1][0])
    b = points[i][1] - slope * points[i][0]
    points.insert(i, [0.100001, slope * 0.100001 + b])
    return points


def calc_auc(pct_d, pct_a):
    pts = interpolate([[d, a] for d, a in zip(pct_d, pct_a)])
    return sum((p[0]-lp[0])/100 * (lp[1]+(p[1]-lp[1])/2)/100
               for p, lp in zip(pts[1:], pts[:-1]))


def calc_logauc(pct_d, pct_a):
    pts = interpolate([[d/100, a/100] for d, a in zip(pct_d, pct_a)
                       if LOGAUC_MIN*100 <= d <= LOGAUC_MAX*100])
    if len(pts) < 2:
        return 0.0
    area = sum((p[1]-lp[1])/np.log(10) +
               (p[1]-(p[1]-lp[1])/(p[0]-lp[0])*p[0]) * (np.log10(p[0])-np.log10(lp[0]))
               for p, lp in zip(pts[1:], pts[:-1]) if p[0]-lp[0] > 1e-6)
    return area / np.log10(LOGAUC_MAX/LOGAUC_MIN) - RANDOM_LOGAUC


def shade(hex_color, factor):
    r, g, b = [int(hex_color.lstrip('#')[i:i+2], 16) for i in (0, 2, 4)]
    return '#{:02x}{:02x}{:02x}'.format(*[int(min(255, max(0, c*factor))) for c in (r,g,b)])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv_files', nargs='+', required=True)
    parser.add_argument('--method_names', nargs='+', required=True)
    parser.add_argument('--score_columns', nargs='+', default=None,
                        help='Score column per CSV file. If not provided, uses --score_column for all.')
    parser.add_argument('--score_column', default='glide_gscore',
                        help='Default score column if --score_columns not provided.')
    parser.add_argument('--name_column', default='ligand_id')
    parser.add_argument('--program', default='glide',
                        choices=['glide', 'vina', 'dock', 'boltz2'])
    parser.add_argument('--hit_identifier', default='ASD')
    parser.add_argument('--target_name', default='Target')
    parser.add_argument('--output_dir', default='.')
    parser.add_argument('--save_plot', default=None)
    args = parser.parse_args()

    if len(args.csv_files) != len(args.method_names):
        raise ValueError("--csv_files and --method_names must have the same number of entries")

    score_cols = args.score_columns if args.score_columns else [args.score_column] * len(args.csv_files)

    if len(score_cols) != len(args.csv_files):
        raise ValueError("--score_columns must have the same number of entries as --csv_files")

    os.makedirs(args.output_dir, exist_ok=True)
    base = COLORS.get(args.program, '#7f8c8d')
    shades = [1.5, 1.0, 0.6, 0.4]
    metrics = []

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    for idx, (csv_file, method, score_col) in enumerate(zip(args.csv_files, args.method_names, score_cols)):
        df = pd.read_csv(csv_file)
        y_true = df[args.name_column].str.startswith(args.hit_identifier).astype(int).values
        scores = df[score_col].values
        if args.program != 'boltz2':
            scores = -scores

        pct_d, pct_a = roc_curve(y_true, scores)
        auc = calc_auc(pct_d, pct_a)
        logauc = calc_logauc(pct_d, pct_a)
        metrics.append({'Method': method, 'AUC (%)': round(auc*100, 3),
                        'logAUC (%)': round(logauc*100, 3)})

        c = shade(base, shades[idx % len(shades)])
        ax1.plot(pct_d, pct_a, color=c, lw=2.0, ls='--', alpha=0.85,
                 label=f'{method} (AUC={auc*100:.2f}%)')
        mask = pct_d <= LOGAUC_MAX * 100
        ax2.plot(pct_d[mask], pct_a[mask], color=c, lw=2.0, ls='--', alpha=0.85,
                 label=f'{method} (logAUC={logauc*100:.2f}%)')

    for ax in (ax1, ax2):
        ax.set_ylabel('% Actives Found', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10, framealpha=0.95)
        ax.grid(True, alpha=0.3)
        ax.set_ylim([0, 100])
        ax.tick_params(labelsize=10)

    ax1.plot([0, 100], [0, 100], 'k--', lw=1.2, alpha=0.4, label='Random')
    ax1.set_xlabel('% Decoys Found', fontsize=12, fontweight='bold')
    ax1.set_title(f'{args.target_name} — ROC Curve', fontsize=13, fontweight='bold')
    ax1.set_xlim([0, 100])

    x_rand = np.logspace(np.log10(LOGAUC_MIN*100), np.log10(LOGAUC_MAX*100), 100)
    ax2.plot(x_rand, x_rand, 'k--', lw=1.2, alpha=0.4, label='Random')
    ax2.set_xlabel('% Decoys Found (log scale)', fontsize=12, fontweight='bold')
    ax2.set_title(f'{args.target_name} — Semilogarithmic ROC Curve',
                  fontsize=13, fontweight='bold')
    ax2.set_xscale('log')
    ax2.set_xlim([LOGAUC_MIN*100, LOGAUC_MAX*100])

    plt.tight_layout()
    if args.save_plot:
        path = os.path.join(args.output_dir, args.save_plot)
        plt.savefig(path, dpi=300, bbox_inches='tight')
        print(f"Saved: {path}")
    plt.show()

    metrics_df = pd.DataFrame(metrics)
    print(f"\nMetrics — {args.target_name} ({args.program})")
    print(metrics_df.to_string(index=False))
    print(f"Random logAUC baseline: {RANDOM_LOGAUC*100:.2f}%")

    out = os.path.join(args.output_dir, f'{args.target_name}_{args.program}_metrics.csv')
    metrics_df.to_csv(out, index=False)
    print(f"Saved: {out}")


if __name__ == '__main__':
    main()