#!/usr/bin/env python3
"""
Ensemble docking enrichment calculator with PMF reweighting (BEmin/BEavg).
Computes EF and EF' from GaMD ensemble docking scores.

Usage:
    python ensemble_enrichment_calculator.py <scores_dir> \
        --pmf_file M2R_ensemble_PMF.xvg \
        --frame_pattern "*cluster{}-HTVS-last1_pv-best.csv" \
        --score_column r_i_glide_gscore \
        --name_column NAME \
        --output_dir results/ \
        --filename_prefix M2R_Glide
"""

import pandas as pd
import numpy as np
import math
import argparse
import os
import glob
from datetime import datetime


def load_pmf(pmf_file):
    pmf = {}
    with open(pmf_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('@'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    frame_int = int(float(parts[0]))
                    if frame_int < 0:
                        continue
                    frame_id = f"{frame_int:02d}"
                    pmf[frame_id] = float(parts[1])
                except ValueError:
                    continue
    print(f"Loaded PMF values for {len(pmf)} frames from {pmf_file}")
    return pmf


def compute_be_scores(input_dir, pmf, frame_pattern, score_col, name_col):
    ligand_scores = {}
    framewise = {}

    for frame_id, pmf_val in pmf.items():
        pattern = frame_pattern.format(frame_id)
        matches = glob.glob(os.path.join(input_dir, pattern))
        if not matches:
            print(f"  No file found for frame {frame_id} (pattern: {pattern})")
            continue

        df = pd.read_csv(matches[0], low_memory=False)
        if score_col not in df.columns or name_col not in df.columns:
            print(f"  Missing columns in {matches[0]}, skipping")
            continue

        for _, row in df.iterrows():
            name = row[name_col]
            score = row[score_col]
            if pd.isna(score):
                continue
            be = score + pmf_val
            ligand_scores.setdefault(name, []).append(be)
            framewise.setdefault(name, {})[frame_id] = be

        print(f"  Frame {frame_id}: {len(df)} compounds (PMF = {pmf_val:.3f} kcal/mol)")

    if not ligand_scores:
        raise ValueError("No valid frames processed. Check --frame_pattern and --input_dir.")

    records = [{'NAME': name,
                'BE_min': min(scores),
                'BE_avg': sum(scores) / len(scores),
                'n_frames': len(scores)}
               for name, scores in ligand_scores.items()]

    combined = pd.DataFrame(records)
    framewise_df = pd.DataFrame.from_dict(framewise, orient='index')
    framewise_df.index.name = 'NAME'
    framewise_df = framewise_df.reset_index()

    return combined, framewise_df


def compute_ef(df, name_col, score_col, hit_id, percentages, ascending=True):
    df = df.sort_values(score_col, ascending=ascending).reset_index(drop=True)
    n_total = len(df)
    n_hits = df[name_col].str.startswith(hit_id).sum()

    if n_hits == 0:
        print(f"No actives found with prefix '{hit_id}'")
        return None, n_total, 0

    df['percentile_rank'] = (df.index + 1) / n_total * 100
    results = {}

    for pct in percentages:
        n_sampled = math.ceil(n_total * pct / 100)
        top = df.head(n_sampled)
        hits = top[top[name_col].str.startswith(hit_id)]
        hits_sampled = len(hits)

        ef = (hits_sampled / n_sampled) / (n_hits / n_total) if hits_sampled > 0 else 0.0

        if hits_sampled > 0:
            apr = hits['percentile_rank'].mean()
            ef_prime = (50.0 / apr) * (hits_sampled / n_hits)
        else:
            apr = 0.0
            ef_prime = 0.0

        results[pct] = {
            'n_sampled': n_sampled,
            'hits_sampled': hits_sampled,
            'EF': ef,
            'EF_prime': ef_prime,
            'APR': apr
        }

    return results, n_total, n_hits


def save_results(ef_results, score_type, output_dir, prefix, n_total, n_hits):
    fname = os.path.join(output_dir, f'{prefix}_enrichment_{score_type}.txt')
    with open(fname, 'w') as f:
        f.write(f"Re-ranking method: {score_type}\n")
        f.write(f"Total compounds: {n_total} | Actives: {n_hits} | Decoys: {n_total - n_hits}\n\n")
        f.write("EF  = (Hits_sampled / N_sampled) / (Hits_total / N_total)\n")
        f.write("EF' = (50 / APR_sampled) * (Hits_sampled / Hits_total)\n\n")
        f.write(f"{'%':<8} {'N_sampled':<12} {'Hits':<8} {'EF':<10} {'EF_prime':<10} {'APR':<8}\n")
        f.write("-" * 56 + "\n")
        for pct, v in ef_results.items():
            f.write(f"{pct:<8.1f} {v['n_sampled']:<12} {v['hits_sampled']:<8} "
                    f"{v['EF']:<10.3f} {v['EF_prime']:<10.3f} {v['APR']:<8.2f}\n")
    print(f"Saved: {fname}")


def main():
    parser = argparse.ArgumentParser(
        description='Ensemble docking EF/EF\' calculator with PMF reweighting.')
    parser.add_argument('input_dir',
                        help='Directory containing per-cluster docking score CSVs')
    parser.add_argument('--pmf_file', required=True,
                        help='PMF file in XVG format (e.g. M2R_ensemble_PMF.xvg)')
    parser.add_argument('--frame_pattern', default='*cluster{}-HTVS-last1_pv-best.csv',
                        help='Glob pattern for cluster CSVs with {} as zero-padded frame ID placeholder')
    parser.add_argument('--score_column', default='r_i_glide_gscore',
                        help='Docking score column name (default: r_i_glide_gscore)')
    parser.add_argument('--name_column', default='NAME',
                        help='Compound ID column name (default: NAME)')
    parser.add_argument('--hit_identifier', default='ASD',
                        help='Prefix identifying active compounds (default: ASD)')
    parser.add_argument('--ascending', default=True,
                        type=lambda x: x.lower() != 'false',
                        help='Sort order: True for lower=better (default), False for higher=better')
    parser.add_argument('--output_dir', default='.',
                        help='Output directory (default: current directory)')
    parser.add_argument('--filename_prefix', default='ensemble',
                        help='Prefix for output filenames (default: ensemble)')
    args = parser.parse_args()

    percentages = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
    os.makedirs(args.output_dir, exist_ok=True)

    pmf = load_pmf(args.pmf_file)
    combined, framewise = compute_be_scores(
        args.input_dir, pmf, args.frame_pattern, args.score_column, args.name_column)

    combined[['NAME', 'BE_min']].sort_values('BE_min').to_csv(
        os.path.join(args.output_dir, 'BE_min_ranked.csv'), index=False)
    combined[['NAME', 'BE_avg']].sort_values('BE_avg').to_csv(
        os.path.join(args.output_dir, 'BE_avg_ranked.csv'), index=False)
    framewise.to_csv(
        os.path.join(args.output_dir, 'BE_framewise_scores.csv'), index=False)
    print(f"\nSaved BE_min_ranked.csv, BE_avg_ranked.csv, BE_framewise_scores.csv")

    for score_type in ['BE_min', 'BE_avg']:
        print(f"\n── {score_type} ──")
        ef_results, n_total, n_hits = compute_ef(
            combined, 'NAME', score_type, args.hit_identifier, percentages, args.ascending)
        if ef_results:
            print(f"{'%':<8} {'EF':<10} {'EF_prime':<10}")
            print("-" * 28)
            for pct, v in ef_results.items():
                print(f"{pct:<8.1f} {v['EF']:<10.3f} {v['EF_prime']:<10.3f}")
            save_results(ef_results, score_type, args.output_dir,
                        args.filename_prefix, n_total, n_hits)


if __name__ == '__main__':
    main()
