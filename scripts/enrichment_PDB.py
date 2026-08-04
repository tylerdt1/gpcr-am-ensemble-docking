#!/usr/bin/env python3
"""
Enrichment Factor (EF) and Early Enrichment Factor (EF') calculator
for GPCR AM retrospective docking benchmarking.

Usage:
    python enrichment_calculator.py <scores.csv> \
        --name_column ligand_id \
        --score_column glide_gscore \
        --output_dir results/ \
        --filename_prefix M2R_PDB_Glide
"""

import pandas as pd
import numpy as np
import math
import argparse
import os
from datetime import datetime


def load_and_sort(csv_file, name_col, score_col, ascending=True):
    df = pd.read_csv(csv_file)
    df = df.sort_values(by=score_col, ascending=ascending).reset_index(drop=True)
    print(f"Loaded {len(df)} compounds from {csv_file}")
    return df


def get_compound_counts(df, name_col, hit_id):
    n_total = len(df)
    n_hits = df[name_col].str.startswith(hit_id).sum()
    print(f"Total compounds: {n_total} | Actives: {n_hits} | Decoys: {n_total - n_hits}")
    return n_total, n_hits


def compute_ef(df, name_col, hit_id, n_total, n_hits, percentages):
    results = {}
    for pct in percentages:
        n_sampled = math.ceil(n_total * pct / 100)
        hits_sampled = df.head(n_sampled)[name_col].str.startswith(hit_id).sum()
        ef = (hits_sampled / n_sampled) / (n_hits / n_total) if n_hits > 0 else 0.0
        results[pct] = {'n_sampled': n_sampled, 'hits_sampled': hits_sampled, 'EF': ef}
    return results


def compute_ef_prime(df, name_col, hit_id, n_total, n_hits, percentages):
    df = df.copy()
    df['percentile_rank'] = (df.index + 1) / n_total * 100
    results = {}
    for pct in percentages:
        n_sampled = math.ceil(n_total * pct / 100)
        top = df.head(n_sampled)
        hits = top[top[name_col].str.startswith(hit_id)]
        hits_sampled = len(hits)
        if hits_sampled > 0:
            apr = hits['percentile_rank'].mean()
            ef_prime = (50.0 / apr) * (hits_sampled / n_hits)
        else:
            apr = 0.0
            ef_prime = 0.0
        results[pct] = {'n_sampled': n_sampled, 'hits_sampled': hits_sampled,
                        'APR': apr, 'EF_prime': ef_prime}
    return results


def save_results(ef_results, ef_prime_results, output_dir, prefix):
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    fname = os.path.join(output_dir, f"{prefix}_enrichment_{timestamp}.txt")

    with open(fname, 'w') as f:
        f.write("EF = (Hits_sampled / N_sampled) / (Hits_total / N_total)\n")
        f.write("EF' = (50 / APR_sampled) * (Hits_sampled / Hits_total)\n\n")

        f.write(f"{'%':<8} {'N_sampled':<12} {'Hits':<8} {'EF':<10} {'EF_prime':<10} {'APR':<8}\n")
        f.write("-" * 56 + "\n")

        for pct in ef_results:
            ef = ef_results[pct]
            efp = ef_prime_results[pct]
            f.write(f"{pct:<8.1f} {ef['n_sampled']:<12} {ef['hits_sampled']:<8} "
                    f"{ef['EF']:<10.3f} {efp['EF_prime']:<10.3f} {efp['APR']:<8.2f}\n")

    print(f"Results saved to: {fname}")


def main():
    parser = argparse.ArgumentParser(
        description='Calculate EF and EF\' for GPCR AM benchmarking.')
    parser.add_argument('csv_file', help='Path to docking scores CSV')
    parser.add_argument('--name_column', default='ligand_id',
                        help='Column containing compound IDs (default: ligand_id)')
    parser.add_argument('--score_column', default='glide_gscore',
                        help='Column containing docking scores (default: glide_gscore)')
    parser.add_argument('--hit_identifier', default='ASD',
                        help='Prefix identifying active compounds (default: ASD)')
    parser.add_argument('--ascending', default=True, type=lambda x: x.lower() != 'false',
                        help='Sort order: True for lower=better (default), False for higher=better')
    parser.add_argument('--output_dir', default='.',
                        help='Output directory (default: current directory)')
    parser.add_argument('--filename_prefix', default='enrichment',
                        help='Prefix for output filename (default: enrichment)')
    args = parser.parse_args()

    percentages = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]

    df = load_and_sort(args.csv_file, args.name_column, args.score_column, args.ascending)
    n_total, n_hits = get_compound_counts(df, args.name_column, args.hit_identifier)

    if n_hits == 0:
        print("No actives found. Check --hit_identifier.")
        return

    ef_results = compute_ef(df, args.name_column, args.hit_identifier,
                             n_total, n_hits, percentages)
    ef_prime_results = compute_ef_prime(df, args.name_column, args.hit_identifier,
                                        n_total, n_hits, percentages)

    print(f"\n{'%':<8} {'EF':<10} {'EF_prime':<10}")
    print("-" * 28)
    for pct in percentages:
        print(f"{pct:<8.1f} {ef_results[pct]['EF']:<10.3f} {ef_prime_results[pct]['EF_prime']:<10.3f}")

    save_results(ef_results, ef_prime_results, args.output_dir, args.filename_prefix)


if __name__ == "__main__":
    main()