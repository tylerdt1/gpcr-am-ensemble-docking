import pandas as pd
import os

# ── File paths ──────────────────────────────────────────────────────────────
SMILES_FILE = '/work/users/t/y/tylerdt/Schrodinger/LigPrep/CCR2-no_dup-ligprep/CCR2-no_dup-combined.smi'
OUT_DIR     = '/work/users/t/y/tylerdt/clean_ensemble/CCR2/'

GLIDE_DIR   = '/work/users/t/y/tylerdt/Schrodinger/Glide-HTVS-clusters/Enrichment-HTVS-clusters/CCR2-clusters-HTVS-lig_25-res-enrichment/Nclusters_10/'
VINA_DIR    = '/users/t/y/tylerdt/GPCR_Retrospective_Docking/AutoDockVina/CCR2/ensemble_predictions/enrichment/Nclusters_10/'

os.makedirs(OUT_DIR, exist_ok=True)

# ── Load SMILES file ─────────────────────────────────────────────────────────
rows = []
with open(SMILES_FILE, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split('\t') if '\t' in line else line.split(' ')
        if len(parts) >= 2:
            rows.append({'SMILES': parts[0], 'ligand_id': parts[1]})

smiles_df = pd.DataFrame(rows)
smiles_df['is_active'] = smiles_df['ligand_id'].str.startswith('ASD').astype(int)
print(f"SMILES loaded: {len(smiles_df)} compounds ({smiles_df['is_active'].sum()} actives)")

# ── Helper: merge SMILES and is_active ───────────────────────────────────────
def add_meta(df, id_col):
    df = df.rename(columns={id_col: 'ligand_id'})
    df = df.merge(smiles_df[['ligand_id', 'SMILES', 'is_active']],
                  on='ligand_id', how='left')
    df['is_active'] = df['is_active'].fillna(0).astype(int)
    return df

# ── Helper: rename cluster columns ───────────────────────────────────────────
def rename_int_clusters(df):
    rename_map = {str(i): f'cluster{i:02d}_BE' for i in range(10)}
    return df.rename(columns=rename_map)

def rename_str_clusters(df):
    rename_map = {f'{i:02d}': f'cluster{i:02d}_CS' for i in range(10)}
    return df.rename(columns=rename_map)

# ════════════════════════════════════════════════════════════════════════════
# GLIDE
# ════════════════════════════════════════════════════════════════════════════
print("\n── Glide ──")

df = pd.read_csv(f'{GLIDE_DIR}BE_framewise_scores.csv')
df = rename_int_clusters(df)
df = add_meta(df, 'NAME')
cols = ['ligand_id', 'SMILES'] + [f'cluster{i:02d}_BE' for i in range(10)] + ['is_active']
df[cols].to_csv(f'{OUT_DIR}CCR2_Ensemble_Glide_framewise_scores.csv', index=False)
print(f"  Framewise: {len(df)} compounds ({df['is_active'].sum()} actives)")

df = pd.read_csv(f'{GLIDE_DIR}BE_min_ranked.csv')
df = add_meta(df, 'NAME')
df[['ligand_id', 'SMILES', 'BE_min', 'is_active']].to_csv(
    f'{OUT_DIR}CCR2_Ensemble_Glide_BEmin_ranked.csv', index=False)
print(f"  BEmin: {len(df)} compounds ({df['is_active'].sum()} actives)")

df = pd.read_csv(f'{GLIDE_DIR}BE_avg_ranked.csv')
df = add_meta(df, 'NAME')
df[['ligand_id', 'SMILES', 'BE_avg', 'is_active']].to_csv(
    f'{OUT_DIR}CCR2_Ensemble_Glide_BEavg_ranked.csv', index=False)
print(f"  BEavg: {len(df)} compounds ({df['is_active'].sum()} actives)")

# ════════════════════════════════════════════════════════════════════════════
# VINA
# ════════════════════════════════════════════════════════════════════════════
print("\n── Vina ──")

df = pd.read_csv(f'{VINA_DIR}BE_framewise_scores.csv')
df = rename_int_clusters(df)
df = add_meta(df, 'BaseName')
cols = ['ligand_id', 'SMILES'] + [f'cluster{i:02d}_BE' for i in range(10)] + ['is_active']
df[cols].to_csv(f'{OUT_DIR}CCR2_Ensemble_Vina_framewise_scores.csv', index=False)
print(f"  Framewise: {len(df)} compounds ({df['is_active'].sum()} actives)")

df = pd.read_csv(f'{VINA_DIR}BE_min_ranked.csv')
df = add_meta(df, 'BaseName')
df[['ligand_id', 'SMILES', 'BE_min', 'is_active']].to_csv(
    f'{OUT_DIR}CCR2_Ensemble_Vina_BEmin_ranked.csv', index=False)
print(f"  BEmin: {len(df)} compounds ({df['is_active'].sum()} actives)")

df = pd.read_csv(f'{VINA_DIR}BE_avg_ranked.csv')
df = add_meta(df, 'BaseName')
df[['ligand_id', 'SMILES', 'BE_avg', 'is_active']].to_csv(
    f'{OUT_DIR}CCR2_Ensemble_Vina_BEavg_ranked.csv', index=False)
print(f"  BEavg: {len(df)} compounds ({df['is_active'].sum()} actives)")

print("\nDone — all ensemble BE CSVs written.")
