import pandas as pd

# ── File paths ──────────────────────────────────────────────────────────────
SMILES_FILE     = '/work/users/t/y/tylerdt/Schrodinger/LigPrep/M4R-no_dup-ligprep/M4R-no_dup-combined-fixed_ASD.smi'
GLIDE_FILE      = '/work/users/t/y/tylerdt/Schrodinger/Glide-HTVS-exp/exp-HTVS-best-res_centroid/M4R-ACh-XY6-exp-HTVS/M4R-ACh-XY6-HTVS_pv-best.csv'
VINA_FILE       = '/users/t/y/tylerdt/GPCR_Retrospective_Docking/AutoDockVina/M4R_ACh_XY6/M4R_ACh_XY6_vina_out_final.csv'

# ── Load SMILES file (mixed tab/space separator fix) ─────────────────────────
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

print(f"Total compounds: {len(smiles_df)}")
print(f"Actives: {smiles_df['is_active'].sum()}")
print(f"Decoys: {(smiles_df['is_active'] == 0).sum()}")

# ── Glide ────────────────────────────────────────────────────────────────────
glide = pd.read_csv(GLIDE_FILE)[['NAME', 'r_i_glide_gscore']].rename(columns={
    'NAME': 'ligand_id',
    'r_i_glide_gscore': 'glide_gscore'
})
glide = glide.merge(smiles_df[['ligand_id', 'SMILES', 'is_active']],
                    on='ligand_id', how='left')
glide['is_active'] = glide['is_active'].fillna(0).astype(int)
glide = glide[['ligand_id', 'SMILES', 'glide_gscore', 'is_active']]
glide.to_csv('M4R_PDB_Glide_scores.csv', index=False)
print(f"Glide: {len(glide)} compounds ({glide['is_active'].sum()} actives)")

# ── Vina ─────────────────────────────────────────────────────────────────────
vina = pd.read_csv(VINA_FILE)[['BaseName', 'docking_score']].rename(columns={
    'BaseName': 'ligand_id',
    'docking_score': 'vina_score'
})
vina = vina.merge(smiles_df[['ligand_id', 'SMILES', 'is_active']],
                  on='ligand_id', how='left')
vina['is_active'] = vina['is_active'].fillna(0).astype(int)
vina = vina[['ligand_id', 'SMILES', 'vina_score', 'is_active']]
vina.to_csv('M4R_PDB_Vina_scores.csv', index=False)
print(f"Vina: {len(vina)} compounds ({vina['is_active'].sum()} actives)")

print("\nDone — two clean CSVs written.")

