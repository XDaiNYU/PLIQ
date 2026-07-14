#!/usr/bin/env python3
"""Run the bundled 7ZU2_DHT model 15 worked example (four PDBs → scored CSV)."""
from __future__ import annotations

from pathlib import Path

from pliq import run_pliq_from_pdbs

HERE = Path(__file__).resolve().parent

df = run_pliq_from_pdbs(
    HERE / "ref" / "7ZU2_DHT_ligand_chainB.pdb",
    HERE / "ref" / "7ZU2_DHT_protein_renum.pdb",
    HERE / "dock" / "posebusters_7ZU2_DHT_model_15_aligned_ligand.pdb",
    HERE / "dock" / "posebusters_7ZU2_DHT_model_15_aligned_protein.pdb",
    HERE / "pliq_result_full.csv",
    pdb_id="7ZU2_DHT",
    af_model_id=15,
)

summary_cols = [
    "pdb_id", "af_model_id",
    "run_otmol_ok", "run_binana_ok", "run_posebusters_ok", "run_tmscore_ok",
    "ligand_rmsd", "interface_rmsd", "fnat_all", "otmol_map_pairs", "binana_h_mode",
    "binana_hbond_tp_n", "binana_hydrophobic_tp_n", "tm_TMscore", "pb_valid_fraction",
    "pliq", "pliq_eDockQ1", "eDockQ6_1", "eDockQ6_2", "eDockQ6_3", "eDockQ6_4",
    "pliq_term_fnat", "pliq_term_kernel_i", "pliq_term_kernel_l",
    "pliq_term_B6", "pliq_term_H", "pliq_term_B6xH", "pliq_term_posebusters", "pliq_term_TM",
]
summary_cols = [c for c in summary_cols if c in df.columns]
df[summary_cols].to_csv(HERE / "pliq_result_summary.csv", index=False)

print("PLIQ example finished.")
print(f"  full CSV:    {HERE / 'pliq_result_full.csv'} ({df.shape[1]} columns)")
print(f"  summary CSV: {HERE / 'pliq_result_summary.csv'} ({len(summary_cols)} columns)")
print(f"  pliq = {df.iloc[0]['pliq']:.4f}  (default headline score = eDockQ6_1)")
