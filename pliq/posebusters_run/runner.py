# -*- coding: utf-8 -*-
"""Run PoseBusters (redock) over CSV cases and save results to CSV."""
import os
import pandas as pd
from pathlib import Path
from tqdm import tqdm

try:
    from posebusters import PoseBusters
except ImportError:
    PoseBusters = None

try:
    from rdkit.Chem.rdmolfiles import MolFromPDBFile
except ImportError:
    MolFromPDBFile = None


def run_posebusters(
    csv_file,
    path_ref_lig,
    path_docked,
    out_csv,
    n_models=20,
    n_cases=0,
    csv_sep="\t",
    subfolder_col="Subfolder",
    use_mol2=True,
    fallback_pdb=True,
):
    """
    Run PoseBusters in redock mode. If use_mol2=True, convert PDB to MOL2 with obabel;
    on failure if fallback_pdb=True, load PDB with RDKit (sanitize=False) and pass Mol objects.
    """
    if PoseBusters is None:
        raise ImportError("posebusters is required. Install with: pip install posebusters")

    df1 = pd.read_csv(csv_file, sep=csv_sep)
    if n_cases and n_cases > 0:
        df1 = df1.iloc[: n_cases]
    df1 = df1.copy()
    if subfolder_col in df1.columns and "pdb_id" not in df1.columns:
        df1["pdb_id"] = df1[subfolder_col]

    buster = PoseBusters(config="redock")
    results = []
    list_fail = []
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)

    for _, row in tqdm(df1.iterrows(), total=len(df1), desc="PoseBusters"):
        case_name = row["pdb_id"]
        true_pdb_path = Path(path_ref_lig) / case_name / f"{case_name}_ligand_chainB.pdb"
        docked_folder = Path(path_docked) / f"posebusters_{case_name}"

        for af_model_id in range(n_models):
            pred_pdb_path = docked_folder / f"posebusters_{case_name}_model_{af_model_id}_aligned_ligand.pdb"
            cond_file = docked_folder / f"posebusters_{case_name}_model_{af_model_id}_aligned_protein.pdb"
            pred_mol2_path = pred_pdb_path.with_suffix(".mol2")
            true_mol2_path = true_pdb_path.with_suffix(".mol2")

            if not pred_pdb_path.exists() or not true_pdb_path.exists() or not cond_file.exists():
                list_fail.append({"pdb_id": case_name, "af_model_id": af_model_id, "reason": "file_missing"})
                continue

            ok = False
            last_err = ""
            if use_mol2:
                try:
                    import subprocess
                    subprocess.run(
                        ["obabel", str(pred_pdb_path), "-O", str(pred_mol2_path), "--ph", "7.4"],
                        check=True, capture_output=True, text=True, timeout=30,
                    )
                    subprocess.run(
                        ["obabel", str(true_pdb_path), "-O", str(true_mol2_path), "--ph", "7.4"],
                        check=True, capture_output=True, text=True, timeout=30,
                    )
                    df_result = buster.bust([str(pred_mol2_path)], str(true_mol2_path), str(cond_file), full_report=True)
                    df_result["pdb_id"] = case_name
                    df_result["af_model_id"] = af_model_id
                    results.append(df_result)
                    ok = True
                except Exception as e:
                    last_err = f"mol2: {str(e)[:200]}"
                finally:
                    for p in [pred_mol2_path, true_mol2_path]:
                        if p.exists():
                            try:
                                os.remove(p)
                            except Exception:
                                pass

            if not ok and fallback_pdb and MolFromPDBFile is not None:
                try:
                    mol_pred = MolFromPDBFile(str(pred_pdb_path), sanitize=False, removeHs=True)
                    mol_true = MolFromPDBFile(str(true_pdb_path), sanitize=False, removeHs=True)
                    mol_cond = MolFromPDBFile(str(cond_file), sanitize=False, removeHs=True)
                    if mol_pred is not None and mol_true is not None and mol_cond is not None:
                        df_result = buster.bust([mol_pred], mol_true, mol_cond, full_report=True)
                        df_result["pdb_id"] = case_name
                        df_result["af_model_id"] = af_model_id
                        results.append(df_result)
                        ok = True
                except Exception as e:
                    last_err = f"{last_err}; pdb: {str(e)[:200]}" if last_err else f"pdb: {str(e)[:200]}"

            if not ok:
                list_fail.append({
                    "pdb_id": case_name,
                    "af_model_id": af_model_id,
                    "reason": last_err or "posebusters_failed",
                })

        if results:
            pd.concat(results, ignore_index=True).to_csv(out_csv, index=False)

    if results:
        pd.concat(results, ignore_index=True).to_csv(out_csv, index=False)
    if list_fail:
        fail_csv = Path(out_csv).with_name(Path(out_csv).stem + "_failures.csv")
        pd.DataFrame(list_fail).to_csv(fail_csv, index=False)
        print(f"PoseBusters: {len(list_fail)} failures recorded -> {fail_csv}")
    return results, list_fail
