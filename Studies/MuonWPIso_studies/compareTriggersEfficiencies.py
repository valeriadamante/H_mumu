import os
import sys
import argparse
import uproot
import numpy as np
import pandas as pd
import json
from itertools import product
import ROOT

# Configurazione ambiente e utility
if "ANALYSIS_PATH" not in os.environ:
    raise RuntimeError("Devi avere ANALYSIS_PATH settato nell'ambiente.")
sys.path.append(os.environ["ANALYSIS_PATH"])

# Importiamo i path e i dizionari dai tuoi utils
from Studies.MuonWPIso_studies.utils import *

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.floating): return float(obj)
        if isinstance(obj, np.integer): return int(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        return super(NpEncoder, self).default(obj)

def get_wp_labels():
    def extract(d):
        names = {}
        for k, v in d.items():
            if isinstance(v, str): names[k] = v
            elif isinstance(v, dict):
                for wp in v.keys(): names[f"{k}_{wp}"] = v[wp]
        return list(names.keys())
    return extract(ID_WPs), extract(Iso_WPs)

def analyze_tree(infile, tree_name, wp_combos):
    print(f"--- Processing Tree: {tree_name}")
    with uproot.open(infile) as f:
        if tree_name not in f:
            return None, None
        data = f[tree_name].arrays(library="np")

    w = data["weight_base"]
    den = w.sum()

    # Pre-calcolo maschere per efficienza
    masks_mu1 = {wp: (data[f"mu1_{wp.split('-')[0]}"] & data[f"mu1_{'_'.join(wp.split('-')[1:])}"]) for wp in wp_combos}
    masks_mu2 = {wp: (data[f"mu2_{wp.split('-')[0]}"] & data[f"mu2_{'_'.join(wp.split('-')[1:])}"]) for wp in wp_combos}

    results = {}
    for wp1, wp2 in product(wp_combos, wp_combos):
        mask = masks_mu1[wp1] & masks_mu2[wp2]
        num = w[mask].sum()
        results[(wp1, wp2)] = num / den if den > 0 else 0

    return results, den

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True, help="Nome del sample (es. DYJets o Signal)")
    parser.add_argument("--year", required=True, help="Anno (es. 2022, 2022EE, 2023)")
    args = parser.parse_args()

    # Costruzione path standard
    period = f"Run3_{args.year}"
    formatted_mini_path = miniTuple_path.format(period=period)
    infile = f"{formatted_mini_path}/mini_{args.sample}.root"

    if not os.path.exists(infile):
        print(f"Errore: Il file {infile} non esiste.")
        sys.exit(1)

    print(f"Analizzando file: {infile}")
    # Creazione della tabella di confronto
    comparison_results = []

    ids, isos = get_wp_labels()
    wp_combos = [f"{i}-{s}" for i, s in product(ids, isos)]

    trees_to_compare = ["singleMu", "singleMuOrDiMuon"]
    tree_data = {}
    trees_for_denum = []

    for t in trees_to_compare:
        effs, den = analyze_tree(infile, t, wp_combos)
        if effs is not None:
            tree_data[t] = effs
        else:
            print(f"Attenzione: Tree {t} non trovato nel file.")
        with uproot.open(infile) as f:
            trees_for_denum.append(f[t].arrays(library="np"))

    t1, t2 = trees_to_compare[0], trees_to_compare[1]


    denum_01 = trees_for_denum[0][("initial_weights")][0]
    denum_11 = trees_for_denum[0][("Pass_OS_JetTagSel_SignalRegion_weights")][0]
    denum_21 = trees_for_denum[0][("Pass_TRG_weights")][0]


    denum_02 = trees_for_denum[1][("initial_weights")][0]
    denum_12 = trees_for_denum[1][("Pass_OS_JetTagSel_SignalRegion_weights")][0]
    denum_22 = trees_for_denum[1][("Pass_TRG_weights")][0]

    comparison_results.append({
            "mu1_WP": 'no requirements, initial events',
            "mu2_WP": 'no requirements, initial events',
            f"yield {t1}": denum_01,
            f"yield {t2}": denum_02,
            "diff": denum_02 - denum_01,
            "rel gain (%)": ((denum_02 - denum_01) / denum_01 * 100) if denum_01 > 0 else 0
        })
    comparison_results.append({
            "mu1_WP": 'no requirements, pass baseline',
            "mu2_WP": 'no requirements, pass baseline',
            f"yield {t1}": denum_11,
            f"yield {t2}": denum_12,
            "diff": denum_12 - denum_11,
            "rel gain (%)": ((denum_12 - denum_11) / denum_11 * 100) if denum_11 > 0 else 0
        })

    comparison_results.append({
            "mu1_WP": 'no requirements, pass TRG',
            "mu2_WP": 'no requirements, pass TRG',
            f"yield {t1}": denum_21,
            f"yield {t2}": denum_22,
            "diff": denum_22 - denum_21,
            "rel gain (%)": ((denum_22 - denum_21) / denum_21 * 100) if denum_21 > 0 else 0
        })

    if len(tree_data) < 2:
        print("Impossibile confrontare: mancano i tree necessari.")
        sys.exit(1)


    for wp1, wp2 in product(wp_combos, wp_combos):
        e1 = tree_data[t1][(wp1, wp2)]
        e2 = tree_data[t2][(wp1, wp2)]

        # Guadagno relativo
        rel_gain = ((e2 - e1) / e1 * 100) if e1 > 0 else 0

        comparison_results.append({
            "mu1_WP": wp1,
            "mu2_WP": wp2,
            f"eff_{t1}": e1,
            f"eff_{t2}": e2,
            "diff": e2 - e1,
            "rel gain (%)": rel_gain
        })

    # Salvataggio nel path delle tabelle definito in utils
    formatted_table_path = table_path.format(period=period)
    os.makedirs(formatted_table_path, exist_ok=True)

    df = pd.DataFrame(comparison_results)
    output_file = f"{formatted_table_path}/triggers_comparison_{args.sample}_{args.year}.csv"
    df.to_csv(output_file, index=False)

    print(f"\n>>> Confronto completato con successo!")
    print(f"Risultati salvati in: {output_file}")