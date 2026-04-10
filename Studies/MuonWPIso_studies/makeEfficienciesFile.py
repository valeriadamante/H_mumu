import os
import sys
import argparse
import uproot
import numpy as np
import pandas as pd
import json
from itertools import product
import ROOT

if "ANALYSIS_PATH" not in os.environ:
    raise RuntimeError("Devi avere ANALYSIS_PATH settato nell'ambiente.")
sys.path.append(os.environ["ANALYSIS_PATH"])

from Studies.MuonWPIso_studies.utils import *
# Aggiungi questa classe all'inizio dello script
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)



ROOT.EnableThreadSafety()

def get_wp_labels():
    def extract(d):
        names = {}
        for k, v in d.items():
            if isinstance(v, str): names[k] = v
            elif isinstance(v, dict):
                for wp in v.keys(): names[f"{k}_{wp}"] = v[wp]
        return list(names.keys())
    return extract(ID_WPs), extract(Iso_WPs)

def run_multi_sample_analysis(args):
    ids, isos = get_wp_labels()
    wp_combos = [f"{i}-{s}" for i, s in product(ids, isos)]
    n_wp = len(wp_combos)

    storage = {combo: {cat: 0.0 for cat in samples_dict.keys()} for combo in product(wp_combos, wp_combos)}
    storage_w2 = {combo: {cat: 0.0 for cat in samples_dict.keys()} for combo in product(wp_combos, wp_combos)}

    totals_den = {cat: 0.0 for cat in samples_dict.keys()}
    totals_den_w2 = {cat: 0.0 for cat in samples_dict.keys()}

    period = f"Run3_{args.year}"

    formatted_mini_path = miniTuple_path.format(period=period)
    formatted_efficiency_path = efficiency_path.format(period=period)
    print(formatted_efficiency_path)
    formatted_table_path = table_path.format(period=period)
    os.makedirs(formatted_efficiency_path, exist_ok=True)
    os.makedirs(formatted_table_path, exist_ok=True)
    for category, samples in samples_dict.items():
        for sample in samples:
            infile = f"{formatted_mini_path}/mini_{sample}.root"
            if not os.path.exists(infile): continue

            print(f"--- Loading {sample} ({category})...")
            with uproot.open(infile) as f:
                data = f[args.treeName].arrays(library="np")


            w = data["weight_base"]
            w2 = w**2
            totals_den[category] += w.sum()
            totals_den_w2[category] += w2.sum()

            masks_mu1 = {wp: (data[f"mu1_{wp.split('-')[0]}"] & data[f"mu1_{'_'.join(wp.split('-')[1:])}"]) for wp in wp_combos}
            masks_mu2 = {wp: (data[f"mu2_{wp.split('-')[0]}"] & data[f"mu2_{'_'.join(wp.split('-')[1:])}"]) for wp in wp_combos}

            for wp1 in wp_combos:
                m1 = masks_mu1[wp1]
                for wp2 in wp_combos:
                    m2 = masks_mu2[wp2]
                    mask = m1 & m2
                    storage[(wp1, wp2)][category] += w[mask].sum()
                    storage_w2[(wp1, wp2)][category] += w2[mask].sum()

    # Creazione Output ROOT per istogrammi 2D
    def create_h2(name, title):
        h = ROOT.TH2D(name, title, n_wp, 0.5, n_wp+0.5, n_wp, 0.5, n_wp+0.5)
        for k, label in enumerate(wp_combos, 1):
            h.GetXaxis().SetBinLabel(k, label)
            h.GetYaxis().SetBinLabel(k, label)
        return h

    # Efficienze e Significatività
    h_s_sqrtb = create_h2("h_s_sqrtb", "S/#sqrt{B}")
    h_eff_sig = create_h2("h_eff_sig", "Signal Efficiency")
    h_eff_bkg_p = create_h2("h_eff_bkg_p", "Prompt Bkg Efficiency")
    h_eff_bkg_np = create_h2("h_eff_bkg_np", "NonPrompt Bkg Efficiency")
    h_eff_bkg_tt = create_h2("h_eff_bkg_tt", "TT Bkg Efficiency")
    h_eff_bkg_tot = create_h2("h_eff_bkg_tot", "Total Bkg Efficiency")

    # Numeratori
    h_num_sig = create_h2("h_num_sig", "Signal Yield")
    h_num_bkg_p = create_h2("h_num_bkg_p", "Prompt Bkg Yield")
    h_num_bkg_np = create_h2("h_num_bkg_np", "NonPrompt Bkg Yield")
    h_num_bkg_tt = create_h2("h_num_bkg_tt", "TT Bkg Yield")
    h_num_bkg_tot = create_h2("h_num_bkg_tot", "Total Bkg Yield")

    # Denominatori
    h_den_sig = create_h2("h_den_sig", "Signal Yield")
    h_den_bkg_p = create_h2("h_den_bkg_p", "Prompt Bkg Yield")
    h_den_bkg_np = create_h2("h_den_bkg_np", "NonPrompt Bkg Yield")
    h_den_bkg_tt = create_h2("h_den_bkg_tt", "TT Bkg Yield")
    h_den_bkg_tot = create_h2("h_den_bkg_tot", "Total Bkg Yield")

    results_full = []
    results_minimal = []

    for i, wp1 in enumerate(wp_combos, start=1):
        h_s_sqrtb.GetXaxis().SetBinLabel(i, wp1)
        for j, wp2 in enumerate(wp_combos, start=1):
            yields = storage[(wp1, wp2)]
            w2 = storage_w2[(wp1, wp2)]

            sig = yields["Signal"]
            sig_err = np.sqrt(w2["Signal"])
            h_num_sig.SetBinContent(i, j, sig)
            h_num_sig.SetBinError(i, j, sig_err)

            bkg_p = yields["Prompt"]
            bkg_p_err = np.sqrt(w2["Prompt"])
            h_num_bkg_p.SetBinContent(i, j, bkg_p)
            h_num_bkg_p.SetBinError(i, j, bkg_p_err)

            bkg_np = yields["NonPrompt"]
            bkg_np_err = np.sqrt(w2["NonPrompt"])
            h_num_bkg_np.SetBinContent(i, j, bkg_np)
            h_num_bkg_np.SetBinError(i, j, bkg_np_err)

            bkg_tt =  yields["TTbar"]
            bkg_tt_err = np.sqrt(w2["TTbar"])
            h_num_bkg_tt.SetBinContent(i, j, bkg_tt)
            h_num_bkg_tt.SetBinError(i, j, bkg_tt_err)

            bkg_tot = yields["Prompt"] + yields["NonPrompt"] + yields["TTbar"]
            bkg_tot_err = np.sqrt(w2["Prompt"] + w2["NonPrompt"] + w2["TTbar"])
            h_num_bkg_tot.SetBinContent(i, j, bkg_tot)
            h_num_bkg_tot.SetBinError(i, j, bkg_tot_err)



            def fill_eff(h, num, num_w2, den, den_w2):
                eff = num / den if den > 0 else 0
                err = np.sqrt(num_w2) / den if den > 0 else 0 # Errore approssimato (den err trascurabile)
                h.SetBinContent(i, j, eff)
                h.SetBinError(i, j, err)
                return eff, err

            den_sig = totals_den["Signal"]
            den_sig_err = np.sqrt(totals_den_w2["Signal"])
            h_den_sig.SetBinContent(i, j, den_sig)
            h_den_sig.SetBinError(i, j, den_sig_err)
            den_bkg_p = totals_den["Prompt"]
            den_bkg_p_err = np.sqrt(totals_den_w2["Prompt"])
            h_den_bkg_p.SetBinContent(i, j, den_bkg_p)
            h_den_bkg_p.SetBinError(i, j, den_bkg_p_err)
            den_bkg_np = totals_den["NonPrompt"]
            den_bkg_np_err = np.sqrt(totals_den_w2["NonPrompt"])
            h_den_bkg_np.SetBinContent(i, j, den_bkg_np)
            h_den_bkg_np.SetBinError(i, j, den_bkg_np_err)
            den_bkg_tt =  totals_den["TTbar"]
            den_bkg_tt_err = np.sqrt(totals_den_w2["TTbar"])
            h_den_bkg_tt.SetBinContent(i, j, den_bkg_tt)
            h_den_bkg_tt.SetBinError(i, j, den_bkg_tt_err)
            den_bkg_tot = totals_den["Prompt"] + totals_den["NonPrompt"] + totals_den["TTbar"]
            den_bkg_tot_err = np.sqrt(totals_den_w2["Prompt"] + totals_den_w2["NonPrompt"] + totals_den_w2["TTbar"])
            h_den_bkg_tot.SetBinContent(i, j, den_bkg_tot)
            h_den_bkg_tot.SetBinError(i, j, den_bkg_tot_err)


            e_s, e_s_err = fill_eff(h_eff_sig, sig, w2["Signal"], den_sig, totals_den_w2["Signal"])
            e_p, e_p_err = fill_eff(h_eff_bkg_p, yields["Prompt"], w2["Prompt"], totals_den["Prompt"], totals_den_w2["Prompt"])
            e_np, e_np_err = fill_eff(h_eff_bkg_np, yields["NonPrompt"], w2["NonPrompt"], totals_den["NonPrompt"], totals_den_w2["NonPrompt"])
            e_tt, e_tt_err = fill_eff(h_eff_bkg_tt, yields["TTbar"], w2["TTbar"], totals_den["TTbar"], totals_den_w2["TTbar"])
            e_bkg_tot, e_bkg_tot_err = fill_eff(h_eff_bkg_tot, bkg_tot, (bkg_tot_err**2), den_bkg_tot, 0)


            s_sqrtb = sig / np.sqrt(bkg_tot) if bkg_tot > 0 else 0
            s_sqrtbkg_err = 0
            if s_sqrtb > 0:
                s_sqrtbkg_err = s_sqrtb * np.sqrt((sig_err/sig)**2 + (bkg_tot_err/(2*bkg_tot))**2)
            h_s_sqrtb.SetBinContent(i, j, s_sqrtb)
            h_s_sqrtb.SetBinError(i, j, s_sqrtbkg_err)

            s_sqrtb_prompt = sig / np.sqrt(bkg_p) if bkg_p > 0 else 0
            s_sqrtbkg_prompt_err = 0
            if s_sqrtb_prompt > 0:
                s_sqrtbkg_prompt_err = s_sqrtbkg_prompt_err * np.sqrt((sig_err/sig)**2 + (bkg_p_err/(2*bkg_p))**2)

            s_sqrtb_nonprompt = sig / np.sqrt(bkg_np) if bkg_np > 0 else 0
            s_sqrtbkg_nonprompt_err = 0
            if s_sqrtb_nonprompt > 0:
                s_sqrtbkg_nonprompt_err = s_sqrtbkg_nonprompt_err * np.sqrt((sig_err/sig)**2 + (bkg_np_err/(2*bkg_np))**2)


            s_sqrtb_tt = sig / np.sqrt(bkg_tt) if bkg_tt > 0 else 0
            s_sqrtbkg_tt_err = 0
            if s_sqrtb_tt > 0:
                s_sqrtbkg_tt_err = s_sqrtbkg_tt_err * np.sqrt((sig_err/sig)**2 + (bkg_tt_err/(2*bkg_tt))**2)


            results_minimal.append({
                "mu1_ID": wp1.split('-')[0],
                "mu1_Iso": wp1.split('-')[1],
                "mu2_ID": wp2.split('-')[0],
                "mu2_Iso": wp2.split('-')[1],
                "s_sqrtB": {"value": s_sqrtb},
                "s_sqrtB_Prompt": {"value": s_sqrtb_prompt},
                "s_sqrtB_NonPrompt": {"value": s_sqrtb_nonprompt},
                "s_sqrtB_TTbar": {"value": s_sqrtb_tt},
                "Signal": {"eff": e_s},
                "Background": {"eff": e_bkg_tot},
                "Prompt": {"eff": e_p},
                "NonPrompt": {"eff": e_np},
                "TTbar": {"eff": e_tt},
            })

            results_full.append({
                "mu1_ID": wp1.split('-')[0],
                "mu1_Iso": wp1.split('-')[1],
                "mu2_ID": wp2.split('-')[0],
                "mu2_Iso": wp2.split('-')[1],
                "s_sqrtB": {"value": s_sqrtb, "error": s_sqrtbkg_err},
                "s_sqrtB_Prompt": {"value": s_sqrtb_prompt, "error": s_sqrtbkg_prompt_err},
                "s_sqrtB_NonPrompt": {"value": s_sqrtb_nonprompt, "error": s_sqrtbkg_nonprompt_err},
                "s_sqrtB_TTbar": {"value": s_sqrtb_tt, "error": s_sqrtbkg_tt_err},
                "Signal": {"yield": sig, "yield_err":sig_err ,"eff": e_s, "err": e_s_err, "den": den_sig,  "den_err": den_sig_err,},
                "Background": {"yield": bkg_tot, "yield_err":bkg_tot_err ,"eff": e_bkg_tot, "err": e_bkg_tot_err, "den": den_bkg_tot,  "den_err": den_bkg_tot_err,},
                "Prompt": {"yield": bkg_p, "yield_err":bkg_p_err ,"eff": e_p, "err": e_p_err, "den": den_bkg_p,  "den_err": den_bkg_p_err,},
                "NonPrompt": {"yield": bkg_np, "yield_err":bkg_np_err ,"eff": e_np, "err": e_np_err, "den": den_bkg_np,  "den_err": den_bkg_np_err,},
                "TTbar": {"yield": bkg_tt, "yield_err":bkg_tt_err ,"eff": e_tt, "err": e_tt_err, "den": den_bkg_tt,  "den_err": den_bkg_tt_err,},
            })

    # --- Salvataggio File ---
    # 1. JSON & CSV
    with open(f"{formatted_efficiency_path}/results_{args.treeName}.json", "w") as f:
        json.dump(results_full, f, indent=4, cls=NpEncoder)
    df_full = pd.json_normalize(results_full, sep='_')
    # --- Estrazione Best WPs ---
    metrics = {
        "Massimo S/sqrt(B)": "s_sqrtB_value",
        "Massima Signal Efficiency": "Signal_eff",
        "Massimo S/sqrt(B_prompt)": "s_sqrtB_Prompt_value"
    }

    print("\n" + "="*50)
    print("RISULTATI OTTIMIZZAZIONE WP")
    print("="*50)

    for title, col in metrics.items():
        if col in df_full.columns:
            # Trova la riga con il valore massimo per la colonna specificata
            best_row = df_full.loc[df_full[col].idxmax()]

            print(f"\n>>> {title}:")
            print(f"    Valore: {best_row[col]:.5f}")
            print(f"    Muone 1: ID = {best_row['mu1_ID']}, Iso = {best_row['mu1_Iso']}")
            print(f"    Muone 2: ID = {best_row['mu2_ID']}, Iso = {best_row['mu2_Iso']}")
        else:
            print(f"\n[!] Errore: Colonna {col} non trovata nel DataFrame.")

    print("\n" + "="*50)
    df_minimal = pd.json_normalize(results_minimal, sep='_')
    target_col = 's_sqrtB_value'
    if target_col in df_full.columns:
        print(f"reordering df full by {target_col}")
        df_full = df_full.sort_values(by=target_col, ascending=False)

    if target_col in df_minimal.columns:
        print(f"reordering df minimal by {target_col}")
        df_minimal = df_minimal.sort_values(by=target_col, ascending=False)

    df_full.to_csv(f"{formatted_table_path}/results_{args.treeName}.csv", index=False)
    df_minimal.to_csv(f"{formatted_table_path}/results_{args.treeName}_CBIds_IsoIds.tsv", sep="\t", index=False)
    # 2. ROOT Efficienze
    feff = ROOT.TFile(f"{formatted_efficiency_path}/efficiencies_{period}_{args.treeName}_CBIds_IsoIds.root", "RECREATE")
    for h in [h_s_sqrtb, h_eff_sig, h_eff_bkg_p, h_eff_bkg_np, h_eff_bkg_tt, h_eff_bkg_tot]: h.Write()
    feff.Close()

    # 3. ROOT Numeratori
    fnum = ROOT.TFile(f"{formatted_efficiency_path}/numerators_{period}_{args.treeName}_CBIds_IsoIds.root", "RECREATE")
    for h in [h_num_sig, h_num_bkg_p, h_num_bkg_np, h_num_bkg_tt, h_num_bkg_tot]: h.Write()
    fnum.Close()

    # 4. ROOT Denominatori
    fden = ROOT.TFile(f"{formatted_efficiency_path}/denumerators_{period}_{args.treeName}_CBIds_IsoIds.root", "RECREATE")
    for h in [h_den_sig, h_den_bkg_tot]: h.Write()
    fden.Close()

    print(f"\n>>> Elaborazione completata per {args.year}. File salvati in {formatted_efficiency_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", required=True)
    parser.add_argument("--treeName", default="singleMu")
    args = parser.parse_args()
    run_multi_sample_analysis(args)


# #!/usr/bin/env python3
# import uproot
# import numpy as np
# import ROOT
# import os, argparse
# import pandas as pd
# from itertools import product
# import sys

# # Setup ambiente FLAF
# if "ANALYSIS_PATH" not in os.environ:
#     raise RuntimeError("Devi avere ANALYSIS_PATH settato nell'ambiente.")
# sys.path.append(os.environ["ANALYSIS_PATH"])

# from Studies.MuonWPIso_studies.utils import *

# ROOT.EnableThreadSafety()

# parser = argparse.ArgumentParser()
# parser.add_argument("--year", required=True)
# parser.add_argument("--samples", required=True, help="Campioni separati da virgola (es. VBFHto2Mu,DY)")
# parser.add_argument("--treeName", default='singleMu')
# args = parser.parse_args()

# # Percorsi dinamici da utils.py
# period = f"Run3_{args.year}"
# efficiency_path = efficiency_path.format(period=period)
# table_path = table_path.format(period=period)
# miniTuple_path = miniTuple_path.format(period=period)

# os.makedirs(efficiency_path, exist_ok=True)
# os.makedirs(table_path, exist_ok=True)

# def getWPLabels():
#     wp_dict_names = {}
#     iso_dict_names = {}
#     for idCut_name, idCut_def in ID_WPs.items():
#         if isinstance(idCut_def, str): wp_dict_names[idCut_name] = idCut_def
#         elif isinstance(idCut_def, dict):
#             for wp, def_ in idCut_def.items(): wp_dict_names[f"{idCut_name}_{wp}"] = def_
#     for isoCut_name, isoCut_def in Iso_WPs.items():
#         if isinstance(isoCut_def, str): iso_dict_names[isoCut_name] = isoCut_def
#         elif isinstance(isoCut_def, dict):
#             for wp, def_ in isoCut_def.items(): iso_dict_names[f"{isoCut_name}_{wp}"] = def_
#     return wp_dict_names, iso_dict_names

# samples = args.samples.split(",")
# wp_dict_names, iso_dict_names = getWPLabels()
# ids = list(wp_dict_names.keys())
# isos = list(iso_dict_names.keys())
# combinations = list(product(ids, isos))
# n_wp = len(combinations)

# for sample in samples:
#     infile = f"{miniTuple_path}/mini_{sample}.root"
#     if not os.path.exists(infile):
#         print(f"File {infile} non trovato, salto...")
#         continue

#     with uproot.open(infile) as f:
#         data = f[args.treeName].arrays(library="np")

#     fout = ROOT.TFile(f"{efficiency_path}/eff_{sample}_{args.treeName}.root", "RECREATE")

#     # Metadati
#     total_weight = data["weight_base"].sum()
#     sumw2_total = np.sum(data["weight_base"]**2)

#     h_eff = ROOT.TH2D("h_eff", ";mu1 ID+Iso WP;mu2 ID+Iso WP", n_wp, 0.5, n_wp + 0.5, n_wp, 0.5, n_wp + 0.5)
#     results = []

#     for i, (id1, iso1) in enumerate(combinations, start=1):
#         mask1 = data[f"mu1_{id1}"] & data[f"mu1_{iso1}"]
#         label1 = f"{id1}_{iso1}"
#         for j, (id2, iso2) in enumerate(combinations, start=1):
#             mask2 = data[f"mu2_{id2}"] & data[f"mu2_{iso2}"]
#             label2 = f"{id2}_{iso2}"

#             w = data["weight_base"][mask1 & mask2]
#             yield_pass = w.sum()
#             eff = yield_pass / max(total_weight, 1e-9)
#             err = np.sqrt(np.sum(w**2)) / max(total_weight, 1e-9)

#             h_eff.SetBinContent(i, j, eff)
#             h_eff.SetBinError(i, j, err)
#             if j == 1: h_eff.GetXaxis().SetBinLabel(i, label1)
#             if i == 1: h_eff.GetYaxis().SetBinLabel(j, label2)

#             results.append({
#                 "mu1_WP": label1, "mu2_WP": label2,
#                 "Eff": eff, "Yield": yield_pass, "Total": total_weight
#             })

#     pd.DataFrame(results).to_csv(f"{table_path}/eff_{sample}.csv", sep="\t", index=False)
#     h_eff.Write()
#     fout.Close()
#     print(f"Completato: {sample}")