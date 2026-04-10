#!/usr/bin/env python3
import json
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import argparse
import os
import sys

# Setup FLAF
if "ANALYSIS_PATH" not in os.environ:
    raise RuntimeError("ANALYSIS_PATH non settata.")
sys.path.append(os.environ["ANALYSIS_PATH"])

from Studies.MuonWPIso_studies.utils import *

hep.style.use("CMS")

parser = argparse.ArgumentParser()
parser.add_argument("--year", required=True)
parser.add_argument("--treeName", default="singleMu")
args = parser.parse_args()

period = f"Run3_{args.year}"

formatted_json_path = efficiency_path.format(period=period)
formatted_plot_path = os.path.join(formatted_json_path, f"plots_{args.treeName}")
os.makedirs(formatted_plot_path, exist_ok=True)

input_json = f"{formatted_json_path}/results_{args.treeName}.json"

# with open(input_json) as f:
#     data = json.load(f)

for muID in ID_WPs.keys():
    for muIso in Iso_WPs.keys():
        with open(input_json) as f:
            data = json.load(f)
        filtered_data = [
            e for e in data
            if muID in e["mu1_ID"]
            and muIso in e["mu1_Iso"]
            and muID in e["mu2_ID"]
            and muIso in e["mu2_Iso"]
        ]

        # 3. Ordinamento dei dati filtrati per significatività
        data = sorted(filtered_data, key=lambda e: e["s_sqrtB"]["value"], reverse=True)

        # 4. Taglio per leggibilità (Top 80 dopo il filtro)
        print(f"Entries totali: {len(data)}")
        if len(data) > 80:
            print(f"len exceeds 80: {len(data)}, thus reducing it to 80..")
            data = data[:80]

        # 5. Creazione Labels (Corretta con le tue chiavi del dizionario)
        # Usiamo \n per mandare a capo le info del secondo muone
        # Usa le chiavi corrette del tuo dict: mu1_ID e mu1_Iso
        labels = [f"mu1: ID = {e['mu1_ID']}, Iso = {e['mu1_Iso'].split("_")[-1]} \n mu2: ID = {e['mu2_ID']}, Iso = {e['mu2_Iso'].split("_")[-1]}" for e in data]
        x = np.arange(len(data))

        def quick_plot(y, yerr, ylabel, name):
            fig_width = max(15, len(x) * 0.4)
            fig, ax = plt.subplots(figsize=(fig_width, 8))
            ax.errorbar(x, y, yerr=yerr, fmt="o", color="black", ecolor="red", capsize=2)
            ax.set_xticks(x)
            ax.set_xlabel(f"Iso={muIso}")
            ax.set_xticklabels(labels, rotation=45, fontsize=15)
            ax.set_ylabel(ylabel)
            # ax.set_title(muIso)
            # ax.text(
            #         0.95,0.95, muIso, transform=ax.transAxes,
            #         fontsize=25, ha="right", va="top"
            #     )
            hep.cms.label('Preliminary', data=False,  year=args.year if args.year != "all" else "2022-2023", lumi=lumi_dict[args.year if args.year != "all" else "2022-2023"], ax=ax, loc=0, com="13.6",)
            plt.savefig(f"{formatted_plot_path}/{name}_{muID}_{muIso}_{args.treeName}.png", bbox_inches='tight')
            plt.close()

        # 1. Plot S/sqrt(B)
        vals = np.array([e["s_sqrtB"]["value"] for e in data])
        errs = np.array([e["s_sqrtB"]["error"] for e in data])
        quick_plot(vals, errs, r"$S/\sqrt{B}$", "s_sqrtB", )

        # 2. Plot Signal Efficiency
        vals = np.array([e["Signal"]["eff"] for e in data])
        errs = np.array([e["Signal"]["err"] for e in data])
        quick_plot(vals, errs, "Signal Efficiency", "eff_signal")

        # 3. Stack Plot dei Background (Yields)
        fig, ax = plt.subplots(figsize=(max(15, len(x)*0.4), 8))
        bottom = np.zeros(len(x))

        # Assicurati che le chiavi corrispondano esattamente a quelle del tuo JSON (es. P maiuscola)
        categories = ["Prompt",  "Signal", "NonPrompt", "TTbar",]
        colors = ["dodgerblue","cyan", "orange", "forestgreen", ]

        for cat, color in zip(categories, colors):
            # ESTRAZIONE CORRETTA: e[cat]["yield"] invece di e[cat]
            yields = np.array([e[cat]["yield"] for e in data])

            ax.bar(x, yields, bottom=bottom, label=cat, color=color, alpha=0.7)
            bottom += yields

        ax.set_xticks(x)
        ax.set_xlabel(f"Iso={muIso}")
        ax.set_xticklabels(labels, rotation=45, fontsize=15)
        ax.set_yscale("log")
        ax.set_ylabel("Yields")
        ax.legend()
        # ax.text(
        #         0.95,0.95, muIso, transform=ax.transAxes,
        #         fontsize=25, ha="right", va="top"
        #     )
        hep.cms.label('Preliminary', data=False,  year=args.year if args.year != "all" else "2022-2023", lumi=lumi_dict[args.year if args.year != "all" else "2022-2023"], ax=ax, loc=0, com="13.6",)

        plt.savefig(f"{formatted_plot_path}/yields_stack_{muID}_{muIso}.png", bbox_inches='tight')
        print(f"--- Stack plot salvato in: {formatted_plot_path}/yields_stack_{muID}_{muIso}.png")
        plt.close()
        print(f"--- Plot pronti in: {formatted_plot_path}")
