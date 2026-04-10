#!/usr/bin/env python3
import os
import sys
import ROOT
import argparse
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import re
import importlib

# Setup stile CMS
plt.style.use(hep.style.CMS)

# Setup ambiente FLAF
if "ANALYSIS_PATH" not in os.environ:
    raise RuntimeError("Devi avere ANALYSIS_PATH settato nell'ambiente.")
sys.path.append(os.environ["ANALYSIS_PATH"])

import FLAF.Common.Utilities as Utilities
from FLAF.Common.Setup import Setup
from Studies.MuonWPIso_studies.utils import *

ROOT.EnableThreadSafety()

# ---------------------------------------------------------------------
def get_wp_value(wp_string):
    """Estrae il valore numerico da espressioni tipo 'mu1_var > 0.15'"""
    if not isinstance(wp_string, str): return None
    match = re.search(r'([<>]=?|==)\s*([-+]?\d*\.?\d+)', wp_string)
    if match:
        return float(match.group(2))
    return None

# ---------------------------------------------------------------------
def plot_id_iso_vars(rdf, var_base_name, wp_subdict, out_dir, process_name, year, suffix):
    """
    Crea istogrammi direttamente da RDF senza usare AsNumpy
    """
    col1 = f"mu1_{var_base_name}"
    col2 = f"mu2_{var_base_name}"

    available_cols = [str(c) for c in rdf.GetColumnNames()]
    if col1 not in available_cols or col2 not in available_cols:
        return

    # Determiniamo il range e il binning
    # Per variabili ID (interi) usiamo binning discreto, per Iso (float) range 0-0.5
    if "Id" in var_base_name or "ID" in var_base_name and not "pfRelIso" in var_base_name:
        # Troviamo min/max in modo dinamico da RDF per i bin
        vmin = rdf.Min(col1).GetValue()
        vmax = rdf.Max(col1).GetValue()
        nbins = int(vmax - vmin + 1)
        h_range = (vmin - 0.5, vmax + 0.5)
    else:
        nbins = 50
        h_range = (0.0, 0.5)
    # Creazione dei modelli di istogramma (Histo1D restituisce un puntatore intelligente ROOT)
    # Sintassi: ("nome", "titolo", nbins, xlow, xup)
    model1 = (f"h1_{var_base_name}", "", nbins, h_range[0], h_range[1])
    model2 = (f"h2_{var_base_name}", "", nbins, h_range[0], h_range[1])
    # rdf.Display(col1).Print()
    # rdf.Display(col2).Print()
    # ROOT.EnableImplicitMT()
    # Azioni RDF: ROOT calcolerà questi istogrammi in un unico loop
    h1_ptr = rdf.Histo1D(model1, col1, "weight_base")
    h2_ptr = rdf.Histo1D(model2, col2, "weight_base")

    # Ora passiamo a Matplotlib usando gli istogrammi già pronti
    fig, ax = plt.subplots(figsize=(12, 10))

    # hep.histplot accetta oggetti TH1 di ROOT
    hep.histplot(h1_ptr.GetValue(), ax=ax, label='Muon 1', color='blue', lw=2)
    hep.histplot(h2_ptr.GetValue(), ax=ax, label='Muon 2', color='red', lw=2)

    # Disegno delle linee dei Working Points
    if isinstance(wp_subdict, dict):
        for wp_name, wp_def in wp_subdict.items():
            val = get_wp_value(wp_def)
            if val is not None:
                ax.axvline(val, linestyle='--', color='gray', alpha=0.7)
                # Posizionamento dinamico del testo sopra l'asse
                ax.text(val, ax.get_ylim()[1]*0.6, wp_name, rotation=90, verticalalignment='center', fontsize=10)
    elif isinstance(wp_subdict, str):
        val = get_wp_value(wp_subdict)
        if val is not None:
            ax.axvline(val, linestyle='--', color='gray', alpha=0.7)

    ax.set_yscale('log')
    ax.set_xlabel(var_base_name.replace("_", " "))
    ax.set_ylabel("Weighted Events")
    ax.legend(loc='upper right')

    lumi_dict = {"2022": "7.98", "2022EE": "26.67", "2023": "18.06", "2023BPix": "9.69"}
    lumi = lumi_dict.get(year, "N/A")
    hep.cms.label("Preliminary", data=("Data" in process_name), lumi=lumi, year=year, ax=ax)

    out_name = f"{out_dir}/{process_name}_{var_base_name}_{suffix}.png"
    fig.savefig(out_name, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"--- Salvato: {out_name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", default="2022")
    parser.add_argument("--samples", required=False, default="")
    args = parser.parse_args()
    period = f"Run3_{args.year}"
    setup = Setup.getGlobal(os.environ["ANALYSIS_PATH"], period, "")

    # Importazione dinamica analisi
    analysis_import = setup.global_params["analysis_import"]
    analysis = importlib.import_module(analysis_import)

    formatted_plot_path = plot_path.format(period=period)
    os.makedirs(formatted_plot_path, exist_ok=True)

    # Caricamento campioni tramite GetAnaTuplePathDict
    samples_dict = GetAnaTuplePathDict(period, args.samples)

    for process_name, samples_list in samples_dict.items():
        print(f"\n>>> Processing Process: {process_name}")
        files = expand_filelist(samples_list)
        if not files: continue

        rdf = ROOT.RDataFrame("Events", Utilities.ListToVector(files))

        # Inizializzazione correzioni e selezioni
        Utilities.InitializeCorrections(setup, find_dataset_name(setup.datasets, process_name), stage="HistTuple")
        corrections = Corrections.getGlobal()
        bTagAlgo = analysis.Taggers_branchesNames[setup.global_params.get("bTagAlgo", "particleNet")]
        bTagWPDict = corrections.btag.getWPValues()

        # Selezione finale
        df_singleMu, df_singleAndDiMu, columns = GetAllDfStuff(rdf, [], bTagAlgo, bTagWPDict)

        all_WPs = {**ID_WPs, **Iso_WPs}

        for var_name, wp_subdict in all_WPs.items():
            print(f"Plotting variable: {var_name}")
            plot_id_iso_vars(df_singleMu, var_name, wp_subdict,
                             formatted_plot_path, process_name, args.year, "singleMuOnly")
            plot_id_iso_vars(df_singleAndDiMu, var_name, wp_subdict,
                             formatted_plot_path, process_name, args.year, "singleMuOrDiMu")
# #!/usr/bin/env python3
# import os
# import sys
# import ROOT
# import argparse
# import numpy as np
# import matplotlib.pyplot as plt
# import mplhep as hep
# import re
# import importlib

# # Setup stile CMS
# plt.style.use(hep.style.CMS)

# # Setup ambiente FLAF
# if "ANALYSIS_PATH" not in os.environ:
#     raise RuntimeError("Devi avere ANALYSIS_PATH settato nell'ambiente.")
# sys.path.append(os.environ["ANALYSIS_PATH"])

# import FLAF.Common.Utilities as Utilities
# from FLAF.Common.Setup import Setup
# from Studies.MuonWPIso_studies.utils import *
# from Analysis.H_mumu import * # Importa le funzioni di selezione specifiche

# # Carichiamo le definizioni dei WP e i path dalle tue utility
# from Studies.MuonWPIso_studies.utils import *

# ROOT.EnableThreadSafety()
# ROOT.EnableImplicitMT()

# # ---------------------------------------------------------------------
# def get_wp_value(wp_string):
#     """Estrae il valore numerico da espressioni tipo 'mu1_var > 0.15' o '== 3'"""
#     if not isinstance(wp_string, str): return None
#     match = re.search(r'([<>]=?|==)\s*([-+]?\d*\.?\d+)', wp_string)
#     if match:
#         return float(match.group(2))
#     return None

# # ---------------------------------------------------------------------
# def plot_id_iso_vars(rdf, var_base_name, wp_subdict, out_dir, process_name, year, suffix):
#     """
#     Estrae i dati da RDF e plotta le distribuzioni mu1/mu2 con i WP
#     """
#     col1 = f"mu1_{var_base_name}"
#     col2 = f"mu2_{var_base_name}"

#     # Controlliamo se le colonne esistono nell'RDF
#     available_cols = [str(c) for c in rdf.GetColumnNames()]
#     if col1 not in available_cols or col2 not in available_cols:
#         print(f"[SKIP] Colonne {col1}/{col2} non trovate nell'RDF.")
#         return

#     # Convertiamo in Numpy solo le colonne necessarie per il plot
#     # Nota: Usiamo un filtro per non caricare troppi dati in memoria se non serve
#     data = rdf.AsNumpy(columns=[col1, col2, "weight_base"])

#     fig, ax = plt.subplots(figsize=(12, 10))

#     # Uniamo i valori per determinare il binning
#     all_vals = np.concatenate([data[col1], data[col2]])
#     vmin, vmax = np.nanmin(all_vals), np.nanmax(all_vals)

#     # Binning logico per variabili intere (ID) o continue (Iso)
#     if "Id" in var_base_name and not "pfRelIso" in var_base_name:
#         bins = np.arange(int(vmin)-0.5, int(vmax)+1.5, 1)
#     else:
#         bins = np.linspace(0, min(vmax, 0.5), 50) # Limitiamo il range iso per visibilità

#     # Plot Muon 1 e Muon 2
#     for col, label, color in zip([col1, col2], ['Muon 1', 'Muon 2'], ['blue', 'red']):
#         ax.hist(data[col], bins=bins, weights=data["weight_base"],
#                 histtype='step', lw=2, label=label, color=color)

#     # Disegno delle linee dei Working Points
#     if isinstance(wp_subdict, dict):
#         for wp_name, wp_def in wp_subdict.items():
#             val = get_wp_value(wp_def)
#             if val is not None:
#                 ax.axvline(val, linestyle='--', color='gray', alpha=0.7)
#                 ax.text(val, ax.get_ylim()[1]*0.5, wp_name, rotation=90, verticalalignment='center')
#     elif isinstance(wp_subdict, str):
#         val = get_wp_value(wp_subdict)
#         if val is not None:
#             ax.axvline(val, linestyle='--', color='gray', alpha=0.7)

#     ax.set_yscale('log')
#     ax.set_xlabel(var_base_name.replace("_", " "))
#     ax.set_ylabel("Weighted Events")
#     ax.legend(loc='upper right')

#     # Label CMS
#     lumi = "N/A" # Puoi recuperarla da setup se necessario
#     hep.cms.label("Preliminary", data=False, lumi=lumi, year=year, ax=ax)

#     out_name = f"{out_dir}/{process_name}_{var_base_name}_{suffix}.png"
#     fig.savefig(out_name, dpi=300, bbox_inches='tight')
#     plt.close(fig)
#     print(f"--- Salvato: {out_name}")


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser.add_argument("--year", default="2022")
#     parser.add_argument("--samples", required=False, default="")
#     args = parser.parse_args()
#     period = f"Run3_{args.year}"
#     setup = Setup.getGlobal(os.environ["ANALYSIS_PATH"], period)

#     print(args.samples)
#     print(args.year)


#     # Importazione dinamica del modulo di analisi (es. H_mumu)
#     analysis_import = setup.global_params["analysis_import"]
#     analysis = importlib.import_module(analysis_import)

#     # Path di output
#     formatted_plot_path = plot_path.format(period=period)
#     os.makedirs(formatted_plot_path, exist_ok=True)

#     # Loop sui processi definiti nel framework

#     for process_name,samples_list in GetAnaTuplePathDict(period, args.samples).items():
#         print(f"\n>>> Processing Process: {process_name}")

#         files = expand_filelist(samples_list)
#         print(files)
#         if not files: continue

#         rdf = ROOT.RDataFrame("Events", Utilities.ListToVector(files))

#         # Applichiamo le selezioni base del framework (GetAllDfStuff, selezioni trigger, ecc.)
#         # Nota: GetAllDfStuff di solito restituisce RDF filtrati
#         bTagAlgo = analysis.Taggers_branchesNames[setup.global_params.get("bTagAlgo", "particleNet")]
#         # Qui simuliamo la logica dello script che hai mandato
#         columns = []
#         Utilities.InitializeCorrections(setup, find_dataset_name(setup.datasets, process_name), stage="HistTuple")
#         print("corrections initialized")
#         corrections=Corrections.getGlobal()
#         print("corrections get global")
#         bTagAlgo = analysis.Taggers_branchesNames[
#             setup.global_params.get("bTagAlgo", "particleNet")
#         ]
#         print(f"btagalgo = {bTagAlgo}")
#         bTagWPDict = corrections.btag.getWPValues()
#         print(f"bTagWPDict = {bTagWPDict}")
#         df_singleMu, df_singleAndDiMu, columns = GetAllDfStuff(rdf, columns, bTagAlgo, bTagWPDict)

#         # Uniamo i dizionari di ID e Iso per il loop di plotting
#         all_WPs = {**ID_WPs, **Iso_WPs}

#         for var_name, wp_subdict in all_WPs.items():
#             # Filtro opzionale se vuoi debuggare solo una variabile
#             # if "pfRelIso" not in var_name: continue

#             print(f"Plotting variable: {var_name}")

#             plot_id_iso_vars(df_singleMu, var_name, wp_subdict,
#                              formatted_plot_path, process_name, args.year, "singleMuOnly")

#             plot_id_iso_vars(df_singleAndDiMu, var_name, wp_subdict,
#                              formatted_plot_path, process_name, args.year, "singleMuOrDiMu")

