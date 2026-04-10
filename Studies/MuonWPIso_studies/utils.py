
import os
import sys
import ROOT
import glob
import argparse
import pprint

sys.path.append(os.environ["ANALYSIS_PATH"])

import FLAF.Common.Utilities as Utilities
from FLAF.Common.Setup import Setup
from Analysis.MuonRelatedFunctions import GetMuMuP4Observables,GetAllMuMuCorrectedPtRelatedObservables
from Analysis.JetRelatedFunctions import JetCollectionDef
import Analysis.H_mumu as analysis
from Corrections.Corrections import Corrections

samples_dict = {
    "Signal": ["GluGluHto2Mu", "VBFHto2Mu_m125_powheg"],
    "Prompt": ["DY", "EWK"],
    "NonPrompt": ["W_NJets"],
    "TTbar": ["TT"]
}

ID_WPs = {
    "looseId":"mu{idx}_looseId",
    "mediumId":"mu{idx}_mediumId",
    "tightId":"mu{idx}_tightId",
    # "mvaMuID_WP":
    #     {
    #         "tight":"mu{idx}_mvaMuID_WP >= 2",
    #         "medium":"mu{idx}_mvaMuID_WP >= 1",
    #     },
    "mvaLowPt":
        {
        "loose":"mu{idx}_mvaLowPt > -0.6",
        "medium":"mu{idx}_mvaLowPt > -0.2",
        "tight":"mu{idx}_mvaLowPt > 0.15",
        },
}

Iso_WPs = {
    "pfRelIso04_all": {
        "loose":  "mu{idx}_pfRelIso04_all <= 0.25",
        "medium": "mu{idx}_pfRelIso04_all <= 0.20",
        "tight":  "mu{idx}_pfRelIso04_all <= 0.15",
    },
    "pfIsoId": {
        "loose":  "mu{idx}_pfIsoId >= 2",
        "medium": "mu{idx}_pfIsoId >= 3",
        "tight":  "mu{idx}_pfIsoId >= 4",
    },

    # "tkRelIso": {
    #     "loose": "mu{idx}_tkRelIso < 0.10",
    #     "tight": "mu{idx}_tkRelIso < 0.05",
    # },
    # "tkIsoId": {
    #     "loose":  "mu{idx}_tkIsoId >= 1",
    #     "tight":  "mu{idx}_pfIsoId >= 2",
    # },

    "miniPFRelIso_all": {
        "loose":  "mu{idx}_miniPFRelIso_all <= 0.4",
        "medium": "mu{idx}_miniPFRelIso_all <= 0.2",
        "tight":  "mu{idx}_miniPFRelIso_all <= 0.1",
    },

    "miniIsoId": {
        "loose":  "mu{idx}_miniIsoId >= 1",
        "medium": "mu{idx}_miniIsoId >= 2",
        "tight":  "mu{idx}_miniIsoId >= 3",
    },

    # "multiIsoId": {
    #     "loose":  "mu{idx}_multiIsoId >= 1",
    #     "medium": "mu{idx}_multiIsoId >= 2",
    # },
}


# Iso_WPs = {
#     "pfRelIso_loose": "mu{idx}_pfRelIso04_all < 0.25 || mu{idx}_pfIsoId == 2 ", # pfIsoId == 2
#     "pfRelIso_medium": "mu{idx}_pfRelIso04_all < 0.2  || mu{idx}_pfIsoId == 3",
#     "pfRelIso_tight": "mu{idx}_pfRelIso04_all < 0.15  || mu{idx}_pfIsoId == 4",
#     "tkRelIso_loose":"mu{idx}_tkRelIso < 0.1 || mu{idx}_tkIsoId==1",
#     "tkRelIso_tight":"mu{idx}_tkRelIso < 0.05|| mu{idx}_tkIsoId==2",
#     "miniPFRelIso_loose":"mu{idx}_miniPFRelIso_all < 0.4 mu{idx}_miniIsoId==1",
#     "miniPFRelIso_loose":"mu{idx}_miniPFRelIso_all < 0.4 mu{idx}_miniIsoId==1",
#     "miniPFRelIso_medium":"mu{idx}_miniPFRelIso_all < 0.2 mu{idx}_miniIsoId==2",
#     "miniPFRelIso_tight":"mu{idx}_miniPFRelIso_all < 0.1 mu{idx}_miniIsoId==3",
#     "multiIso_loose":"mu{idx}_multiIsoId == 1",
#     "multiIso_medium":"mu{idx}_multiIsoId == 2",
# }

def find_dataset_name(datasets, process_name):
    for process_from_Dataset,process_dict in datasets.items():
        if datasets[process_from_Dataset]['process_name'] == process_name:
            return process_from_Dataset
        if process_name=="DYto2L_InclusivePlusBinned": return "DYto2L_M_50_amcatnloFXFX"
    return process_name
# ---------------------------------------------------------------------
def expand_filelist(patterns):
    files = []
    for p in patterns:
        files.extend(glob.glob(p))
    files = sorted(set(files))
    # if "/eos/home-v/vdamante/H_mumu/anaTuples//v0_FromCentralNanoAOD_AndShifted_NO_massCut_andNOCutsOnIDIso/Run3_2023BPix/DYto2L_M_50_amcatnloFXFX/anaTuple_107.root" in files:
    #     files.remove("/eos/home-v/vdamante/H_mumu/anaTuples//v0_FromCentralNanoAOD_AndShifted_NO_massCut_andNOCutsOnIDIso/Run3_2023BPix/DYto2L_M_50_amcatnloFXFX/anaTuple_107.root")
    # if "/eos/home-v/vdamante/H_mumu/anaTuples//v0_FromCentralNanoAOD_AndShifted_NO_massCut_andNOCutsOnIDIso/Run3_2022/DYto2L_M_50_amcatnloFXFX/anaTuple_13.root" in files:
    #     files.remove("/eos/home-v/vdamante/H_mumu/anaTuples//v0_FromCentralNanoAOD_AndShifted_NO_massCut_andNOCutsOnIDIso/Run3_2022/DYto2L_M_50_amcatnloFXFX/anaTuple_13.root")
    # # print("/eos/home-v/vdamante/H_mumu/anaTuples//v0_FromCentralNanoAOD_AndShifted_NO_massCut_andNOCutsOnIDIso/Run3_2022/DYto2L_M_50_amcatnloFXFX/anaTuple_13.root" in files)

    # if "/eos/home-v/vdamante/H_mumu/anaTuples//v0_FromCentralNanoAOD_AndShifted_NO_massCut_andNOCutsOnIDIso/Run3_2022EE/DYto2L_M_50_amcatnloFXFX/anaTuple_445.root" in files:
    #     files.remove("/eos/home-v/vdamante/H_mumu/anaTuples//v0_FromCentralNanoAOD_AndShifted_NO_massCut_andNOCutsOnIDIso/Run3_2022EE/DYto2L_M_50_amcatnloFXFX/anaTuple_445.root")
    # print("/eos/home-v/vdamante/H_mumu/anaTuples//v0_FromCentralNanoAOD_AndShifted_NO_massCut_andNOCutsOnIDIso/Run3_2022EE/DYto2L_M_50_amcatnloFXFX/anaTuple_445.root" in files)

    return files

#### tutti i path


base_path = '/eos/home-v/vdamante/H_mumu/muon_WPStudies'
miniTuple_path = base_path+"/miniTuples/{period}"
table_path = base_path+"/tables/{period}"
efficiency_path = base_path+"/efficiencies/{period}"
plot_path = base_path+"/plots/{period}"

anaTuple_path = "/eos/user/v/vdamante/H_mumu/testPR/AnaTuples/"
# anaTuple_path = "/eos/home-v/vdamante/H_mumu/NewCustomProd_MuonWPStudies/AnaTuples/"
anaTuples_commonPath = anaTuple_path+"/{period}"

def GetAnaTuplePathDict(period,samples=None,custom_path=None):
    anaPath = anaTuples_commonPath
    if custom_path:
        anaPath = custom_path
    # print(f"considering anaPath {anaPath}")
    total_dict = {
        "VBFHto2Mu_m125_powheg": [
            f"{anaPath.format(period=period)}/VBFHto2Mu_m125_powheg/anaTuple*.root",
            f"{anaPath.format(period=period)}/VBFHto2Mu_M125_powheg/anaTuple*.root",
            f"{anaPath.format(period=period)}/VBFHto2Mu/anaTuple*.root"

        ],
        "VBFHto2Mu_M125_powheg": [
            f"{anaPath.format(period=period)}/VBFHto2Mu_m125_powheg/anaTuple*.root",
            f"{anaPath.format(period=period)}/VBFHto2Mu_M125_powheg/anaTuple*.root",
            f"{anaPath.format(period=period)}/VBFHto2Mu/anaTuple*.root"

        ],
        "GluGluHto2Mu": [
            f"{anaPath.format(period=period)}/GluGluHto2Mu/anaTuple*.root"
        ],
        "DYto2L_InclusivePlusBinned": [
            f"{anaPath.format(period=period)}/DYto2L_M_50_amcatnloFXFX/anaTuple*.root",
            f"{anaPath.format(period=period)}/DYto2L_M_50_0J_amcatnloFXFX/anaTuple*.root",
            f"{anaPath.format(period=period)}/DYto2L_M_50_1J_amcatnloFXFX/anaTuple*.root",
            f"{anaPath.format(period=period)}/DYto2L_M_50_2J_amcatnloFXFX/anaTuple*.root",
            f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_100to200_1J_amcatnloFXFX/anaTuple*.root",
            f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_100to200_2J_amcatnloFXFX/anaTuple*.root",
            f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_200to400_1J_amcatnloFXFX/anaTuple*.root",
            f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_200to400_2J_amcatnloFXFX/anaTuple*.root",
            f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_400to600_1J_amcatnloFXFX/anaTuple*.root",
            f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_400to600_2J_amcatnloFXFX/anaTuple*.root",
            f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_40to100_1J_amcatnloFXFX/anaTuple*.root",
            f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_40to100_2J_amcatnloFXFX/anaTuple*.root",
            f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_600_1J_amcatnloFXFX/anaTuple*.root",
            f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_600_2J_amcatnloFXFX/anaTuple*.root"

            # f"{anaPath.format(period=period)}/DYto2E_M_50_0J_amcatnloFXFX/anaTuple*.root",
            # f"{anaPath.format(period=period)}/DYto2E_M_50_1J_amcatnloFXFX/anaTuple*.root",
            # f"{anaPath.format(period=period)}/DYto2E_M_50_2J_amcatnloFXFX/anaTuple*.root",
            # f"{anaPath.format(period=period)}/DYto2Mu_M_50_0J_amcatnloFXFX/anaTuple*.root",
            # f"{anaPath.format(period=period)}/DYto2Mu_M_50_1J_amcatnloFXFX/anaTuple*.root",
            # f"{anaPath.format(period=period)}/DYto2Mu_M_50_2J_amcatnloFXFX/anaTuple*.root",
            # f"{anaPath.format(period=period)}/DYto2Tau_M_50_0J_amcatnloFXFX/anaTuple*.root",
            # f"{anaPath.format(period=period)}/DYto2Tau_M_50_1J_amcatnloFXFX/anaTuple*.root",
            # f"{anaPath.format(period=period)}/DYto2Tau_M_50_2J_amcatnloFXFX/anaTuple*.root",

        ],
        # "DYto2L": [
        #     f"{anaPath.format(period=period)}/DYto2L_M_50_amcatnloFXFX/anaTuple*.root",
        # ],
        # "DYto2L_jet_binned":[
        #     f"{anaPath.format(period=period)}/DYto2L_M_50_0J_amcatnloFXFX/anaTuple*.root",
        #     f"{anaPath.format(period=period)}/DYto2L_M_50_1J_amcatnloFXFX/anaTuple*.root",
        #     f"{anaPath.format(period=period)}/DYto2L_M_50_2J_amcatnloFXFX/anaTuple*.root"
        # ],
        # "DYto2L_PtBinned": [
        #     f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_100to200_1J_amcatnloFXFX/anaTuple*.root",
        #     f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_100to200_2J_amcatnloFXFX/anaTuple*.root",
        #     f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_200to400_1J_amcatnloFXFX/anaTuple*.root",
        #     f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_200to400_2J_amcatnloFXFX/anaTuple*.root",
        #     f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_400to600_1J_amcatnloFXFX/anaTuple*.root",
        #     f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_400to600_2J_amcatnloFXFX/anaTuple*.root",
        #     f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_40to100_1J_amcatnloFXFX/anaTuple*.root",
        #     f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_40to100_2J_amcatnloFXFX/anaTuple*.root",
        #     f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_600_1J_amcatnloFXFX/anaTuple*.root",
        #     f"{anaPath.format(period=period)}/DYto2L_M_50_PTLL_600_2J_amcatnloFXFX/anaTuple*.root"
        # ],
        "TT": [
            f"{anaPath.format(period=period)}/TTto2L2Nu/anaTuple*.root",
            f"{anaPath.format(period=period)}/TTto4Q/anaTuple*.root",
            f"{anaPath.format(period=period)}/TTtoLNu2Q/anaTuple*.root"
        ],
        "W_NJets": [
            f"{anaPath.format(period=period)}/WtoLNu_0J_amcatnloFXFX/anaTuple*.root",
            f"{anaPath.format(period=period)}/WtoLNu_1J_amcatnloFXFX/anaTuple*.root",
            f"{anaPath.format(period=period)}/WtoLNu_2J_amcatnloFXFX/anaTuple*.root"
        ],
        "EWK":[
            f"{anaPath.format(period=period)}/EWK_2L2J_madgraph_herwig/anaTuple*.root",
        ],
        "DY_all": [
            f"{anaPath.format(period=period)}/DY*/anaTuple*.root",
        ],
        "data":[
            f"{anaPath.format(period=period)}/data/anaTuple*.root",
        ]
    }
    samples_list = samples.split(',') if samples!=None else []
    # print(f"samples are {samples}, and sample list is {samples_list}")
    dict_to_return = total_dict if not samples_list else {}
    for sample_item in samples_list:
        dict_to_return.update(
            {sample_item:total_dict[sample_item]}
        )
    # print(dict_to_return)
    return dict_to_return


#### lumi dict


lumi_dict = {
    "2022":"7.9804",
    "2022EE":"26.6717",
    "2023":"17.794",
    "2023BPix":"9.451",
    "2024":"106.1",
    "2025":"110.73",
    "2022-2023":"61.897100",

}
#### selezioni

mu_pt_for_selection = "pt_Central"
### definitions needed for DFs
trigger_sel = f"( HLT_singleMu && ( (mu1_{mu_pt_for_selection} > 26 && mu1_HasMatching_singleMu) || (mu2_{mu_pt_for_selection} > 26 && mu2_HasMatching_singleMu) ) ) "
trigger_sel_withDiMuon = f"( HLT_singleMu && ( (mu1_{mu_pt_for_selection} > 26 && mu1_HasMatching_singleMu) || (mu2_{mu_pt_for_selection} > 26 && mu2_HasMatching_singleMu) ) ) || ( HLT_diMu && ( (mu1_{mu_pt_for_selection} > 19 && mu1_HasMatching_diMu) && (mu2_{mu_pt_for_selection} > 10 && mu2_HasMatching_diMu) ) )"
SignalRegion_def = "(m_mumu < 135 && m_mumu > 115)"

def DefineAllMuIdIsoVars(df_trgsel, columns,ID_WPs,Iso_WPs):
    for mu_idx in [1,2]:
        for idCut_name,idCut_def in ID_WPs.items():
            if isinstance(idCut_def, str):
                columns.append(f"mu{mu_idx}_{idCut_name}")
                if f"mu{mu_idx}_{idCut_name}" in df_trgsel.GetColumnNames(): continue
                df_trgsel = df_trgsel.Define(f"mu{mu_idx}_{idCut_name}", idCut_def.format(idx=mu_idx))
            elif isinstance(idCut_def, dict):
                for idCut_WP,idCut_WP_def in idCut_def.items():
                    columns.append(f"mu{mu_idx}_{idCut_name}_{idCut_WP}")
                    if f"mu{mu_idx}_{idCut_name}_{idCut_WP}" in df_trgsel.GetColumnNames(): continue
                    df_trgsel = df_trgsel.Define(f"mu{mu_idx}_{idCut_name}_{idCut_WP}", idCut_WP_def.format(idx=mu_idx))
            else:
                raise RuntimeError("cannot define the type of the cut name!!")
        for isoCut_name,isoCut_def in Iso_WPs.items():
            if isinstance(isoCut_def, str):
                if f"mu{mu_idx}_{isoCut_name}" in df_trgsel.GetColumnNames():
                    columns.append(f"mu{mu_idx}_{isoCut_name}")
                    continue
                df_trgsel = df_trgsel.Define(f"mu{mu_idx}_{isoCut_name}", isoCut_def.format(idx=mu_idx))
            elif isinstance(isoCut_def, dict):
                for isoCut_WP,isoCut_WP_def in isoCut_def.items():
                    columns.append(f"mu{mu_idx}_{isoCut_name}_{isoCut_WP}")
                    if f"mu{mu_idx}_{isoCut_name}_{isoCut_WP}" in df_trgsel.GetColumnNames(): continue
                    df_trgsel = df_trgsel.Define(f"mu{mu_idx}_{isoCut_name}_{isoCut_WP}", isoCut_WP_def.format(idx=mu_idx))
    return df_trgsel

def GetAllDfStuff(input_df,columns,bTagAlgo, bTagWPDict):
    print("defining muon p4 observables")
    print(f"before any filtering, yield is {input_df.Sum("weight_base").GetValue()}, corresponding to {input_df.Count().GetValue()} unweighted events")
    input_df = GetMuMuP4Observables(input_df)
    input_df = GetAllMuMuCorrectedPtRelatedObservables(input_df)

    input_df = JetCollectionDef(input_df,
        bTagAlgo,
        bTagWPDict[Utilities.WorkingPointsbTag.Loose],
        bTagWPDict[Utilities.WorkingPointsbTag.Medium])
    # -------------------------
    # Define all WP combinations
    # -------------------------
    wp_columns = []
    columns.extend([
        "weight_base",
    ])
    input_df = input_df.Define("initial_events", f"{input_df.Count().GetValue()}")
    columns.append("initial_events")
    print(f"appending initial events")
    print(f"before filtering, yield is {input_df.Sum("weight_base").GetValue()}, corresponding to {input_df.Count().GetValue()} unweighted events")
    input_df = input_df.Define("initial_weights", f"""{input_df.Sum("weight_base").GetValue()}""")
    columns.append("initial_weights")
    print(f"appending initial weights")

    input_df = input_df.Define("OS", "mu1_charge * mu2_charge < 0").Filter(f"OS && JetTagSel && {SignalRegion_def}")
    print(f"after filtering for OS, jet tag sel and {SignalRegion_def}, there are {input_df.Count().GetValue()} events")
    input_df = input_df.Define("Pass_OS_JetTagSel_SignalRegion_events", f"{input_df.Count().GetValue()}")
    print(f"after filtering for OS && JetTagSel && {SignalRegion_def}, yield is {input_df.Sum("weight_base").GetValue()}, corresponding to {input_df.Count().GetValue()} unweighted events")
    columns.append("Pass_OS_JetTagSel_SignalRegion_events")
    input_df = input_df.Define("Pass_OS_JetTagSel_SignalRegion_weights", f"""{input_df.Sum("weight_base").GetValue()}""")
    columns.append("Pass_OS_JetTagSel_SignalRegion_weights")

    df_singleMuonOnly = input_df.Filter(f"({trigger_sel})")
    print(f"after filtering {trigger_sel}, yield is {df_singleMuonOnly.Sum("weight_base").GetValue()}, corresponding to {df_singleMuonOnly.Count().GetValue()} unweighted events")
    df_singleMuonOnly = df_singleMuonOnly.Define("Pass_TRG_events", f"{df_singleMuonOnly.Count().GetValue()}")
    df_singleMuonOnly = df_singleMuonOnly.Define("Pass_TRG_weights", f"""{df_singleMuonOnly.Sum("weight_base").GetValue()}""")
    print(f"after filtering for singleMu, there are {df_singleMuonOnly.Count().GetValue()} events")
    df_singleAndDiMuonTrg = input_df.Filter(f"({trigger_sel_withDiMuon})")
    print(f"after filtering {trigger_sel}, yield is {df_singleAndDiMuonTrg.Sum("weight_base").GetValue()}, corresponding to {df_singleAndDiMuonTrg.Count().GetValue()} unweighted events")
    df_singleAndDiMuonTrg = df_singleAndDiMuonTrg.Define("Pass_TRG_events", f"{df_singleAndDiMuonTrg.Count().GetValue()}")
    df_singleAndDiMuonTrg = df_singleAndDiMuonTrg.Define("Pass_TRG_weights", f"""{df_singleAndDiMuonTrg.Sum("weight_base").GetValue()}""")
    columns.append("Pass_TRG_events")
    columns.append("Pass_TRG_weights")

    df_singleMuonOnly = DefineAllMuIdIsoVars(df_singleMuonOnly, columns,ID_WPs,Iso_WPs)
    df_singleAndDiMuonTrg = DefineAllMuIdIsoVars(df_singleAndDiMuonTrg, columns,ID_WPs,Iso_WPs)
    return df_singleMuonOnly,df_singleAndDiMuonTrg,list(set(columns))