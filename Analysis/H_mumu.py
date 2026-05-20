import os
import sys
import ROOT

if __name__ == "__main__":
    sys.path.append(os.environ["ANALYSIS_PATH"])

from FLAF.Common.Utilities import *
from FLAF.Common.Setup import *
from FLAF.Common.HistHelper import *

from Analysis.GetTriggerWeights import *
from Analysis.MuonRelatedFunctions import *
from Analysis.JetRelatedFunctions import *

from Corrections.Corrections import Corrections

# =========================
# HEADERS
# =========================
HEADERS = [
    "FLAF/include/Utilities.h",
    "include/Helper.h",
    "include/HmumuCore.h",
    "FLAF/include/AnalysisTools.h",
    "FLAF/include/AnalysisMath.h",
]

for header in HEADERS:
    DeclareHeader(f"{os.environ['ANALYSIS_PATH']}/{header}")


# =========================
# BTAG MAP
# =========================
BTAG_BRANCHES = {
    "particleNet": "PNetB",
    "deepJet": "DeepFlavB",
    "UParTAK4": "UParTAK4B",
}


# =========================
# CHANNELS / REGIONS / CATEGORIES
# =========================
def getChannelsRegionsAndCategoriesNames(global_params):

    channels = global_params["channels"]["selection"]

    categories = list(global_params.get("categories", {}).keys())

    custom_categories = global_params.get("custom_categories", [])
    if isinstance(custom_categories, str):
        custom_categories = global_params.get(custom_categories, [])

    categories = categories + custom_categories
    categories_to_select = global_params.get("categories_to_select", [])
    if not categories_to_select:
        categories_to_select = categories

    custom_subcategories = global_params.get("custom_subcategories", [])

    regions_block = global_params.get("regions", {})
    regions_name = list(regions_block.keys())[0] if regions_block else None
    regions = regions_block.get(regions_name, {}) if regions_name else {}

    return {
        "channels": channels,
        "categories": categories,
        "categories_to_select": categories_to_select,
        "subcategories": custom_subcategories,
        "regions": regions,
    }


# =========================
# KEY FILTER BUILDER
# =========================
def createKeyFilterDict(global_params, period):

    filter_dict = {}

    reg_cat = getChannelsRegionsAndCategoriesNames(global_params)

    channels = reg_cat["channels"]
    regions = reg_cat["regions"]
    categories = reg_cat["categories"]
    subcategories = reg_cat["subcategories"]
    categories_to_select = reg_cat["categories_to_select"]

    triggers_dict = global_params["hist_triggers"]

    for ch in channels:

        trig_cfg = triggers_dict[ch]
        triggers = trig_cfg.get(period, trig_cfg["default"])

        for reg_name, reg_cut in regions.items():

            for cat in categories_to_select:

                base = f"( {ch} && {triggers} && {reg_name} && {cat} )"

                if subcategories:
                    for sub in subcategories:
                        key = (ch, reg_name, cat, sub)
                        filter_dict[key] = f"( {base} && {sub} )"
                else:
                    key = (ch, reg_name, cat)
                    filter_dict[key] = base

    return filter_dict


# =========================
# NN INPUT VARS
# =========================
def SaveVarsForNNInput(variables):

    mumu_vars = [
        "pt_mumu_ScaRe_FSR",
        "y_mumu_ScaRe_FSR",
        "eta_mumu_ScaRe_FSR",
        "phi_mumu_ScaRe_FSR",
        "m_mumu_ScaRe_FSR",
        "dR_mumu_ScaRe_FSR",
        "cosTheta_CS_ScaRe_FSR",
        "mu1_pt_rel_ScaRe_FSR",
        "mu2_pt_rel_ScaRe_FSR",
        "mu1_eta",
        "mu2_eta",
        "phi_CS_ScaRe_FSR",
    ]

    jj_vars = [
        "j1_pt",
        "j1_eta",
        "j2_pt",
        "j2_eta",
        "HasVBF",
        "m_jj",
        "delta_eta_jj",
    ]

    mumu_jj_vars = [
        "Zeppenfeld_Var",
        "R_pt",
        "pt_centrality",
        "minDeltaPhi",
        "minDeltaEta",
        "minDeltaEtaSigned",
    ]  # , "pT_all_sum","pT_jj_sum",
    softJets_vars = [
        "N_softJet",
        "SoftJet_energy",
        "SoftJet_Et",
        "SoftJet_HtCh_fraction",
        "SoftJet_HtNe_fraction",
        "SoftJet_HtHF_fraction",
    ]

    global_vars = [
        "entryIndex",
        "luminosityBlock",
        "run",
        "event",
        "sample_type",
        "sample_name",
        "period",
        "isData",
        "nJet",
    ]

    for v in global_vars + mumu_vars + jj_vars + mumu_jj_vars + softJets_vars:
        variables.append(v)

    return variables


# =========================
# WEIGHTS
# =========================
def GetWeight(
    channel, process_name, muID_WP_for_SF, muIso_WP_for_SF, enable_ID, enable_trigger
):

    weights = ["weight_base"]

    trg = {"muMu": ["weight_TrgSF_singleMu_IsoMu24Central"]}

    ID = {
        "muMu": [
            f"weight_mu1_MuonID_SF_{muID_WP_for_SF}ID_TrkCentral",
            f"weight_mu1_MuonID_SF_{muIso_WP_for_SF}PFIso_{muID_WP_for_SF}IDCentral",
            f"weight_mu2_MuonID_SF_{muID_WP_for_SF}ID_TrkCentral",
            f"weight_mu2_MuonID_SF_{muIso_WP_for_SF}PFIso_{muID_WP_for_SF}IDCentral",
        ]
    }

    if enable_ID:
        weights += ID[channel]
    if enable_trigger:
        weights += trg[channel]

    return "*".join(weights)


# =========================
# DATAFRAME BUILDER
# =========================
class DataFrameBuilderForHistograms(DataFrameBuilderBase):

    def defineTriggers(self):

        for ch in self.config["channels"]["selection"]:
            for trg in self.config["triggers"][ch]:

                name = "HLT_" + trg
                self.colToSave.append(name)

                if name not in self.df.GetColumnNames():
                    self.df = self.df.Define(name, "1")

    def defineRegions(self):

        regions_block = self.config.get("regions", {})
        if not regions_block:
            return

        region_defs = list(regions_block.values())[0]

        for reg_name, reg_cut in region_defs.items():
            self.df = self.df.Define(reg_name, reg_cut)
            self.colToSave.append(reg_name)
        bigOR_string = " || ".join(reg for reg in self.ch_reg_cat_dict["regions"])
        self.df = self.df.Filter(bigOR_string)

    def SignRegionDef(self):

        self.df = self.df.Define("OS", "mu1_charge*mu2_charge < 0")
        self.df = self.df.Define("SS", "!OS")

        self.colToSave += ["OS", "SS"]

    def defineCategories(self):

        singleMuTh = self.config["singleMu_th"][self.period]
        WP = self.config["muons"]["ID_WP"]
        pt_sel = self.config["muons"]["pt_for_selection"]

        for name, expr in self.config["categories"].items():

            expr = expr.format(
                MuPtTh=singleMuTh,
                WP_to_use=WP,
                mu_pt_for_selection=pt_sel,
            )

            self.df = self.df.Define(name, expr)
            self.colToSave.append(name)
        bigOR_string = " || ".join(
            cat for cat in self.ch_reg_cat_dict["categories_to_select"]
        )
        self.df = self.df.Filter(bigOR_string)

    def defineChannels(self):
        for ch in self.config["channels"]["selection"]:
            self.df = self.df.Define(ch, self.config["channels"]["definition"][ch])
            self.colToSave.append(ch)
        bigOR_string = " || ".join(ch for ch in self.ch_reg_cat_dict["channels"])
        self.df = self.df.Filter(bigOR_string)

    def __init__(
        self,
        df,
        config,
        period,
        corrections,
        isData=False,
        wantTriggerSFErrors=False,
        colToSave=[],
        is_not_Cache=False,
    ):

        super().__init__(df)

        self.config = config
        self.isData = isData
        self.period = period
        self.colToSave = colToSave

        self.wantTriggerSFErrors = wantTriggerSFErrors
        self.corrections = corrections

        self.bTagAlgo = BTAG_BRANCHES[
            self.config.get("btag", {}).get("tagger", "particleNet")
        ]
        self.bTagWPDict = corrections.btag.getWPValues()
        self.ch_reg_cat_dict = getChannelsRegionsAndCategoriesNames(self.config)


def PrepareDFBuilder(dfBuilder):
    dfBuilder.df = GetMuonP4Observables(dfBuilder.df)
    if "muScaRe" in dfBuilder.corrections.to_apply:
        dfBuilder.df = dfBuilder.corrections.muScaRe.getP4VariationsForLegs(
            dfBuilder.df
        )
    dfBuilder.df = GetAllMuonsObservablesNew(dfBuilder.df)

    dfBuilder.defineChannels()
    dfBuilder.defineTriggers()
    dfBuilder.SignRegionDef()

    dfBuilder.df = JetCollectionDef(
        dfBuilder.df,
        dfBuilder.bTagAlgo,
        dfBuilder.bTagWPDict[WorkingPointsbTag.Loose],
        dfBuilder.bTagWPDict[WorkingPointsbTag.Medium],
    )
    dfBuilder.df = JetObservablesDef(dfBuilder.df)
    dfBuilder.df = VBFJetSelection(dfBuilder.df)
    dfBuilder.df = VBFJetMuonsObservables(dfBuilder.df)

    dfBuilder.defineRegions()
    dfBuilder.defineCategories()

    return dfBuilder
