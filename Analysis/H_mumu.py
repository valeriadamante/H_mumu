import ROOT
import sys
import os

if __name__ == "__main__":
    sys.path.append(os.environ["ANALYSIS_PATH"])


from FLAF.Common.Utilities import *
from FLAF.Common.Setup import *
import FLAF.Common.triggerSel as Triggers
from FLAF.Common.HistHelper import *
from Analysis.GetTriggerWeights import *
from Analysis.MuonRelatedFunctions import *
from Analysis.JetRelatedFunctions import *
from Corrections.Corrections import Corrections

for header in [
    "FLAF/include/Utilities.h",
    "include/Helper.h",
    "include/HmumuCore.h",
    "FLAF/include/AnalysisTools.h",
    "FLAF/include/AnalysisMath.h",
]:
    DeclareHeader(os.environ["ANALYSIS_PATH"] + "/" + header)

Taggers_branchesNames = {
    "particleNet": "PNetB",
    "deepJet": "DeepFlavB",
    "UParTAK4": "UParTAK4B",
}


def createKeyFilterDict(global_params, period):
    filter_dict = {}
    filter_str = ""
    channels_to_consider = global_params["channels_to_consider"]
    # sign_regions_to_consider = global_params["MuMuMassRegions"]
    categories = []

    if "categories" in global_params.keys():
        if isinstance(global_params["categories"], dict):
            categories = list(global_params["categories"].keys())
        elif isinstance(global_params["categories"], list):
            categories=global_params['categories']
        elif global_params["categories"]=="all":
            categories = list(global_params["category_definition"].keys())
        else:
            raise RuntimeError(f"global params has a key named categories, but it's neither a list nor a string=all nor a dict, it's: ", global_params["categories"])
    else:
        categories = list(global_params["category_definition"].keys())

    ### add custom categories eventually:
    custom_categories = []
    custom_categories_name = global_params.get(
        "custom_categories", None
    )  # can be extended to list of names
    if custom_categories_name:
        custom_categories = list(global_params.get(custom_categories_name, []))
        if not custom_categories:
            print("No custom categories found")

    ### regions
    custom_regions = []
    custom_regions_name = global_params.get(
        "custom_regions", None
    )  # can be extended to list of names, if for example adding QCD regions + other control regions
    if custom_regions_name:
        custom_regions = list(global_params.get(custom_regions_name, []))
        if not custom_regions:
            print("No custom regions found")

    all_categories = categories + custom_categories
    custom_subcategories = list(global_params.get("custom_subcategories", []))
    triggers_dict = global_params["hist_triggers"]
    for ch in channels_to_consider:
        triggers = triggers_dict[ch]["default"]
        if period in triggers_dict[ch].keys():
            triggers = triggers_dict[ch][period]
        for reg in custom_regions:
            for cat in all_categories:
                filter_base = f" ( {ch} && {triggers} && {reg} && {cat} ) "
                if custom_subcategories:
                    for subcat in custom_subcategories:
                        # filter_base += f"&& {custom_subcat}"
                        filter_str = f"(" + filter_base + f" && {subcat}"
                        filter_str += ")"
                        key = (ch, reg, cat, subcat)
                        filter_dict[key] = filter_str
                else:
                    filter_str = f"(" + filter_base
                    filter_str += ")"
                    key = (ch, reg, cat)
                    filter_dict[key] = filter_str
    return filter_dict


def GetBTagWeight(global_cfg_dict, cat, applyBtag=False):
    btag_weight = "1"
    btagshape_weight = "1"
    global_cfg_dict["ApplyBweight"]
    if global_cfg_dict["ApplyBweight"]:
        if applyBtag:
            if global_cfg_dict["btag_wps"][cat] != "":
                btag_weight = f"weight_bTagSF_{btag_wps[cat]}_Central"
        else:
            if cat not in global_cfg_dict["boosted_categories"] and not cat.startswith(
                "baseline"
            ):
                btagshape_weight = "weight_bTagShape_Central"
    return f"{btag_weight}*{btagshape_weight}"


def SaveVarsForNNInput(variables):
    mumu_vars = [
        "pt_mumu",
        "y_mumu",
        "eta_mumu",
        "phi_mumu",
        "m_mumu",
        "dR_mumu",
        "cosTheta_CS",
        "mu1_pt_rel",
        "mu2_pt_rel",
        "mu1_eta",
        "mu2_eta",
        "phi_CS",
    ]  # , "Ebeam"
    jj_vars = [
        "j1_pt",
        "j1_eta",
        "j2_pt",
        "j2_eta",
        "HasVBF",
        "m_jj",
        "delta_eta_jj",
    ]  # ,"j1_idx","j1_y","j1_phi","delta_phi_jj"
    mumu_jj_vars = [
        "Zepperfield_Var",
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
    ]  # ATTENTION: THESE ARE VECTORS, NOT FLAT OBSERVABLES
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
    ]  # ,"PV_npvs"
    # global_vars = ["FullEventId","luminosityBlock", "run","event", "sample_type", "sample_name", "period", "isData", "nJet"] # ,"PV_npvs"
    for var in global_vars + mumu_vars + jj_vars + mumu_jj_vars + softJets_vars:
        variables.append(var)
    return variables


def GetWeight(
    channel, process_name, muID_WP_for_SF, muIso_WP_for_SF, enable_ID, enable_trigger
):
    weights_to_apply = [
        "weight_base",
    ]
    # quick fix for DY weights. In future should pass the full dataset and process info to DefineWeightForHistograms
    # if process_name.startswith("DY"):
    #     weights_to_apply.extend(
    #         [
    #             "weight_EWKCorr_VptCentral",
    #             "weight_DYw_DYWeightCentral",
    #         ]
    #     )

    trg_weights_dict = {"muMu": ["weight_TrgSF_singleMu_IsoMu24Central"]}

    ID_weights_dict = {
        "muMu": [
            f"weight_mu1_MuonID_SF_{muID_WP_for_SF}ID_TrkCentral",
            f"weight_mu1_MuonID_SF_{muIso_WP_for_SF}PFIso_{muID_WP_for_SF}IDCentral",
            f"weight_mu2_MuonID_SF_{muID_WP_for_SF}ID_TrkCentral",
            f"weight_mu2_MuonID_SF_{muIso_WP_for_SF}PFIso_{muID_WP_for_SF}IDCentral",
        ]
    }

    # should be moved to config
    if enable_ID:
        weights_to_apply.extend(ID_weights_dict[channel])
    if enable_trigger:
        weights_to_apply.extend(trg_weights_dict[channel])

    total_weight = "*".join(weights_to_apply)
    return total_weight


class DataFrameBuilderForHistograms(DataFrameBuilderBase):
    def defineTriggers(self):
        for ch in self.config["channelSelection"]:
            for trg in self.config["triggers"][ch]:
                trg_name = "HLT_" + trg
                self.colToSave.append(trg_name)
                if trg_name not in self.df.GetColumnNames():
                    print(f"{trg_name} not present in colNames")
                    self.df = self.df.Define(trg_name, "1")

    # def defineSampleType(self):
    #     self.df = self.df.Define(
    #         f"sample_type",
    #         f"""std::string process_name = "{self.config["process_name"]}"; return process_name;""",
    #     )

    def defineRegions(self):
        region_defs = self.config["MuMuMassRegions"]
        for reg_name, reg_cut in region_defs.items():
            self.df = self.df.Define(reg_name, reg_cut)
            self.colToSave.append(reg_name)

    def SignRegionDef(self):
        self.df = self.df.Define("OS", "mu1_charge*mu2_charge < 0")
        self.colToSave.append("OS")
        self.df = self.df.Define("SS", "!OS")
        self.colToSave.append("SS")

    def defineCategories(self):  # at the end
        singleMuTh = self.config["singleMu_th"][self.period]
        WP_to_use = self.config["WP_to_use"]
        mu_pt_for_selection = self.config["mu_pt_for_selection"]
        for mu_idx in [1,2]:
            if f"mu{mu_idx}_genPartIdx" not in self.df.GetColumnNames():
                self.df = self.df.Define(f"mu{mu_idx}_genPartIdx", "100")
            if f"mu{mu_idx}_gen_kind" not in self.df.GetColumnNames():
                self.df = self.df.Define(f"mu{mu_idx}_gen_kind", "2")

        for category_to_def in self.config["category_definition"].keys():
            category_name = category_to_def
            cat_str = self.config["category_definition"][category_to_def].format(
                MuPtTh=singleMuTh,
                WP_to_use=WP_to_use,
                mu_pt_for_selection=mu_pt_for_selection,
            )
            self.df = self.df.Define(category_to_def, cat_str)
            self.colToSave.append(category_to_def)

    def defineChannels(self):
        self.df = self.df.Define(f"muMu", f"return true;")
        self.colToSave.append("muMu")

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
        super(DataFrameBuilderForHistograms, self).__init__(df)
        self.config = config
        self.isData = isData
        self.period = period
        self.colToSave = colToSave
        self.wantTriggerSFErrors = wantTriggerSFErrors
        self.corrections = corrections
        self.bTagAlgo = Taggers_branchesNames[
            self.config.get("bTagAlgo", "particleNet")
        ]
        self.bTagWPDict = corrections.btag.getWPValues()


def PrepareDFBuilder(dfBuilder):
    print("Preparing DFBuilder...")
    dfBuilder.df = GetMuMuP4Observables(dfBuilder.df)
    if (
        "muScaRe" in dfBuilder.corrections.to_apply
        and (dfBuilder.config["corrections"]["muScaRe"]["stage"] == "HistTuple" or "m_mumu_resolution" in dfBuilder.config["variables"])
    ):
        dfBuilder.df = dfBuilder.corrections.muScaRe.getP4VariationsForLegs(
            dfBuilder.df
        )

    dfBuilder.df = GetAllMuMuCorrectedPtRelatedObservables(
        dfBuilder.df, suff=dfBuilder.config["mu_pt_for_definitions"]
    ) 
    dfBuilder.defineChannels()
    dfBuilder.defineTriggers()
    dfBuilder.SignRegionDef()

    from FLAF.Common.Utilities import WorkingPointsbTag

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
