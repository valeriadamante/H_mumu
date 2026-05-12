import importlib
from FLAF.Common.Utilities import *
from FLAF.Common.HistHelper import *
from Corrections.Corrections import Corrections
from Corrections.CorrectionsCore import getSystName, central
from Analysis.GetTriggerWeights import defineTriggerWeights, defineTriggerWeightsErrors
from Analysis.MuonRelatedFunctions import *

initialized = False
analysis = None


def Initialize():
    global initialized
    if not initialized:
        headers_dir = os.path.dirname(os.path.abspath(__file__))
        ROOT.gROOT.ProcessLine(f".include {os.environ['ANALYSIS_PATH']}")
        ROOT.gInterpreter.Declare(f'#include "FLAF/include/HistHelper.h"')
        ROOT.gInterpreter.Declare(f'#include "FLAF/include/Utilities.h"')
        ROOT.gROOT.ProcessLine('#include "FLAF/include/AnalysisTools.h"')
        ROOT.gROOT.ProcessLine('#include "FLAF/include/AnalysisMath.h"')
        ROOT.gInterpreter.Declare(
            f'#include "include/Helper.h"'
        )  # not related to FullEvtId definition but needed for analysis specific purpose. At a certain point it will be moved to analysis specific section.
        initialized = True


def analysis_setup(setup):
    global analysis
    analysis_import = setup.global_params["analysis_import"]
    analysis = importlib.import_module(f"{analysis_import}")


def GetDfw(df, setup, dataset_name):
    global_params = setup.global_params
    corrections = Corrections.getGlobal()
    period = global_params["era"]
    kwargset = (
        {}
    )  # here go the customisations for each analysis eventually extrcting stuff from the global params
    kwargset["isData"] = global_params["process_group"] == "data"
    kwargset["wantTriggerSFErrors"] = (
        global_params["compute_rel_weights"]
        and "trigger" in corrections.to_apply.keys()
    )
    kwargset["colToSave"] = []
    dfw = analysis.DataFrameBuilderForHistograms(
        df, global_params, period, corrections, **kwargset, is_not_Cache=True
    )

    new_dfw = analysis.PrepareDFBuilder(dfw)
    further_cuts = global_params.get("further_cuts", {})
    if further_cuts.keys():
        for key in global_params["further_cuts"].keys():
            vars_to_add = global_params["further_cuts"][key][0]
            for var_to_add in vars_to_add:
                if var_to_add not in new_dfw.colToSave:
                    new_dfw.colToSave.append(var_to_add)
    return new_dfw


central_df_weights_computed = False


def DefineWeightForHistograms(
    *,
    dfw,
    isData,
    uncName,
    uncScale,
    unc_cfg_dict,
    hist_cfg_dict,
    global_params,
    final_weight_name,
    df_is_central,
):
    is_central = uncName == central
    global central_df_weights_computed
    corrections = Corrections.getGlobal()
    if not isData and (not central_df_weights_computed or not df_is_central):
        lepton_legs = ["mu1", "mu2"]
        offline_legs = ["mu1", "mu2"]
        triggers_to_use = set()
        channels = global_params["channelSelection"]
        for channel in channels:
            trigger_list = global_params.get("triggers", {}).get(channel, [])
            for trigger in trigger_list:
                if trigger not in corrections.trigger_dict.keys():
                    raise RuntimeError(
                        f"Trigger does not exist in triggers.yaml, {trigger}"
                    )
                triggers_to_use.add(trigger)
        syst_name = getSystName(uncName, uncScale)
        is_central = uncName == central

        dfw.df, all_weights = corrections.getNormalisationCorrections(
            dfw.df,
            lepton_legs=lepton_legs,
            offline_legs=offline_legs,
            trigger_names=triggers_to_use,
            unc_source=uncName,
            unc_scale=uncScale,
            ana_caches=None,
            return_variations=is_central and global_params["compute_unc_histograms"],
            use_genWeight_sign_only=True,
        )
        print("we should define the trigger weights:")
        defineTriggerWeights(
            dfw, global_params["muons"].get("pt_for_TrgSFEvaluation", "pt")
        )

        if df_is_central and global_params["compute_unc_histograms"]:
            print("we should define the trigger weights errors too!!:")
        if df_is_central and global_params["compute_unc_histograms"]:
            defineTriggerWeightsErrors(
                dfw,
                global_params["muons"].get("pt_for_TrgSFEvaluation", "pt"),
            )
        if df_is_central:
            central_df_weights_computed = True

    categories = global_params["categories"]
    process_group = global_params["process_group"]
    process_name = global_params["process_name"]
    isCentral = uncName == "Central"
    muID_WP_for_SF = global_params["muons"].get("ID_WP", "Medium")
    muIso_WP_for_SF = global_params["muons"].get("Iso_WP", "Loose")
    enable_trigger = False # "trigger" in corrections.to_apply.keys()
    enable_ID = "mu" in corrections.to_apply.keys()
    total_weight_expression = (
        analysis.GetWeight(
            "muMu",
            process_name,
            muID_WP_for_SF,
            muIso_WP_for_SF,
            enable_trigger=enable_trigger,
            enable_ID=enable_ID,
        )
        if process_group != "data"
        else "1"
    )  # are we sure?
    # print(f"the total weight expression is {total_weight_expression}")
    weight_name = "final_weight"
    if weight_name not in dfw.df.GetColumnNames():
        dfw.df = dfw.df.Define(weight_name, total_weight_expression)

    if not isCentral:
        if (
            uncName in unc_cfg_dict["norm"].keys()
            and "expression" in unc_cfg_dict["norm"][uncName].keys()
            and process_name
            in unc_cfg_dict["norm"][uncName].get("processes", [process_name])
        ):
            weight_name = unc_cfg_dict["norm"][uncName]["expression"].format(
                scale=uncScale,
                muID_WP_for_SF=muID_WP_for_SF,
                muIso_WP_for_SF=muIso_WP_for_SF,
            )
    # print(f"Defining final weight: {final_weight_name} as {weight_name}")
    dfw.df = dfw.df.Define(final_weight_name, weight_name)
