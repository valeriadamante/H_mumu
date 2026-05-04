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
    isData = dataset_name == "data"
    period = global_params["era"]
    dfw = analysis.DataFrameBuilderForHistograms(df, global_params, period)
    new_dfw = analysis.PrepareDfForHistograms(dfw, isData)

    full_res_vars = []
    flavor = global_params.get("histTuple_flavor", "default")
    if flavor == "default":
        full_res_vars = global_params.get("histTuple_fullResolution_variables", [])
    else:
        flavor_entry = global_params["histTuple_flavors"][flavor]
        full_res_vars = flavor_entry.get("fullResolution_variables", [])
        variables = flavor_entry.get("variables", [])
        global_params["variables"] = variables
        global_params["histTuple_fullResolution_variables"] = full_res_vars

    for var in full_res_vars:
        new_dfw.colToSave.append(var)

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
            dfw, global_params.get("mu_pt_for_triggerMatchingAndSF", "pt")
        )

        if df_is_central and global_params["compute_unc_histograms"]:
            print("we should define the trigger weights errors too!!:")
        if df_is_central and global_params["compute_unc_histograms"]:
            defineTriggerWeightsErrors(
                dfw,
                global_params.get("mu_pt_for_triggerMatchingAndSF", "pt"),
            )
        if df_is_central:
            central_df_weights_computed = True

    categories = global_params["categories"]
    process_group = global_params["process_group"]
    process_name = global_params["process_name"]
    isCentral = uncName == "Central"
    muID_WP_for_SF = global_params.get("muIDWP", "Loose")
    muIso_WP_for_SF = global_params.get("muIsoWP", "Medium")
    enable_trigger = "trigger" in corrections.to_apply.keys()
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
