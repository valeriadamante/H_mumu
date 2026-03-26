import AnaProd.baseline as AnaBaseline
import FLAF.Common.BaselineSelection as CommonBaseline
from Corrections.Corrections import Corrections
from AnaProd.observables import GetObservablesCols

lepton_legs = ["mu1", "mu2"]
offline_legs = ["mu1", "mu2"]
loadTF = False
loadHHBtag = False


def getDefaultColumnsToSave(isData):
    colToSave = GetObservablesCols("default", isData)
    return list(set(colToSave))


def addAllVariables(
    dfw,
    syst_name,
    isData,
    trigger_class,
    lepton_legs,
    isSignal,
    applyTriggerFilter,
    global_params,
    channels,
    dataset_cfg,
):
    #### baseline cuts (lepton selection + Jet Veto Map definition / application)
    dfw.Apply(AnaBaseline.LeptonsSelection)
    dfw.Apply(Corrections.getGlobal().jet.getEnergyResolution)
    dfw.Apply(Corrections.getGlobal().btag.getWPid, "Jet")
    dfw.Apply(Corrections.getGlobal().JetVetoMap.GetJetVetoMap)
    dfw.Apply(
        CommonBaseline.ApplyJetVetoMap,
        apply_filter=False,
        defineElectronCleaning=True,
        isV12=global_params["nano_version"] == "v12",
    )

    ##### Muons observables ####
    n_legs = 2
    for leg_idx in range(n_legs):

        def LegVar(
            var_name,
            var_expr,
            var_type=None,
            var_cond=None,
            default=0,
        ):
            cond = var_cond
            define_expr = (
                f"static_cast<{var_type}>({var_expr})" if var_type else var_expr
            )
            if cond:
                define_expr = f"{cond} ? ({define_expr}) : {default}"
            dfw.DefineAndAppend(f"mu{leg_idx+1}_{var_name}", define_expr)

        LegVar("legType", "Leg::mu")

        Muon_observables = GetObservablesCols(
            "Muon", isData, global_params["nano_version"]
        )
        for muon_obs in list(set(Muon_observables)):
            muon_obs_split = muon_obs.split("_")
            mu_obs = "_".join(muon_obs_split[1:])
            LegVar(
                mu_obs,
                f"Muon_{mu_obs}.at(mu{leg_idx+1}_idx)",
                var_cond=f"mu{leg_idx+1}_idx>=0",
                default="-100000.f",
            )
        for var in ["pt", "eta", "phi", "mass"]:
            LegVar(
                var,
                f"Muon_p4[mu{leg_idx+1}_idx].{var}()",
                var_cond=f"mu{leg_idx+1}_idx>=0",
                var_type="float",
                default="-1000.f",
            )

        LegVar(
            "pt_nano",
            f"Muon_p4_nano.at(mu{leg_idx+1}_idx).Pt()",
            var_cond=f"mu{leg_idx+1}_idx>=0",
            var_type="float",
            default="-100000.f",
        )
        if not isData:
            LegVar(
                "gen_kind",
                f"genLeptons.at(mu{leg_idx+1}_genMatchIdx).kind()",
                var_type="int",
                var_cond=f"mu{leg_idx+1}_idx>=0 && mu{leg_idx+1}_genMatchIdx>=0",
                default="static_cast<int>(GenLeptonMatch::NoMatch)",
            )
            LegVar(
                "gen_charge",
                f"genLeptons.at(mu{leg_idx+1}_genMatchIdx).charge()",
                var_type="int",
                var_cond=f"mu{leg_idx+1}_idx>=0 && mu{leg_idx+1}_genMatchIdx>=0",
                default="-10",
            )
        else:
            LegVar(
                "gen_kind",
                f"static_cast<int>(GenLeptonMatch::NoMatch)",
                var_cond=f"mu{leg_idx+1}_idx>=0",
                var_type="int",
                default="static_cast<int>(GenLeptonMatch::NoMatch)",
            )
            LegVar(
                "gen_charge",
                f"-10",
                var_cond=f"mu{leg_idx+1}_idx>=0",
                var_type="int",
                default="-10",
            )

        # defining each leg p4 for FindMatching from Muon_p4

        for suffix in ["p4_bsConstrainedPt", "p4_nano", "p4"]:
            if f"mu{leg_idx+1}_{suffix}" not in dfw.df.GetColumnNames():
                dfw.df = dfw.df.Define(
                    f"mu{leg_idx+1}_{suffix}",
                    f"mu{leg_idx+1}_idx >= 0 ? Muon_{suffix}[mu{leg_idx+1}_idx] : LorentzVectorM(0.,0.,0.,0.)",
                )

    dfw.Apply(
        AnaBaseline.LowerMassCut,
        p4_cols=["p4", "p4_nano", "p4_bsConstrainedPt"],
        cut_value=50,
    )

    ##### Jets observables ###
    JetObservables = GetObservablesCols(
        "Jet", isData, global_params["nano_version"]
    ) + ["Jet_idx"]
    dfw.colToSave.extend(list(set(JetObservables)))
    for jvar in ["pt"]:
        jet_obs_name = f"v_ops::{var}(Jet_p4_nano)"
        if f"{jet_obs_name}" in dfw.df.GetColumnNames():
            dfw.DefineAndAppend(f"Jet_{jvar}_nano", jet_obs_name)

    #### LHE weights (special case) ####
    LHE_weights_special = GetObservablesCols(
        "LHEWeight_special", isData, global_params["nano_version"]
    )
    if not isData:
        for var in LHE_weights_special:
            dfw.DefineAndAppend(var[1], var[0])

    #### vars to store: BS, LHE, LHEPart,  ####
    for obs_name in ["BeamSpot", "PV", "SoftActivityJet"]:
        dfw.colToSave.extend(
            GetObservablesCols(obs_name, isData, global_params["nano_version"])
        )
    if not isData:
        for obs_name in ["LHE", "LHEPart"]:
            dfw.colToSave.extend(
                GetObservablesCols(obs_name, isData, global_params["nano_version"])
            )

    if trigger_class is not None:
        hltBranches = dfw.Apply(
            trigger_class.ApplyTriggers,
            lepton_legs,
            isData,
            applyTriggerFilter,
            global_params.get("extraFormat_for_triggerMatchingAndSF", {}),
        )
        dfw.colToSave.extend(hltBranches)
