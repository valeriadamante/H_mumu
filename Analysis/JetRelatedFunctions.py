import ROOT

if __name__ == "__main__":
    sys.path.append(os.environ["ANALYSIS_PATH"])


from FLAF.Common.Utilities import *

JetObservables = [
    "PNetRegPtRawCorr",
    "PNetRegPtRawCorrNeutrino",
    "PNetRegPtRawRes",
    "UParTAK4RegPtRawCorr",
    "UParTAK4RegPtRawCorrNeutrino",
    "UParTAK4RegPtRawRes",
    "area",
    "btagDeepFlavB",
    "btagDeepFlavCvB",
    "btagDeepFlavCvL",
    "btagDeepFlavQG",
    "btagPNetB",
    "btagPNetCvB",
    "btagPNetCvL",
    "btagPNetCvNotB",
    "btagPNetQvG",
    "btagPNetTauVJet",
    "btagUParTAK4B",
    "btagUParTAK4CvB",
    "btagUParTAK4CvL",
    "btagUParTAK4CvNotB",
    "btagUParTAK4QvG",
    "btagUParTAK4TauVJet",
    "chEmEF",
    "chHEF",
    "chMultiplicity",
    "electronIdx1",
    "electronIdx2",
    "eta",
    "genJetIdx",
    "hfEmEF",
    "hfHEF",
    "hfadjacentEtaStripsSize",
    "hfcentralEtaStripSize",
    "hfsigmaEtaEta",
    "hfsigmaPhiPhi",
    "jetId",
    "mass",
    "muEF",
    "muonIdx1",
    "muonIdx2",
    "muonSubtrFactor",
    "nConstituents",
    "nElectrons",
    "nMuons",
    "nSVs",
    "neEmEF",
    "neHEF",
    "neMultiplicity",
    "phi",
    "pt",
    "puIdDisc",
    "puId_beta",
    "puId_dR2Mean",
    "puId_frac01",
    "puId_frac02",
    "puId_frac03",
    "puId_frac04",
    "puId_jetR",
    "puId_jetRchg",
    "puId_majW",
    "puId_minW",
    "puId_nCharged",
    "puId_ptD",
    "puId_pull",
    "rawFactor",
    "svIdx1",
    "svIdx2",
    "ptRes",
    "vetoMap",
    "passJetIdTight",
    "passJetIdTightLepVeto",
    "isInsideVetoRegion",
]

JetObservablesMC = ["hadronFlavour", "partonFlavour", "genJetIdx"]


def JetCollectionDef(df, bTagAlgo, LooseWPValue, MediumWPValue, mu_suff="ScaRe_FSR"):
    if "Jet_idx" not in df.GetColumnNames():
        print("Jet_idx not in df.GetColumnNames")
        df = df.Define(f"Jet_idx", f"CreateIndexes(Jet_pt.size())")
    df = df.Define(
        f"Jet_p4",
        f"GetP4(Jet_pt, Jet_eta, Jet_phi, Jet_mass, Jet_idx)",
    )

    #### Jet PreSelection ####
    df = df.Define(
        "Jet_preSel",
        f"""v_ops::pt(Jet_p4) > 25 && abs(v_ops::eta(Jet_p4))< 4.7 && (Jet_passJetIdTight) """,
    )
    # ed on “loose” selection: pT > 15 GeV and |η|<4.7 and passTightLepVetoId and (chEmEF + neEmEF) < 0.9)
    df = df.Define(
        "Jet_preSel_andDeadZoneVetoMap",
        "Jet_preSel && !Jet_vetoMap",
    )

    df = df.Define(
        f"Jet_NoOverlapWithMuons",
        f"RemoveOverlaps(Jet_p4, Jet_preSel_andDeadZoneVetoMap, {{{{mu1_p4_{mu_suff}, mu2_p4_{mu_suff}}}}}, 2, 0.4)",
    )

    df = df.Define(
        "Jet_IsInsideHornVetoRegion",
        "( abs(v_ops::eta(Jet_p4)) > 2.5 && v_ops::pt(Jet_p4) < 50 ) ",
    )  # questo va x era

    df = df.Define(
        "Jet_IsOutsideOfHornVetoRegion",
        "( abs(v_ops::eta(Jet_p4)) < 2.5 || v_ops::pt(Jet_p4) > 50 ) ",
    )  # questo va x era

    df = df.Define(
        "goodJet",
        "Jet_NoOverlapWithMuons && !Jet_IsInsideHornVetoRegion",
    )  # questo va x era
    # exclude completely the jets in Horn region
    df = df.Define(
        f"SelectedJet_p4",
        f"Jet_p4[goodJet]",
    )
    df = df.Define(
        f"SelectedJet_index",
        f"Jet_idx[goodJet]",
    )

    df = df.Define(f"N_SelectedJets", "SelectedJet_index.size()")

    #### Final state definitions: removing bTagged jets ####

    df = df.Define(
        "Jet_btag_Veto_loose",
        f"Jet_btag{bTagAlgo} >= {LooseWPValue}  && abs(v_ops::eta(Jet_p4))< 2.5 ",
    )
    df = df.Define(
        "Jet_btag_Veto_medium",
        f"Jet_btag{bTagAlgo} >= {MediumWPValue} && abs(v_ops::eta(Jet_p4))< 2.5 ",
    )
    df = df.Define(
        "JetTagSel",
        "Jet_p4[goodJet && Jet_btag_Veto_medium].size() < 1  && Jet_p4[goodJet && Jet_btag_Veto_loose].size() < 2 ",
    )
    return df


def JetObservablesDef(df):
    jet_names = {
        0: "leading",
        1: "subleading",
        2: "third",
        3: "fourth",
    }
    for jet_idx, jet_type in jet_names.items():
        for jet_obs in ["pt", "eta", "phi", "rapidity"]:

            df = df.Define(
                f"{jet_type}jet_{jet_obs}",
                f"if (SelectedJet_index.size()>{jet_idx}) return static_cast<float>(v_ops::{jet_obs}(SelectedJet_p4)[{jet_idx}]); else return -1000.f;",
            )
        df = df.Define(
            f"{jet_type}jet_p4",
            f"if (SelectedJet_index.size()>{jet_idx}) return SelectedJet_p4.at({jet_idx}); else return ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(0.,0.,0.,0.);",
        )
    # define Jet HT:
    if "SelectedJet_pt" not in df.GetColumnNames():
        df = df.Define("SelectedJet_pt", "v_ops::pt(SelectedJet_p4)")
    df = df.Define(
        f"SelectedJets_HT",
        "float SelectedJet_HT; for(size_t jet_idx = 0; jet_idx < SelectedJet_pt.size() ; jet_idx++){SelectedJet_HT+=SelectedJet_pt[jet_idx];} return SelectedJet_HT;",
    )
    # df.Display({"SelectedJets_HT"}).Print()

    df = df.Define(f"delta_eta_jj_ls", "std::abs(leadingjet_eta - subleadingjet_eta)")
    df = df.Define(f"m_jj_ls", "(leadingjet_p4+subleadingjet_p4).M()")

    return df


def VBFNetJetCollectionDef(df, max_jets=4):
    # Define jets for VBF selector network
    df = df.Define(
        "VBFCandJet_selection",
        "Jet_NoOverlapWithMuons && Jet_pt > 20 && ((ROOT::VecOps::abs(Jet_eta) < 2.5) || (ROOT::VecOps::abs(Jet_eta) > 3.0) || (Jet_pt > 50));",
    )
    # Add the desired variables
    jet_vars = [
        "pt",
        "eta",
        "phi",
        # "btagPNetB",
        # "btagPNetCvB",
        "btagPNetCvL",
        # "btagPNetCvNotB",
        "btagPNetQvG",
        # "btagPNetTauVJet",
        # "puIdDisc",
    ]
    for var in jet_vars:
        df = df.Define(f"FilteredJet_{var}_vec", f"Jet_{var}[VBFCandJet_selection]")
        for i in range(max_jets):
            df = df.Define(
                f"FilteredJet_{var}_{i+1}",
                f"static_cast<float>(FilteredJet_{var}_vec.size()>{i} ? FilteredJet_{var}_vec[{i}] : 0.0)",
            )

    return df


def VBFJetSelection(df):
    df = df.Define("VBFJetCand", "FindVBFJets(Jet_p4, goodJet)")
    df = df.Define("HasVBF", "return static_cast<bool>(VBFJetCand.isVBF) ")

    df = df.Define(
        "m_jj",
        "if (HasVBF) return static_cast<float>(VBFJetCand.m_inv); return -1000.f",
    )
    df = df.Define(
        "delta_eta_jj",
        "if (HasVBF) return static_cast<float>(VBFJetCand.eta_separation); return -1000.f",
    )

    df = df.Define(
        "abs_delta_eta_jj",
        "if (HasVBF) return std::abs(delta_eta_jj); return -1000.f",
    )
    df = df.Define(
        "j1_idx",
        "if (HasVBF) return static_cast<int>(VBFJetCand.leg_index[0]); return -1000; ",
    )
    df = df.Define(
        "j2_idx",
        "if (HasVBF) return static_cast<int>(VBFJetCand.leg_index[1]); return -1000; ",
    )
    df = df.Define(
        "j1_pt",
        "if (HasVBF) return static_cast<float>(VBFJetCand.leg_p4[0].Pt()); return -1000.f; ",
    )
    df = df.Define(
        "j2_pt",
        "if (HasVBF) return static_cast<float>(VBFJetCand.leg_p4[1].Pt()); return -1000.f; ",
    )
    df = df.Define(
        "j1_eta",
        "if (HasVBF) return static_cast<float>(VBFJetCand.leg_p4[0].Eta()); return -1000.f; ",
    )
    df = df.Define(
        "j2_eta",
        "if (HasVBF) return static_cast<float>(VBFJetCand.leg_p4[1].Eta()); return -1000.f; ",
    )
    df = df.Define(
        "j1_phi",
        "if (HasVBF) return static_cast<float>(VBFJetCand.leg_p4[0].Phi()); return -1000.f; ",
    )
    df = df.Define(
        "j2_phi",
        "if (HasVBF) return static_cast<float>(VBFJetCand.leg_p4[1].Phi()); return -1000.f; ",
    )
    df = df.Define(
        "j1_y",
        "if (HasVBF) return static_cast<float>(VBFJetCand.leg_p4[0].Rapidity()); return -1000.f; ",
    )
    df = df.Define(
        "j2_y",
        "if (HasVBF) return static_cast<float>(VBFJetCand.leg_p4[1].Rapidity()); return -1000.f; ",
    )
    df = df.Define(
        "delta_phi_jj",
        "if (HasVBF) return static_cast<float>(ROOT::Math::VectorUtil::DeltaPhi( VBFJetCand.leg_p4[0], VBFJetCand.leg_p4[1] ) ); return -1000.f;",
    )

    df = df.Define(f"pt_jj", "(VBFJetCand.leg_p4[0]+VBFJetCand.leg_p4[1]).Pt()")
    df = df.Define(
        "VBFjets_pt",
        f"RVecF void_pt {{}} ; if (HasVBF) return v_ops::pt(VBFJetCand.legs_p4); return void_pt;",
    )
    df = df.Define(
        "VBFjets_eta",
        f"RVecF void_eta {{}} ; if (HasVBF) return v_ops::eta(VBFJetCand.legs_p4); return void_eta;",
    )
    df = df.Define(
        "VBFjets_phi",
        f"RVecF void_phi {{}} ; if (HasVBF) return v_ops::phi(VBFJetCand.legs_p4); return void_phi;",
    )
    df = df.Define(
        "VBFjets_y",
        f"RVecF void_y {{}} ; if (HasVBF) return v_ops::rapidity(VBFJetCand.legs_p4); return void_y;",
    )
    for var in JetObservables:
        if f"Jet_{var}" not in df.GetColumnNames():
            continue
        if f"j1_{var}" not in df.GetColumnNames():
            df = df.Define(
                "j1_" + var,
                f"if (HasVBF && j1_idx >= 0) return static_cast<float>(Jet_{var}[j1_idx]); return -1000.f;",
            )
        if f"j2_{var}" not in df.GetColumnNames():
            df = df.Define(
                "j2_" + var,
                f"if (HasVBF && j2_idx >= 0) return static_cast<float>(Jet_{var}[j2_idx]); return -1000.f;",
            )

    return df


def SoftJetCollectionCleaningInVBF(df, mu_suff="ScaRe_FSR"):
    if "SoftActivityJet_idx" not in df.GetColumnNames():
        print("SoftActivityJet_idx not in df.GetColumnNames")
        df = df.Define(
            f"SoftActivityJet_idx", f"CreateIndexes(SoftActivityJet_pt.size())"
        )
    if f"SoftActivityJet_mass" not in df.GetColumnNames():
        df = df.Define(
            f"SoftActivityJet_mass",
            "RVecF SoftActivityJet_mass(SoftActivityJet_idx.size,0.); return SoftActivityJet_mass;",
        )
    df = df.Define(
        f"SoftActivityJet_p4",
        f"GetP4(SoftActivityJet_pt, SoftActivityJet_eta, SoftActivityJet_phi, SoftActivityJet_mass, SoftActivityJet_idx)",
    )

    df = df.Define(
        f"SoftJetJet_NoOverlapWithMuons",
    )
    df = df.Define(
        f"SoftJetJet_NoOverlapWithMuons",
        f"RemoveOverlaps(Jet_p4, Jet_preSel_andDeadZoneVetoMap, {{{{mu1_p4_{mu_suff}, mu2_p4_{mu_suff}, VBFJetCand.leg_p4[0], VBFJetCand.leg_p4[1]}}}}, 4, 0.4)",
    )
    df = df.Define(
        f"SoftJetCleanedActivity_pt",
        "v_ops::pt(SoftActivityJet_p4[SoftJetJet_NoOverlapWithMuons])",
    )
    df = df.Define(
        f"SoftJetCleanedActivity_eta",
        "v_ops::eta(SoftActivityJet_p4[SoftJetJet_NoOverlapWithMuons])",
    )
    df = df.Define(
        f"SoftJetCleanedActivity_N",
        "SoftActivityJet_p4[SoftJetJet_NoOverlapWithMuons].size()",
    )
    df = df.Define(
        f"SoftJetCleanedActivity_ptSum",
        "float sum=0.; for(size_t sj_idx=0; sj_idx<SoftJetCleanedActivity_pt.size();sj_idx++){{sum+=SoftJetCleanedActivity_pt[sj_idx];}} return sum;",
    )

    df = df.Define(
        "SoftJetJet_NoOverlapWithMuonsAndEtaCleaning",
        "SoftJetJet_NoOverlapWithMuons && SoftJetCleanedActivity_eta < std::max(j1_eta, j2_eta) && SoftJetCleanedActivity_eta > std::min(j1_eta, j2_eta)",
    )


def VBFJetMuonsObservables(df, mu_suff="ScaRe_FSR"):
    df = df.Define(
        "Zeppenfeld_Var",
        f"if (HasVBF) return static_cast<float>((y_mumu_{mu_suff} - 0.5*(j1_y+j2_y))/std::abs(j1_y - j2_y)); return -10000.f;",
    )
    df = df.Define(
        "pT_all_sum",
        f"if(HasVBF) return static_cast<float>((mu1_p4_{mu_suff}+mu2_p4_{mu_suff}+VBFJetCand.leg_p4[0]+VBFJetCand.leg_p4[1]).Pt()); return -10000.f;",
    )
    df = df.Define(
        "pT_single_sum",
        f"if(HasVBF) return static_cast<float>(mu1_p4_{mu_suff}.Pt()+mu2_p4_{mu_suff}.Pt()+VBFJetCand.leg_p4[0].Pt()+VBFJetCand.leg_p4[1].Pt()); return -10000.f;",
    )
    df = df.Define(
        "R_pt",
        f"if(HasVBF) return static_cast<float>(pT_all_sum/pT_single_sum); return -10000.f;",
    )
    df = df.Define(
        "pT_jj_sum",
        f"if(HasVBF) return static_cast<float>((VBFJetCand.leg_p4[0]+VBFJetCand.leg_p4[1]).Pt()); return -10000.f;",
    )
    df = df.Define(
        "pT_jj_diff",
        f"if(HasVBF) return static_cast<float>((VBFJetCand.leg_p4[0]-VBFJetCand.leg_p4[1]).Pt()); return -10000.f;",
    )
    df = df.Define(
        "pT_mumu_sum",
        f"return static_cast<float>((mu1_p4_{mu_suff}+mu2_p4_{mu_suff}).Pt());",
    )
    df = df.Define(
        "pt_centrality",
        f"if(HasVBF) return static_cast<float>(( (pT_mumu_sum-0.5*(pT_jj_sum)) / pT_jj_diff)); return -10000.f;",
    )

    df = df.Define(
        "minDeltaPhi",
        "if(HasVBF) return static_cast<float>(std::min(ROOT::Math::VectorUtil::DeltaPhi( (mu1_p4+mu2_p4), VBFJetCand.leg_p4[0]), ROOT::Math::VectorUtil::DeltaPhi((mu1_p4+mu2_p4), VBFJetCand.leg_p4[1]) ) )  ; return -10000.f;",
    )
    df = df.Define(
        "minDeltaEta",
        f"if(HasVBF) return static_cast<float>(std::min(std::abs(eta_mumu_{mu_suff} - j1_eta),std::abs(eta_mumu_{mu_suff} - j2_eta))) ; return -10000.f;",
    )
    df = df.Define(
        "minDeltaEtaSigned",
        f"if(HasVBF) return static_cast<float>(std::min((eta_mumu_{mu_suff} - j1_eta),(eta_mumu_{mu_suff} - j2_eta))) ; return -10000.f;",
    )

    return df
