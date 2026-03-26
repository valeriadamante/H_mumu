from FLAF.Common.Utilities import *

channels = ["muMu"]
leg_names = ["Electron", "Muon", "Tau"]


def LowerMassCut(df, p4_cols=["p4"], cut_value=50):
    for p4_col in p4_cols:
        df = df.Define(f"m_mumu_{p4_col}", f"(mu1_{p4_col}+mu2_{p4_col}).M()")
    masses_cut = " || ".join([f"m_mumu_{p4_col} > {cut_value}" for p4_col in p4_cols])
    df = df.Filter(masses_cut, masses_cut)
    return df


def LeptonsSelection(df):
    ### muon selection: pt > 15 GeV, abs(eta) < 2.4, medium ID, loose PF Iso ###
    df = df.Define(
        "Muon_acceptanceSel",
        "v_ops::pt(Muon_p4) > 15 && abs(v_ops::eta(Muon_p4)) < 2.4",
    )
    df = df.Define(
        "Muon_idIsoSel",
        "Muon_mediumId && Muon_pfIsoId >= 2",
    )
    df = df.Define("Muon_selectedIdx", "Muon_idx[Muon_acceptanceSel && Muon_idIsoSel]")
    df = df.Filter("Muon_selectedIdx.size()==2", "n_muons=2")
    df = df.Define(
        "Muon_selectedIdxSorted",
        """
                    auto indices = Muon_selectedIdx;
                    if(Muon_p4[indices[1]].pt() > Muon_p4[indices[0]].pt())
                        std::swap(indices[0], indices[1]);
                    return indices; """,
    )
    df = df.Define("mu1_idx", "Muon_selectedIdxSorted[0]")
    df = df.Define("mu2_idx", "Muon_selectedIdxSorted[1]")
    # df = df.Filter("Muon_charge[mu1_idx]*Muon_charge[mu2_idx]<0", "OS muons")

    ### electron veto ###
    df = df.Define(
        "Electron_B0_veto",
        "v_ops::pt(Electron_p4) > 20 && abs(v_ops::eta(Electron_p4)) < 2.5  && Electron_mvaIso_WP90",
    )
    # && abs(Electron_dz) < 0.2 && abs(Electron_dxy) < 0.024 --> to add?
    df = df.Filter("Electron_idx[Electron_B0_veto].size() == 0", "No extra electrons")
    return df


def LeptonsSelection_dev(df):
    df = df.Define(
        "Muon_B0", f"""(v_ops::pt(Muon_p4) > 15 && abs(v_ops::eta(Muon_p4)) < 2.4)"""
    )
    big_ID_OR = "(Muon_looseId || Muon_mvaLowPt > -0.6)"
    big_Iso_OR = "(Muon_pfIsoId >= 2 || Muon_miniIsoId >= 1 ||  Muon_miniPFRelIso_all <= 0.4 || Muon_pfRelIso04_all <= 0.25)"
    big_ID_Iso_OR = f"({big_ID_OR} || {big_Iso_OR})"
    df = df.Define("Big_ID_Iso_OR", big_ID_Iso_OR)
    df = df.Filter(
        "Muon_idx[Muon_B0 && Big_ID_Iso_OR].size()==2",
        "Consider events with exactly 2 muons",
    )
    df = df.Define("mu1_idx", "Muon_idx[Muon_B0 && Big_ID_Iso_OR][0]")
    df = df.Define("mu2_idx", "Muon_idx[Muon_B0 && Big_ID_Iso_OR][1]")

    ### electron veto ###
    df = df.Define(
        "Electron_B0_veto",
        f"""v_ops::pt(Electron_p4) > 20 && abs(v_ops::eta(Electron_p4)) < 2.5  && ( Electron_mvaIso_WP90 == true )""",
    )  # && abs(Electron_dz) < 0.2 && abs(Electron_dxy) < 0.024 --> to add?
    df = df.Filter("Electron_idx[Electron_B0_veto].size() == 0", "No extra electrons")
