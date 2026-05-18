from FLAF.Common.Utilities import *

channels = ["muMu"]


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

    # df = df.Filter("Muon_charge[mu1_idx]*Muon_charge[mu2_idx]<0", "OS muons") # this filter can be applied later too.

    ### electron veto ###
    df = df.Define(
        "Electron_B0_veto",
        "v_ops::pt(Electron_p4) > 20 && abs(v_ops::eta(Electron_p4)) < 2.5  && Electron_mvaIso_WP90",
    )
    # && abs(Electron_dz) < 0.2 && abs(Electron_dxy) < 0.024 --> to add?
    df = df.Filter("Electron_idx[Electron_B0_veto].size() == 0", "No extra electrons")
    return df


def DiMuonMassCut(df, p4_cols=["p4"], cut_value=50):
    for p4_col in p4_cols:
        df = df.Define(f"m_mumu_{p4_col}", f"(mu1_{p4_col}+mu2_{p4_col}).M()")
    masses_cut = ""
    if len(p4_cols) > 1:
        masses_cut = " || ".join(
            [f"m_mumu_{p4_col} > {cut_value}" for p4_col in p4_cols]
        )
    elif len(p4_cols) == 1:
        masses_cut = f"m_mumu_{p4_cols[0]} > {cut_value}"
    df = df.Filter(masses_cut, masses_cut)
    return df
