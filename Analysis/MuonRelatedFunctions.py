import ROOT

if __name__ == "__main__":
    sys.path.append(os.environ["ANALYSIS_PATH"])


from FLAF.Common.Utilities import *


def GetMuMuMassResolution(df, pt_to_use):
    sigma_pt = {
        "scare": "mu{0}_ptErr/mu{0}_pt",
        "nano": "mu{0}_ptErr/mu{0}_pt",
        "scare_reapplied": "mu{0}_ptErr/mu{0}_pt",
        "BS": "mu{0}_bsConstrainedPtErr/mu{0}_bsConstrainedPt",
        "BS_scare": "mu{0}_bsConstrainedPtErr/mu{0}_bsConstrainedPt",
        "RoccoR": "mu{0}_ptErr/mu{0}_pt",
        "BS_RoccoR": "mu{0}_bsConstrainedPtErr/mu{0}_bsConstrainedPt",
    }
    sigma_scaleandresol = {
        "scare": "0.5*(((mu{0}_pt_1_corr_up-mu{0}_pt_1_scale_corr)*(mu{0}_pt_1_corr_up-mu{0}_pt_1_scale_corr))+((mu{0}_pt_1_corr_dn-mu{0}_pt_1_scale_corr)*(mu{0}_pt_1_corr_dn-mu{0}_pt_1_scale_corr)))",
        "nano": "0.",
        "scare_reapplied": "0.5*(((mu{0}_pt_1_corr_up-mu{0}_pt_1_scale_corr)*(mu{0}_pt_1_corr_up-mu{0}_pt_1_scale_corr))+((mu{0}_pt_1_corr_dn-mu{0}_pt_1_scale_corr)*(mu{0}_pt_1_corr_dn-mu{0}_pt_1_scale_corr)))",
        "BS": "0.",
        "BS_scare": "0.5*(((mu{0}_BS_pt_1_corr_up-mu{0}_BS_pt_1_corr)*(mu{0}_BS_pt_1_corr_up-mu{0}_BS_pt_1_corr))+((mu{0}_BS_pt_1_corr_dn-mu{0}_BS_pt_1_corr)*(mu{0}_BS_pt_1_corr_dn-mu{0}_BS_pt_1_corr)))",
        "RoccoR": "0.",
        "BS_RoccoR": "0.",
    }
    for mu_idx in [1, 2]:
        # print(sigma_pt[pt_to_use].format(mu_idx))
        # print(sigma_scaleandresol[pt_to_use].format(mu_idx))
        df = df.Define(f"sigma_mu{mu_idx}_pt_rel", sigma_pt[pt_to_use].format(mu_idx))
        # df.Display({f"sigma_mu{mu_idx}_pt_rel"}).Print()
        # df=df.Define(f"sigma_mu{mu_idx}_pt_rel", f"sigma_mu{mu_idx}_pt_rel/mu{mu_idx}_pt") # is it alreadt relative??
        df = df.Define(
            f"sigma_mu{mu_idx}_scaleresolution",
            sigma_scaleandresol[pt_to_use].format(mu_idx),
        )
        df = df.Define(
            f"sigma_mu{mu_idx}_scaleresolution_rel",
            f"sigma_mu{mu_idx}_scaleresolution/mu{mu_idx}_pt",
        )
        # df.Display({f"sigma_mu{mu_idx}_scaleresolution"}).Print()
        df = df.Define(
            f"sigma_mu{mu_idx}_total_pt_rel",
            f"sqrt( pow(sigma_mu{mu_idx}_pt_rel,2) + sigma_mu{mu_idx}_scaleresolution_rel )",
        )
        # df.Display({f"sigma_mu{mu_idx}_total_pt_rel"}).Print()
    delta_mu_expr = "0.5*sqrt({0}*{0}+{1}*{1}) "
    # delta_mu_expr = "sqrt( 0.5 * (pow( ({0}/{1}), 2) + pow( ({2}/{3}), 2) ) ) "
    df = df.Define(
        "m_mumu_resolution",
        delta_mu_expr.format("sigma_mu1_total_pt_rel", "sigma_mu2_total_pt_rel"),
    )
    # df.Display({"m_mumu_resolution"}).Print()
    # df = df.Define(
    #     "m_mumu_resolution_nano",
    #     delta_mu_expr.format(
    #         "mu1_ptErr",
    #         "mu1_pt_nano",
    #         "mu2_ptErr",
    #         "mu2_pt_nano",
    #     ),
    # )
    # df = df.Define(
    #     "m_mumu_resolution",
    #     delta_mu_expr.format(
    #         "(mu1_pt-mu1_pt_nano)/mu1_pt",
    #         "mu1_pt",
    #         "(mu2_pt-mu2_pt_nano)/mu2_pt",
    #         "mu2_pt",
    #     ),
    # )

    # df = df.Define(
    #     "m_mumu_resolution_BS",
    #     delta_mu_expr.format(
    #         "mu1_bsConstrainedPtErr",
    #         "mu1_bsConstrainedPt",
    #         "mu2_bsConstrainedPtErr",
    #         "mu2_bsConstrainedPt",
    #     ),
    # )
    # df = df.Define(
    #     "m_mumu_resolution_BS_ScaRe",
    #     delta_mu_expr.format(
    #         "(mu1_BS_pt_1_corr-mu1_bsConstrainedPt)",
    #         "mu1_BS_pt_1_corr",
    #         "(mu2_BS_pt_1_corr-mu2_bsConstrainedPt)",
    #         "mu2_BS_pt_1_corr",
    #     ),
    # )

    return df


def GetMuonP4Observables(df):
    for pt_suffix in [
        "",
        "_bsc_scare",
        "_nano_scare",
        "_nano",
        "_bsConstrainedPt",
    ]:
        for mu_idx in [1, 2]:
            mu_pt_name = (
                f"mu{mu_idx}_pt{pt_suffix}"
                if pt_suffix != "_bsConstrainedPt"
                else f"mu{mu_idx}{pt_suffix}"
            )
            if f"mu{mu_idx}_p4{pt_suffix}" in df.GetColumnNames():
                continue
            if mu_pt_name not in df.GetColumnNames():
                continue
            df = df.Define(
                f"mu{mu_idx}_p4{pt_suffix}",
                f"ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>({mu_pt_name},mu{mu_idx}_eta,mu{mu_idx}_phi,mu{mu_idx}_mass)",
            )
    return df


def GetAllMuonsObservablesNew(df):
    df = df.Define("Ebeam", "13600.0/2")

    dimu_obs = {
        "pt_mumu": "{dimu}.Pt()",
        "m_mumu": "{dimu}.M()",
        "y_mumu": "{dimu}.Rapidity()",
        "eta_mumu": "{dimu}.Eta()",
        "phi_mumu": "{dimu}.Phi()",
        "dR_mumu": "ROOT::Math::VectorUtil::DeltaR({mu1p4}, {mu2p4})",
        "cosTheta_Phi_CS": "ComputeCosThetaPhiCS({mu1p4}, {mu2p4}, Ebeam)",
        "cosTheta_CS": "static_cast<float>(std::get<0>(cosTheta_Phi_CS{suff}))",
        "phi_CS": "static_cast<float>(std::get<1>(cosTheta_Phi_CS{suff}))",
    }
    for pt_suffix in [
        "_nano",
        "_bsConstrainedPt",
        "",  # should be same than bsc_scare
        "_bsc_scare",
        "_nano_scare",
        "_FSR_nano_scare",
        "_FSR_bsc_scare",
    ]:
        for mu_idx in [1, 2]:
            mu_pt_name = (
                f"mu{mu_idx}_pt{pt_suffix}"
                if pt_suffix != "_bsConstrainedPt"
                else f"mu{mu_idx}{pt_suffix}"
            )
            if (
                mu_pt_name in df.GetColumnNames()
                and f"mu{mu_idx}_p4{pt_suffix}" not in df.GetColumnNames()
            ):
                df = df.Define(
                    f"mu{mu_idx}_p4{pt_suffix}",
                    f"ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>({mu_pt_name},mu{mu_idx}_eta,mu{mu_idx}_phi,mu{mu_idx}_mass)",
                )
        p4_dimu = f"(mu1_p4{pt_suffix}+mu2_p4{pt_suffix})"
        p4_dimu_list = [f"mu1_p4{pt_suffix}", f"mu2_p4{pt_suffix}"]
        for obs, expr in dimu_obs.items():
            if pt_suffix == "":
                continue
            df = df.Define(
                f"{obs}{pt_suffix}",
                expr.format(
                    dimu=p4_dimu,
                    mu1p4=p4_dimu_list[0],
                    mu2p4=p4_dimu_list[1],
                    suff=pt_suffix,
                ),
            )
    for mu_idx in [1, 2]:
        df = df.Define(
            f"mu{mu_idx}_p4_noCorr",
            f"mu{mu_idx}_bsConstrainedChi2 < 30 ? mu{mu_idx}_p4_bsConstrainedPt : mu{mu_idx}_p4_nano",
        )
        df = df.Define(
            f"mu{mu_idx}_p4_ScaRe",
            f"mu{mu_idx}_bsConstrainedChi2 < 30 ? mu{mu_idx}_p4_bsc_scare : mu{mu_idx}_p4_nano_scare",
        )
        df = df.Define(
            f"mu{mu_idx}_p4_ScaRe_FSR",
            f"mu{mu_idx}_bsConstrainedChi2 < 30 ? mu{mu_idx}_p4_FSR_bsc_scare : mu{mu_idx}_p4_FSR_nano_scare",
        )
    for newsuff in ["noCorr", "ScaRe", "ScaRe_FSR"]:
        df = df.Define(f"mu1_pt_{newsuff}", f"mu1_p4_{newsuff}.pt()")
        df = df.Define(f"mu2_pt_{newsuff}", f"mu2_p4_{newsuff}.pt()")
        p4_dimu_system = f"(mu1_p4_{newsuff}+mu2_p4_{newsuff})"
        p4_dimu_system_list = [f"mu1_p4_{newsuff}", f"mu2_p4_{newsuff}"]
        for obs, expr in dimu_obs.items():
            df = df.Define(
                f"{obs}_{newsuff}",
                expr.format(
                    dimu=p4_dimu_system,
                    mu1p4=p4_dimu_system_list[0],
                    mu2p4=p4_dimu_system_list[1],
                    suff=f"_{newsuff}",
                ),
            )
    return df
