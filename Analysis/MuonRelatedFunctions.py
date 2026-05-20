import ROOT

if __name__ == "__main__":
    sys.path.append(os.environ["ANALYSIS_PATH"])


from FLAF.Common.Utilities import *


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
        df = df.Define(
            f"mu1_pt_rel_{newsuff}", f"mu1_p4_{newsuff}.pt()/m_mumu_{newsuff}"
        )
        df = df.Define(
            f"mu2_pt_rel_{newsuff}", f"mu2_p4_{newsuff}.pt()/m_mumu_{newsuff}"
        )

    pt_variants = {
        "_nano": {
            "pt_err_template": "mu{0}_ptErr/mu{0}_pt",
            "pt_name_template": "mu{0}_pt_nano",
        },
        "_nano_scare": {
            "pt_err_template": "mu{0}_ptErr/mu{0}_pt",
            "pt_name_template": "mu{0}_pt_nano_scare",
        },
        "_nano_scare_FSR": {
            "pt_err_template": "mu{0}_ptErr/mu{0}_pt",
            "pt_name_template": "mu{0}_pt_nano_scare_FSR",
        },
        "_bsConstrainedPt": {
            "pt_err_template": "mu{0}_bsConstrainedPtErr/mu{0}_bsConstrainedPt",
            "pt_name_template": "mu{0}_bsConstrainedPt",
        },
        "_bsc_scare": {
            "pt_err_template": "mu{0}_bsConstrainedPtErr/mu{0}_bsConstrainedPt",
            "pt_name_template": "mu{0}_pt_bsc_scare",
        },
        "": {
            "pt_err_template": "mu{0}_bsConstrainedPtErr/mu{0}_bsConstrainedPt",
            "pt_name_template": "mu{0}_pt_bsc_scare",
        },
        "_bsc_scare_FSR": {
            "pt_err_template": "mu{0}_bsConstrainedPtErr/mu{0}_bsConstrainedPt",
            "pt_name_template": "mu{0}_pt_bsc_scare_FSR",
        },
    }  # TO BE FIXED WITH THE SCARE UNC INCLUSION!! --> add friend ttree for scare

    for pt_suffix, pt_info in pt_variants.items():
        # Check if both muons have the required pT columns
        mu1_pt_name = pt_info["pt_name_template"].format(1)
        mu2_pt_name = pt_info["pt_name_template"].format(2)

        if (
            mu1_pt_name not in df.GetColumnNames()
            or mu2_pt_name not in df.GetColumnNames()
        ):
            continue

        # Calculate relative pT errors for each muon
        for mu_idx in [1, 2]:
            sigma_expr = pt_info["pt_err_template"].format(mu_idx)
            df = df.Define(f"sigma_mu{mu_idx}_pt_rel{pt_suffix}", sigma_expr)

        # Calculate m_mumu_resolution: sqrt(0.5*(sigma1^2 + sigma2^2))
        # According to the formula: Δm_μμ^rel = sqrt(1/2 * ((Δpt(u1)/pt(u1))^2 + (Δpt(u2)/pt(u2))^2))
        resolution_expr = f"sqrt(0.5*(pow(sigma_mu1_pt_rel{pt_suffix},2) + pow(sigma_mu2_pt_rel{pt_suffix},2)))"
        resolution_name = (
            f"m_mumu_resolution{pt_suffix}" if pt_suffix != "" else "m_mumu_resolution"
        )
        df = df.Define(resolution_name, resolution_expr)

    return df
