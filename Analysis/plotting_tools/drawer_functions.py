import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import ROOT
import re
if __name__ == "__main__":
    sys.path.append(os.environ["ANALYSIS_PATH"])
from Analysis.plotting_tools.HelpersForHistograms import findBinEntry, resolve_text_positions, folder_names, period_dict, get_hist_arrays
hep.style.use("CMS")


def is_2d_histogram(axes_cfg_dict, variable):
    var_entry = findBinEntry(axes_cfg_dict, variable)
    if var_entry in axes_cfg_dict and "var_list" in axes_cfg_dict[var_entry]:
        return True
    return False


def get_bin_edges_widths(hist):
    nbins = hist.GetNbinsX()
    bin_edges = np.array([hist.GetBinLowEdge(i) for i in range(1, nbins + 2)])
    bin_widths = np.array([hist.GetBinWidth(i) for i in range(1, nbins + 1)])
    return bin_edges, bin_widths


def integral(hist, divide_by_bin_width=False):
    vals, _, _, bin_widths = get_hist_arrays(hist, divide_by_bin_width)
    return float(np.sum(vals * bin_widths) if divide_by_bin_width else np.sum(vals))


def compute_total_mc_and_stat_err(mc_hists, divide_by_bin_width=False):
    if not mc_hists: return None, None
    first_hist = next(iter(mc_hists.values()))
    nbins = first_hist.GetNbinsX()
    total_vals = np.zeros(nbins, dtype=float)
    total_errs2 = np.zeros(nbins, dtype=float)
    for h in mc_hists.values():
        vals, errs, _, _ = get_hist_arrays(h, divide_by_bin_width)
        total_vals += vals
        total_errs2 += errs ** 2
    return total_vals, np.sqrt(total_errs2)


def choose_reference_binning(histograms_dict):
    for h in histograms_dict.values():
        if h is not None: return h
    return None


def order_mc_contributions(mc_hists, divide_by_bin_width=False):
    names = list(mc_hists.keys())
    in_order = []
    remaining = [n for n in names if n not in in_order]
    return in_order + list(reversed(remaining))


def draw_mc_stack(ax, mc_hists, processes_dict, bin_edges, divide_by_bin_width, page_cfg_dict):
    if not mc_hists: return None, None
    order = order_mc_contributions(mc_hists, divide_by_bin_width)
    mc_vals, mc_labels, mc_colors = [], [], []
    for name in order:
        h = mc_hists[name]
        vals, _, _, _ = get_hist_arrays(h, divide_by_bin_width)
        mc_vals.append(vals)
        cfg = processes_dict.get(name, {})
        mc_labels.append(cfg.get("name", name) + f'[{h.Integral():.2f}]')
        mc_colors.append(cfg.get("color_mplhep", "gray"))

    total_mc_vals, total_mc_errs = compute_total_mc_and_stat_err(mc_hists, divide_by_bin_width)
    hep.histplot(mc_vals, bins=bin_edges, stack=True, histtype="fill", label=mc_labels,
                 facecolor=mc_colors, edgecolor="black", linewidth=0.5, ax=ax)
    hep.histplot(total_mc_vals, bins=bin_edges, histtype="step", color="black", linewidth=0.5, ax=ax)

    bkg_unc_cfg = page_cfg_dict.get('bkg_unc_hist', {})
    unc_hatch = '//' if bkg_unc_cfg.get('fill_style') == 3013 else None
    unc_alpha = bkg_unc_cfg.get('alpha', 0.35)
    y_up, y_dn = total_mc_vals + total_mc_errs, np.maximum(total_mc_vals - total_mc_errs, 0.0)
    ax.fill_between(bin_edges[:-1], y_dn, y_up, step="post", facecolor="none", edgecolor="black",
                    hatch=unc_hatch, alpha=unc_alpha, linewidth=0.8, label=bkg_unc_cfg.get('legend_title', 'Bkg. unc.'))
    return total_mc_vals, total_mc_errs


def draw_signals(ax, signal_hists, processes_dict, bin_edges, divide_by_bin_width, wantSignal):
    if not wantSignal or not signal_hists: return
    for name, h in signal_hists.items():
        cfg = processes_dict.get(name, {})
        scale = cfg.get("scale", 1.0)
        vals, _, _, _ = get_hist_arrays(h, divide_by_bin_width, scale)
        label = cfg.get("name", name)
        if scale != 1.: label += f"x{scale}"
        label += f"[{h.Integral():.2f}]"
        hep.histplot(vals, bins=bin_edges, histtype="step", label=label, color=cfg.get("color_mplhep", "red"),
                     linestyle="--", linewidth=1.5, ax=ax)


def draw_data(ax, data_hist, bin_edges, divide_by_bin_width, wantData=True, blind_region=None):
    if not wantData or data_hist is None: return None, None
    vals, errs, _, _ = get_hist_arrays(data_hist, divide_by_bin_width)
    if blind_region:
        x_min, x_max = blind_region
        mask = (bin_edges[:-1] >= x_min) & (bin_edges[:-1] < x_max)
        vals[mask], errs[mask] = 0.0, 0.0
    label = f"data [{data_hist.Integral():.2f}]"
    hep.histplot(vals, bins=bin_edges, yerr=errs, histtype="errorbar", label=label, color="black", ax=ax)
    return vals, errs


def draw_ratio(ax_ratio, bin_edges, data_vals, data_errs, total_mc_vals, total_mc_errs, x_label, blind_region):
    if data_vals is None or total_mc_vals is None: return
    ratio = np.divide(data_vals, total_mc_vals, out=np.zeros_like(data_vals), where=total_mc_vals != 0)
    ratio_err = np.abs(np.divide(data_errs, total_mc_vals, out=np.zeros_like(data_errs), where=total_mc_vals != 0))
    mc_rel_unc = np.divide(total_mc_errs, total_mc_vals, out=np.zeros_like(total_mc_errs), where=total_mc_vals != 0)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    y_up, y_dn = 1.0 + mc_rel_unc, np.maximum(1.0 - mc_rel_unc, 0.0)
    mask = np.ones_like(ratio, dtype=bool)
    if blind_region and len(blind_region) == 2:
        x_min, x_max = blind_region
        mask = ~((bin_centers >= x_min) & (bin_centers <= x_max))
        ratio[~mask] = y_dn[~mask] = y_up[~mask] = 0.0
    ax_ratio.fill_between(bin_centers, y_dn, y_up, where=y_dn > 0, step="mid", facecolor="ghostwhite",
                           edgecolor="black", hatch='//', alpha=0.5, zorder=1)
    ax_ratio.errorbar(bin_centers[mask], ratio[mask], yerr=ratio_err[mask], fmt='.', color='black', markersize=10, zorder=2)
    ax_ratio.axhline(1.0, color="black", linestyle="--", linewidth=1.0)
    delta = np.abs(ratio[mask] - 1).mean() if len(ratio[mask]) else 0.5
    ax_ratio.set_ylim(round(1 - delta, 2) * 0.9, round(1 + delta, 2) * 1.1)
    ax_ratio.set_ylabel("Data/MC"); ax_ratio.set_xlabel(x_label)


def draw_ratio_comparison(ax_ratio, bin_edges, ref_vals, ref_errs, reg_vals, reg_errs, x_label, y_label, blind_region, color):
    if ref_vals is None or reg_vals is None: return
    ratio = np.divide(reg_vals, ref_vals, out=np.zeros_like(ref_vals), where=ref_vals != 0)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    mask = np.ones_like(ratio, dtype=bool)
    ax_ratio.errorbar(bin_centers[mask], ratio[mask], yerr=None, fmt='.', color=color, markersize=10, zorder=2)
    ax_ratio.axhline(1.0, color="black", linestyle="--", linewidth=1.0)
    delta = np.abs(ratio[mask] - 1).mean() if len(ratio[mask]) else 0.5
    ax_ratio.set_ylim(round(1 - delta, 2) * 0.8, round(1 + delta, 2) * 1.2)
    ax_ratio.set_ylabel(y_label, fontsize=11); ax_ratio.set_xlabel(x_label)


def plot_2d_histogram_from_config(variable, histograms_dict, phys_model_dict, processes_dict,
                                 axes_cfg_dict, page_cfg_dict, page_cfg_custom_dict, filename_base, period,
                                 wantLogX=False, wantLogY=False, wantData=False, wantSignal=False,
                                 category=None, channel=None):

    var_entry = findBinEntry(axes_cfg_dict, variable)
    hist_cfg = axes_cfg_dict.get(var_entry, {})

    canvas_size = page_cfg_dict['page_setup'].get('canvas_size', [1000, 800])
    fig = plt.figure(figsize=(canvas_size[0] / 80, canvas_size[1] / 100))
    ax = fig.add_subplot(1, 1, 1)
    fig.subplots_adjust(top=0.85, right=0.85)

    mc_hists, signal_hists, data_hist = {}, {}, None
    y_maxes = []

    for contrib, hist in histograms_dict.items():
        if hist is None: continue
        if not hist.InheritsFrom("TH2"): continue
        if contrib in phys_model_dict.get('data', []) + ['data']:
            data_hist = hist
        elif contrib in phys_model_dict.get('signals', []):
            signal_hists[contrib] = hist
        elif contrib in phys_model_dict.get('backgrounds', []):
            mc_hists[contrib] = hist

    all_hists = {**mc_hists}
    if wantSignal: all_hists.update(signal_hists)
    if wantData and data_hist: all_hists['data'] = data_hist

    if not all_hists:
        print("[plot_2d] No valid 2D histograms."); return

    n_plots = len(all_hists)
    n_cols = min(3, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols

    fig.set_size_inches(canvas_size[0] / 80 * n_cols, canvas_size[1] / 100 * n_rows)
    gs = fig.add_gridspec(n_rows, n_cols, hspace=0.3, wspace=0.3)

    colors_2d = plt.cm.viridis(np.linspace(0, 1, 256))

    for idx, (name, hist) in enumerate(all_hists.items()):
        row = idx // n_cols
        col = idx % n_cols
        ax = fig.add_subplot(gs[row, col]) if n_plots > 1 else fig.add_subplot(gs[0, 0])

        h = hist
        h.SetOption("COLZ")

        x_bins = h.GetNbinsX()
        y_bins = h.GetNbinsY()
        x_edges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, x_bins + 2)])
        y_edges = np.array([h.GetYaxis().GetBinLowEdge(i) for i in range(1, y_bins + 2)])

        vals = np.array([[h.GetBinContent(ix, iy) for ix in range(1, x_bins + 1)] for iy in range(1, y_bins + 1)])

        im = ax.pcolormesh(x_edges, y_edges, vals, cmap='viridis', shading='flat')
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("Events", fontsize=10)

        cfg = processes_dict.get(name, {})
        plot_label = cfg.get("name", name)
        ax.set_title(plot_label, fontsize=12)

        var_list = hist_cfg.get("var_list", [])
        if len(var_list) >= 2:
            x_entry = findBinEntry(axes_cfg_dict, var_list[0])
            y_entry = findBinEntry(axes_cfg_dict, var_list[1])
            x_label = axes_cfg_dict[x_entry].get("x_title", var_list[0])
            y_label = axes_cfg_dict[y_entry].get("x_title", var_list[1])
        else:
            x_label = hist_cfg.get("x_title", variable)
            y_label = hist_cfg.get("y_title", variable)
        ax.set_xlabel(x_label, fontsize=10)
        ax.set_ylabel(y_label, fontsize=10)

    legend_cfg = page_cfg_dict.get("legend_mplhep", {})
    text_box_names = page_cfg_dict["page_setup"].get("text_boxes_mplhep", [])
    text_box_cfg = {name: page_cfg_dict.get(name, {}) for name in text_box_names}

    try:
        resolved_positions = resolve_text_positions(text_box_cfg)
    except NameError:
        resolved_positions = {}

    for name in text_box_names:
        cfg = text_box_cfg.get(name, {})
        pos = resolved_positions.get(name, cfg.get("pos", [0.02, 1.05]))
        year = period.split('_')[1]
        if year == "all": year = "2022-2023"
        if cfg.get("type") == "cms_mplhep":
            hep.cms.label(label="Preliminary", data=False,
                          ax=ax, loc=0, com=cfg.get("com", "13.6 TeV"), lumi=cfg.get("lumi", period_dict.get(period, "")),
                          year=year, fontsize=cfg.get("text_size", 12))

    plt.savefig(f"{filename_base}.pdf", bbox_inches="tight")
    print(f"2D Plot saved to {filename_base}.pdf")
    plt.savefig(f"{filename_base}.png", bbox_inches="tight")
    print(f"2D Plot saved to {filename_base}.png")
    plt.close()


def plot_histogram_from_config(variable, histograms_dict, phys_model_dict, processes_dict,
                                axes_cfg_dict, page_cfg_dict, page_cfg_custom_dict, filename_base, period,
                                stacked=True, compare_mode=False, compare_vars_mode=False,
                                wantLogX=False, wantLogY=False, wantData=False, wantSignal=False,
                                wantRatio=False, category=None, channel=None,
                                group_minor_contributions=False, minor_fraction=0.001, ref_region=""):

    if is_2d_histogram(axes_cfg_dict, variable):
        return plot_2d_histogram_from_config(
            variable, histograms_dict, phys_model_dict, processes_dict,
            axes_cfg_dict, page_cfg_dict, page_cfg_custom_dict, filename_base, period,
            wantLogX, wantLogY, wantData, wantSignal, category, channel
        )

    var_entry = findBinEntry(axes_cfg_dict, variable)
    hist_cfg = axes_cfg_dict.get(var_entry, {})
    blind_region = hist_cfg.get("blind_region", [])
    divide_by_bin_width = bool(hist_cfg.get("divide_by_bin_width", False))
    plot_opts = dict(page_cfg_dict.get("plot_options", {}))
    enable_kde = False

    canvas_size = page_cfg_dict['page_setup'].get('canvas_size', [1000, 800])
    ratio_plot = bool(wantData and wantRatio and stacked and not compare_mode) or bool(wantRatio and compare_mode)
    fig = plt.figure(figsize=(canvas_size[0] / 80, canvas_size[1] / 100))
    gs = fig.add_gridspec(2 if ratio_plot else 1, 1, height_ratios=[3, 1] if ratio_plot else [2],
                          hspace=0.05 if ratio_plot else 0.25)
    ax = fig.add_subplot(gs[0])
    fig.subplots_adjust(top=0.85)
    ax_ratio = fig.add_subplot(gs[1], sharex=ax) if ratio_plot else None

    mc_hists, signal_hists, data_hist = {}, {}, None
    y_maxes = []

    if compare_mode:
        linestyle_cycle, color_cycle = ['-', '--', ':', '-.'], ['cornflowerblue', 'black', 'red', 'orange', 'gray', 'green', 'cyan', 'blue', 'magenta', 'purple']
        regions = list(histograms_dict.keys())
        ref_region = regions[0] if not ref_region else ref_region
        ref_hist = choose_reference_binning(histograms_dict[ref_region])
        if ref_hist is None: print("[plot] No valid hist for binning."); return
        _, _, bin_edges, _ = get_hist_arrays(ref_hist, False)
        processes = list(set(phys_model_dict.get('backgrounds', []) + phys_model_dict.get('signals')))

        ref_vals_dict, ref_vals_err_dict = {}, {}
        for proc in processes:
            h = histograms_dict[ref_region].get(proc)
            if h is None: continue
            scale = processes_dict.get(proc, {}).get("scale", 1.0)
            vals, val_errs, _, _ = get_hist_arrays(h, divide_by_bin_width, scale)
            ref_vals_dict[proc], ref_vals_err_dict[proc] = vals, val_errs

        for i_region, region in enumerate(regions):
            color = color_cycle[i_region % len(color_cycle)]
            hist_dict = histograms_dict[region]
            for proc in processes:
                if proc not in hist_dict or hist_dict[proc] is None: continue
                h = hist_dict[proc]
                scale = processes_dict.get(proc, {}).get("scale", 1.0)
                plot_vals, plot_errs, _, _ = get_hist_arrays(h, divide_by_bin_width, scale)
                if plot_vals is None: continue
                y_maxes.append(np.max(plot_vals))
                s_label = processes_dict.get(proc, {}).get('name', proc)
                if scale != 1.0: s_label += f"x{scale}"
                plot_label = f"{s_label.split('_')[0]}: {folder_names.get(region, region)}"
                linestyle = linestyle_cycle[i_region % len(linestyle_cycle)]
                hep.histplot(plot_vals, bins=bin_edges, histtype="step", color=color, linestyle=linestyle,
                             linewidth=2, label=plot_label, ax=ax)

                if ratio_plot and proc in ref_vals_dict and region != ref_region:
                    region_unique = '_'.join(list(set(region.split("_")) - set(ref_region.split("_")))) or region
                    ref_region_unique = '_'.join(list(set(ref_region.split("_")) - set(region.split("_")))) or ref_region
                    y_label_ratio = f"$\\frac{{{region_unique}}}{{{ref_region_unique}}}$"
                    x_label = hist_cfg.get("x_title", variable)
                    draw_ratio_comparison(ax_ratio, bin_edges, ref_vals_dict[proc], ref_vals_err_dict[proc],
                                          plot_vals, plot_errs, x_label, y_label_ratio, None, color)

    elif compare_vars_mode:
        linestyle_cycle, color_cycle = ['-', '--', ':', '-.', ':'], ['blue', 'green', 'red', 'cyan']
        var_styles = {var: (linestyle_cycle[i % len(linestyle_cycle)], color_cycle[i % len(color_cycle)])
                      for i, var in enumerate(histograms_dict.keys())}
        first_hist = next(iter(histograms_dict.values()))
        _, _, bin_edges, _ = get_hist_arrays(first_hist, False)
        for var, total_hist in histograms_dict.items():
            if total_hist is None: continue
            values = np.array([total_hist.GetBinContent(i + 1) for i in range(total_hist.GetNbinsX())])
            y_maxes.append(np.max(values))
            style = var_styles.get(var)
            hep.histplot(values, bins=bin_edges, histtype="step", color=style[1], linewidth=len(color_cycle) / 2,
                         label=var, ax=ax, alpha=0.5)

    else:
        for contrib, hist in histograms_dict.items():
            if hist is None: continue
            if contrib in phys_model_dict.get('data', []) + ['data']: data_hist = hist
            elif contrib in phys_model_dict.get('signals', []): signal_hists[contrib] = hist
            elif contrib in phys_model_dict.get('backgrounds', []): mc_hists[contrib] = hist

        ref_hist = choose_reference_binning({**mc_hists, **({'data': data_hist} if data_hist else {}), **signal_hists})
        if ref_hist is None: print("[plot] No valid hist."); return
        _, _, bin_edges, _ = get_hist_arrays(ref_hist, False)

        if group_minor_contributions and mc_hists:
            integrals = {c: mc_hists[c].Integral() for c in mc_hists}
            threshold = minor_fraction * sum(integrals.values())
            minor_contribs = [c for c, val in integrals.items() if val < threshold]
            if minor_contribs:
                objsToMerge = ROOT.TList()
                other_hist = mc_hists[minor_contribs[0]].Clone(f"Other_{ref_hist.GetName()}")
                for minor_contrib in minor_contribs[1:]: objsToMerge.Add(mc_hists[minor_contrib])
                other_hist.Merge(objsToMerge)
                mc_hists = {c: mc_hists[c] for c in mc_hists if c not in minor_contribs}
                mc_hists["Other"] = other_hist

        if stacked:
            total_mc_vals, total_mc_errs = draw_mc_stack(ax, mc_hists, processes_dict, bin_edges, divide_by_bin_width, page_cfg_dict)
        else:
            total_mc_vals, total_mc_errs = compute_total_mc_and_stat_err(mc_hists, divide_by_bin_width)
            for name, h in mc_hists.items():
                vals, _, _, _ = get_hist_arrays(h, divide_by_bin_width)
                cfg = processes_dict.get(name, {})
                hep.histplot(vals, bins=bin_edges, histtype="step", label=cfg.get("title", name),
                             color=cfg.get("color_mplhep", "black"), linewidth=2, ax=ax)
        print(f"total MC vals is {total_mc_vals}")
        for k in mc_hists.keys(): y_maxes.append(mc_hists[k].GetMaximum())
        draw_signals(ax, signal_hists, processes_dict, bin_edges, divide_by_bin_width, wantSignal)
        if data_hist is not None and wantData:
            data_vals, data_errs = draw_data(ax, data_hist, bin_edges, divide_by_bin_width, wantData, blind_region)

    default_y_name = "$\\frac{\\mathrm{Events}}{\\mathrm{bin\\ width}}, \\left( \\frac{1}{\\mathrm{GeV}}\\right)$" if divide_by_bin_width else "Events"
    ax.set_ylabel(hist_cfg.get("y_title", default_y_name), fontsize=14)
    x_label = hist_cfg.get("x_title", variable)
    for mu_idx in [1,2]:
        if f"mu{mu_idx}" in variable: x_label = x_label.format(mu_idx=mu_idx)
    if not ratio_plot: ax.set_xlabel(x_label, fontsize=14)
    else: ax.get_xaxis().set_visible(False)
    ax.set_yscale("log" if wantLogY else "linear")
    ax.set_xscale("log" if wantLogX else "linear")

    y_max = np.max(y_maxes)
    if y_max is not None and np.isfinite(y_max) and y_max > 0:
        max_factor = hist_cfg.get("max_y_sf", 1.2) if not wantLogY else (10 ** hist_cfg.get("max_y_sf", 1.0))
        ax.set_ylim(top=y_max * max_factor)
        if wantLogY: ax.set_ylim(bottom=min(0.00001, y_max * 1e-5))

    ax.set_xlim(bin_edges[0] * 0.99, bin_edges[-1] * 1.01)
    legend_cfg = page_cfg_dict.get("legend_mplhep", {})
    ax.legend(loc='upper right', facecolor=legend_cfg.get("fill_color", "white"), frameon=True,
              fontsize=legend_cfg.get("text_size", 0.02) * 60, framealpha=0.1, ncol=legend_cfg.get("ncols", 2),
              handleheight=1.5, labelspacing=0.1)

    if ratio_plot and not compare_mode and ( (wantData and data_vals is not None) or total_mc_vals is not None):
        ax.set_xlabel(x_label, fontsize=14)
        draw_ratio(ax_ratio, bin_edges, data_vals, data_errs, total_mc_vals, total_mc_errs, x_label=x_label, blind_region=blind_region)

    text_box_names = page_cfg_dict["page_setup"].get("text_boxes_mplhep", [])
    text_box_cfg = {name: page_cfg_dict.get(name, {}) for name in text_box_names}
    try:
        resolved_positions = resolve_text_positions(text_box_cfg)
    except NameError:
        resolved_positions = {}
    if not compare_mode and not compare_vars_mode:
        ax.text(0.4, 0.95, category, transform=ax.transAxes, fontsize=11, verticalalignment='top', horizontalalignment='right')

    for name in text_box_names:
        cfg = text_box_cfg.get(name, {})
        pos = resolved_positions.get(name, cfg.get("pos", [0.02, 1.05]))
        year = period.split('_')[1]
        if year == "all": year = "2022-2023"
        if cfg.get("type") == "cms_mplhep":
            hep.cms.label(label=f"Preliminary", data=("data" in histograms_dict or "Data_Muon" in histograms_dict),
                          ax=ax, loc=0, com=cfg.get("com", "13.6 TeV"), lumi=cfg.get("lumi", period_dict.get(period, "")),
                          year=year, fontsize=cfg.get("text_size", 12))
        else:
            text_content = cfg.get("text", "").format(category=category, channel=channel, variable=variable)
            ax.text(pos[0], pos[1], text_content, transform=ax.transAxes, fontsize=cfg.get("text_size", 10), ha="left", va="top")

    plt.savefig(f"{filename_base}.pdf", bbox_inches="tight")
    print(f"Plot saved to {filename_base}.pdf")
    plt.savefig(f"{filename_base}.png", bbox_inches="tight")
    print(f"Plot saved to {filename_base}.png")
    plt.close()


__all__ = ["plot_histogram_from_config", "plot_2d_histogram_from_config", "is_2d_histogram", "get_bin_edges_widths", "get_hist_arrays", "integral",
           "compute_total_mc_and_stat_err", "choose_reference_binning", "order_mc_contributions",
           "draw_mc_stack", "draw_signals", "draw_data", "draw_ratio", "draw_ratio_comparison"]