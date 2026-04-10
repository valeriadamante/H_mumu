import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import ROOT
import matplotlib.ticker as ticker
import yaml
import re
import matplotlib.colors as mcolors
from HelpersForHistograms import *  # keep your helpers (resolve_text_positions etc.)
hep.style.use("CMS")

# -------------------------
# Utilities for hist arrays
# -------------------------
def get_bin_edges_widths(hist):
    nbins = hist.GetNbinsX()
    bin_edges = np.array([hist.GetBinLowEdge(i) for i in range(1, nbins + 2)])
    bin_widths = np.array([hist.GetBinWidth(i) for i in range(1, nbins + 1)])
    return bin_edges, bin_widths


def get_hist_arrays(hist, divide_by_bin_width=False, scale=1.0):
    """
    Return: vals, errs, bin_edges, bin_widths
    vals, errs are numpy arrays length = nbins
    """
    bin_edges, bin_widths = get_bin_edges_widths(hist)
    nbins = hist.GetNbinsX()
    vals = np.array([hist.GetBinContent(i + 1) for i in range(nbins)], dtype=float) * scale
    errs = np.array([hist.GetBinError(i + 1) for i in range(nbins)], dtype=float) * scale
    if divide_by_bin_width:
        vals = np.divide(vals, bin_widths, out=np.zeros_like(vals), where=bin_widths != 0)
        errs = np.divide(errs, bin_widths, out=np.zeros_like(errs), where=bin_widths != 0)
    return vals, errs, bin_edges, bin_widths


def integral(hist, divide_by_bin_width=False):
    vals, _, _, bin_widths = get_hist_arrays(hist, divide_by_bin_width)
    if divide_by_bin_width:
        return float(np.sum(vals * bin_widths))
    return float(np.sum(vals))


def compute_total_mc_and_stat_err(mc_hists, divide_by_bin_width=False):
    """
    Sum MC per-bin and combine statistical errors in quadrature.
    mc_hists: dict name->TH1
    """
    if not mc_hists:
        return None, None
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
    """
    Return the first non-None histogram from a dict-like structure.
    """
    for name, h in histograms_dict.items():
        if h is None:
            continue
        return h
    return None


# -------------------------
# stack ordering
# -------------------------
def order_mc_contributions(mc_hists, divide_by_bin_width=False):
    """
    Default heuristic: reverse insertion order (same as before).
    Could be extended to read config stack_order.
    """
    names = list(mc_hists.keys())
    in_order = []
    remaining = [n for n in names if n not in in_order]
    remaining_reversed = list(reversed(remaining))
    return in_order + remaining_reversed

# -----------------------------
# Draw helpers: MC stack, signals, data
# -----------------------------
def draw_mc_stack(ax, mc_hists, processes_dict, bin_edges, divide_by_bin_width, page_cfg_dict):
    """
    Draw stacked MC and return (total_vals, total_errs)
    """
    if not mc_hists:
        return None, None
    order = order_mc_contributions(mc_hists, divide_by_bin_width)
    mc_vals, mc_labels, mc_colors = [], [], []

    for name in order:
        h = mc_hists[name]
        yield_h = h.Integral()
        vals, _, _, _ = get_hist_arrays(h, divide_by_bin_width)
        mc_vals.append(vals)
        cfg = processes_dict.get(name, {})
        mc_labels.append(cfg.get("name", name)+f'[{yield_h:.2f}]')
        mc_colors.append(cfg.get("color_mplhep", "gray"))

    total_mc_vals, total_mc_errs = compute_total_mc_and_stat_err(mc_hists, divide_by_bin_width)

    hep.histplot(
        mc_vals, bins=bin_edges, stack=True, histtype="fill",
        label=mc_labels, facecolor=mc_colors, edgecolor="black", linewidth=0.5, ax=ax
    )

    hep.histplot(
        total_mc_vals, bins=bin_edges, histtype="step",
        color="black", linewidth=0.5, ax=ax
    )

    bkg_unc_cfg = page_cfg_dict.get('bkg_unc_hist', {})
    unc_hatch = '//' if bkg_unc_cfg.get('fill_style') == 3013 else None
    unc_alpha = bkg_unc_cfg.get('alpha', 0.35)

    y_up = total_mc_vals + total_mc_errs
    y_dn = total_mc_vals - total_mc_errs
    y_dn = np.maximum(y_dn, 0.0)

    ax.fill_between(
        bin_edges[:-1], y_dn, y_up, step="post",
        facecolor="none", edgecolor="black", hatch=unc_hatch, alpha=unc_alpha,
        linewidth=0.8, label=bkg_unc_cfg.get('legend_title', 'Bkg. unc.')
    )

    return total_mc_vals, total_mc_errs


def draw_signals(ax, signal_hists, processes_dict, bin_edges, divide_by_bin_width, wantSignal):
    if not wantSignal or not signal_hists:
        return
    for name, h in signal_hists.items():
        cfg = processes_dict.get(name, {})
        scale = cfg.get("scale", 1.0)
        vals, _, _, _ = get_hist_arrays(h, divide_by_bin_width, scale)
        label = cfg.get("name", name)
        if scale != 1.:
            label += f"x{scale}"
        label += f"[{h.Integral():.2f}]"
        hep.histplot(
            vals, bins=bin_edges, histtype="step",
            label=label,
            color=cfg.get("color_mplhep", "red"),
            linestyle="--", linewidth=1.5, ax=ax
        )


def draw_data(ax, data_hist, bin_edges, divide_by_bin_width, wantData=True, blind_region=None):
    if not wantData or data_hist is None:
        return None, None
    vals, errs, _, _ = get_hist_arrays(data_hist, divide_by_bin_width)

    # --- BLIND REGION ---
    if blind_region:
        if len(blind_region) == 2:
            x_min = blind_region[0]
            x_max = blind_region[1]
            mask = (bin_edges[:-1] >= x_min) & (bin_edges[:-1] < x_max)
            vals[mask] = 0.0
            errs[mask] = 0.0
    label = f"data [{data_hist.Integral():.2f}]"
    hep.histplot(vals, bins=bin_edges, yerr=errs, histtype="errorbar",
                 label=label, color="black", ax=ax)
    return vals, errs


# -----------------------------
# Ratio
# -----------------------------
def draw_ratio(ax_ratio, bin_edges, data_vals, data_errs,
               total_mc_vals, total_mc_errs, x_label, blind_region):
    if data_vals is None or total_mc_vals is None:
        return

    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.divide(data_vals, total_mc_vals,
                          out=np.zeros_like(data_vals), where=total_mc_vals != 0)
        ratio_err = np.abs(np.array(np.divide(data_errs, total_mc_vals,
                                              out=np.zeros_like(data_errs), where=total_mc_vals != 0)))
        mc_rel_unc = np.divide(total_mc_errs, total_mc_vals,
                               out=np.zeros_like(total_mc_errs), where=total_mc_vals != 0)

    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # MC band
    y_up = 1.0 + mc_rel_unc
    y_dn = np.maximum(1.0 - mc_rel_unc, 0.0)

    mask = np.ones_like(ratio, dtype=bool)
    if blind_region and len(blind_region) == 2:
        x_min, x_max = blind_region
        mask = ~((bin_centers >= x_min) & (bin_centers <= x_max))
        blind_mask = ~mask
        ratio[blind_mask] = 0.0
        y_dn[blind_mask] = 0.0
        y_up[blind_mask] = 0.0

    ax_ratio.fill_between(bin_centers, y_dn, y_up, where=y_dn > 0,
                          step="mid", facecolor="ghostwhite",
                          edgecolor="black", hatch='//', alpha=0.5, zorder=1)

    ax_ratio.errorbar(bin_centers[mask], ratio[mask], yerr=ratio_err[mask], fmt='.', color='black', markersize=10, zorder=2)

    ax_ratio.axhline(1.0, color="black", linestyle="--", linewidth=1.0)
    delta = 0.5  # fallback
    if len(ratio[mask]):
        delta = np.abs(ratio[mask] - 1).mean()

    y_max = round(1 + delta, 2)
    y_min = round(1 - delta, 2)
    # ax_ratio.set_ylim(0.9, 1.1)
    ax_ratio.set_ylim(y_min * 0.9, y_max * 1.1)
    # ax_ratio.set_ylim(0.00001, y_max * 1.2)
    ax_ratio.set_ylabel("Data/MC")
    ax_ratio.set_xlabel(x_label)

def draw_ratio_comparison(ax_ratio, bin_edges, ref_vals, ref_errs,
               reg_vals, reg_errs, x_label, y_label, blind_region, color):
    if ref_vals is None or reg_vals is None:
        return

    # print(ref_vals, reg_vals)

    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.divide(reg_vals, ref_vals,
                          out=np.zeros_like(ref_vals), where=ref_vals != 0)
        # ratio_err = np.abs(np.array(np.divide(ref_errs, reg_vals,
        #                                       out=np.zeros_like(ref_errs), where=reg_vals != 0)))

    # print(ratio)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])


    mask = np.ones_like(ratio, dtype=bool)
    # if blind_region and len(blind_region) == 2:
    #     x_min, x_max = blind_region
    #     mask = ~((bin_centers >= x_min) & (bin_centers <= x_max))
    #     blind_mask = ~mask
    #     ratio[blind_mask] = 0.0
    ax_ratio.errorbar(bin_centers[mask], ratio[mask], yerr=None, fmt='.', color=color, markersize=10, zorder=2,)
    # ax_ratio.errorbar(bin_centers[mask], ratio[mask], yerr=ratio_err[mask], fmt='.', color=color, markersize=10, zorder=2,)

    ax_ratio.axhline(1.0, color="black", linestyle="--", linewidth=1.0)
    delta = 0.5  # fallback
    if len(ratio[mask]):
        delta = np.abs(ratio[mask] - 1).mean()

    y_max = round(1 + delta, 2)
    y_min = round(1 - delta, 2)
    # print(y_max, y_min)
    ax_ratio.set_ylim(y_min * 0.8, y_max * 1.2)
    # print(f"setting y limits to {y_min * 0.1}, {y_max * 1.2}")
    # ax_ratio.set_ylim(0.00001, y_max * 1.2)
    ax_ratio.set_ylabel(y_label, fontsize=11)
    # print(f"setting x_label in {x_label}")
    ax_ratio.set_xlabel(x_label)

# -----------------------
# Main plotting function
# -----------------------
def plot_histogram_from_config(
    variable,
    histograms_dict,
    phys_model_dict,
    processes_dict,
    axes_cfg_dict,
    page_cfg_dict,
    page_cfg_custom_dict,
    filename_base,
    period,
    stacked=True,
    compare_mode=False,
    compare_vars_mode=False,
    wantLogX=False,
    wantLogY=False,
    wantData=False,
    wantSignal=False,
    wantRatio=False,
    category=None,
    channel=None,
    group_minor_contributions=False,
    minor_fraction=0.001,
    ref_region = ""
):

    var_entry = findBinEntry(axes_cfg_dict,variable)

    hist_cfg = axes_cfg_dict.get(var_entry, {})
    blind_region = hist_cfg.get("blind_region", [])
    divide_by_bin_width = bool(hist_cfg.get("divide_by_bin_width", False))
    # plot options (page-level override, then var-level)
    plot_opts = dict(page_cfg_dict.get("plot_options", {}))
    if hist_cfg.get("plot_options", None):
        plot_opts.update()

    ## KDE options, no more used
    enable_kde = False # bool(plot_opts.get("enable_kde", True))
    kde_scope = plot_opts.get("kde_scope", "total_mc")  # default
    kde_bw = plot_opts.get("kde_bw", None)
    kde_points = int(plot_opts.get("kde_points", 500))

    # Setup canvas and ratio
    canvas_size = page_cfg_dict['page_setup'].get('canvas_size', [1000, 800])
    ratio_plot = bool(wantData and wantRatio and stacked and not compare_mode) or bool(wantRatio and compare_mode)
    fig = plt.figure(figsize=(canvas_size[0] / 80, canvas_size[1] / 100))
    gs = fig.add_gridspec(
        2 if ratio_plot else 1, 1,
        height_ratios=[3, 1] if ratio_plot else [2],
        hspace=0.05 if ratio_plot else 0.25
    )
    ax = fig.add_subplot(gs[0])
    fig.subplots_adjust(top=0.85)
    ax_ratio = fig.add_subplot(gs[1], sharex=ax) if ratio_plot else None



    mc_hists = {}
    signal_hists = {}
    data_hist = None
    data_vals = data_errs = total_mc_vals = total_mc_errs = None
    y_max_comp = None
    y_maxes = []
    y_max = None


    if compare_mode:
        linestyle_cycle = [ '-', '--', ':', '-.']
        color_cycle = ['cornflowerblue', 'black', 'red', 'orange', 'gray', 'green', 'cyan', 'blue', 'magenta', 'purple']
        regions = list(histograms_dict.keys())
        ref_region = regions[0] if not ref_region else ref_region
        # print(f"using ref_region {ref_region}")
        # first_region = next(iter(histograms_dict.values()))
        ref_hist = choose_reference_binning(histograms_dict[ref_region])

        if ref_hist is None:
            print("[plot_histogram_from_config] Nessun istogramma valido per il binning in compare_mode.")
            return
        _, _, bin_edges, _ = get_hist_arrays(ref_hist, False)
        processes = list(set(phys_model_dict.get('backgrounds', []) + phys_model_dict.get('signals')))

        # --- salva i valori della regione di riferimento
        ref_vals_dict = {}
        ref_vals_err_dict = {}

        for proc in processes:

            h = histograms_dict[ref_region].get(proc)

            if h is None:
                continue

            scale = processes_dict.get(proc, {}).get("scale", 1.0)

            vals, val_errs, _, _ = get_hist_arrays(
                h,
                divide_by_bin_width,
                scale
            )
            ref_vals_dict[proc] = vals
            ref_vals_err_dict[proc] = val_errs


        for i_region, region in enumerate(regions):
            color = color_cycle[i_region % len(color_cycle)]
            region_name = folder_names.get(region, region)
            hist_dict = histograms_dict[region]

            for i_proc, proc in enumerate(processes):

                if proc not in hist_dict or hist_dict[proc] is None:
                    continue

                h = hist_dict[proc]

                scale = processes_dict.get(proc, {}).get("scale", 1.0)

                plot_vals, plot_errs, _, _ = get_hist_arrays(
                    h,
                    divide_by_bin_width,
                    scale
                )
                if plot_vals is None:
                    continue
                y_maxes.append(np.max(plot_vals))

                s_label = processes_dict.get(proc, {}).get('name', proc)
                if scale != 1.0:
                    s_label += f"x{scale}"
                s_label_name = s_label.split('_')[0]
                plot_label = f"{s_label_name}: {region_name}"

                linestyle = linestyle_cycle[i_region % len(linestyle_cycle)]

                # color = processes_dict.get(proc, {}).get("color_mplhep", "red")

                # -----------------------
                # MAIN PLOT
                # -----------------------

                if not enable_kde:

                    hep.histplot(
                        plot_vals,
                        bins=bin_edges,
                        histtype="step",
                        color=color,
                        linestyle=linestyle,
                        linewidth=2,
                        label=plot_label,
                        ax=ax
                    )

                else:

                    kde_x = None
                    kde_y = None

                    try:
                        kde_x, kde_y = compute_kde_from_hist(plot_vals, bin_edges)

                        max_hist = np.max(plot_vals)
                        max_kde = np.max(kde_y)

                        if max_kde > 0:
                            kde_y *= max_hist / max_kde

                    except Exception as exc:
                        print(f"[KDE] Warning: KDE failed for {proc} in region {region}: {exc}")

                    if kde_x is not None:

                        ax.plot(
                            kde_x,
                            kde_y,
                            color=color,
                            linestyle=linestyle,
                            linewidth=2,
                            label=plot_label
                        )

                    else:

                        hep.histplot(
                            plot_vals,
                            bins=bin_edges,
                            histtype="step",
                            color=color,
                            linestyle=linestyle,
                            linewidth=2,
                            label=plot_label,
                            ax=ax
                        )

                # -----------------------
                # RATIO PLOT
                # -----------------------

                if ratio_plot and proc in ref_vals_dict and region != ref_region:
                    # print(f"Drawing ratio comparison for process {proc} in region {region} vs ref_region {ref_region}")
                    ref_vals = ref_vals_dict[proc]
                    ref_vals_errs = ref_vals_err_dict[proc]
                    region_split = region.split("_")
                    ref_region_split = ref_region.split("_")
                    region_unique = '_'.join(k for k in list(set(region_split).difference(ref_region_split)))
                    ref_region_unique = '_'.join(k for k in list(set(ref_region_split).difference(region_split)))
                    if not ref_region_unique: ref_region_unique=ref_region
                    y_label_ratio = f"$\\frac{{{region_unique}}}{{{ref_region_unique}}}$"
                    x_label=hist_cfg.get("x_title", variable)
                    for mu_idx in [1,2]:
                        if f"mu{mu_idx}" in variable:
                            x_label = x_label.format(mu_idx=mu_idx)
                    draw_ratio_comparison(ax_ratio, bin_edges, ref_vals, ref_vals_errs,plot_vals, plot_errs, x_label, y_label_ratio, None, color)



    # -------------------------
    # Compare variables mode
    # -------------------------
    elif compare_vars_mode:
        linestyle_cycle = ['-', '--', ':', '-.', ':']
        color_cycle = ['blue', 'green', 'red', 'cyan']
        var_styles = {
            var: (linestyle_cycle[i % len(linestyle_cycle)], color_cycle[i % len(color_cycle)])
            for i, var in enumerate(histograms_dict.keys())
        }
        first_hist = next(iter(histograms_dict.values()))
        _, _, bin_edges, _ = get_hist_arrays(first_hist, False)
        linewidth = len(color_cycle) / 2
        alpha = 0.5
        all_plotted_vals = []

        for var, total_hist in histograms_dict.items():
            if total_hist is None:
                continue
            values = np.array([total_hist.GetBinContent(i + 1) for i in range(total_hist.GetNbinsX())])
            y_maxes.append(np.max(values))
            all_plotted_vals.append(values)
            style = var_styles.get(var)
            hep.histplot(
                values, bins=bin_edges, histtype="step",
                color=style[1], linewidth=linewidth, label=var, ax=ax, alpha=alpha,
            )
            linewidth -= 0.2 if linewidth > 0 else 0.1
            alpha += 0.25 / len(color_cycle)


    # -------------------------
    # Stack / standard mode
    # -------------------------
    else:
        for contrib, hist in histograms_dict.items():
            if hist is None:
                continue
            if contrib in phys_model_dict.get('data', []) + ['data']:
                data_hist = hist
            elif contrib in phys_model_dict.get('signals', []):
                signal_hists[contrib] = hist
            elif contrib in phys_model_dict.get('backgrounds', []):
                mc_hists[contrib] = hist
            else:
                print(f"ref not found for {contrib}")

        # reference binning
        ref_hist = None
        if mc_hists:
            ref_hist = next(iter(mc_hists.values()))
        elif data_hist is not None:
            ref_hist = data_hist
        elif signal_hists:
            ref_hist = next(iter(signal_hists.values()))
        else:
            ref_hist = choose_reference_binning(histograms_dict)

        if ref_hist is None:
            print("[plot_histogram_from_config] Nessun istogramma valido trovato.")
            return

        _, _, bin_edges, _ = get_hist_arrays(ref_hist, False)

        total_mc_vals = total_mc_errs = None

        # group minor contributions if requested
        mc_hists_withMinor = mc_hists.copy()
        mc_order_withMinor = []

        if group_minor_contributions and mc_hists:
            integrals = {c: mc_hists[c].Integral() for c in mc_hists}
            total = sum(integrals.values())
            threshold = minor_fraction * total
            minor_contribs = [c for c, val in integrals.items() if val < threshold]
            major_contribs = [c for c in mc_hists.keys() if c not in minor_contribs]

            if minor_contribs:
                objsToMerge = ROOT.TList()
                other_hist = mc_hists[minor_contribs[0]].Clone(f"Other_{ref_hist.GetName()}")
                for minor_contrib in minor_contribs[1:]:
                    objsToMerge.Add(mc_hists[minor_contrib])
                other_hist.Merge(objsToMerge)
                # build new dict
                mc_hists_withMinor = {c: mc_hists[c] for c in mc_hists if c not in minor_contribs}
                mc_hists_withMinor["Other"] = other_hist
                mc_order_withMinor = ["Other"] + major_contribs

        # draw MC
        if stacked:
            total_mc_vals, total_mc_errs = draw_mc_stack(
                ax, mc_hists_withMinor, processes_dict, bin_edges, divide_by_bin_width, page_cfg_dict
            )
        else:
            total_mc_vals, total_mc_errs = compute_total_mc_and_stat_err(mc_hists_withMinor, divide_by_bin_width)
            # hep.histplot(total_mc_vals, bins=bin_edges, histtype="step", alpha=0.35, ax=ax)
            for name, h in mc_hists_withMinor.items():
                vals, _, _, _ = get_hist_arrays(h, divide_by_bin_width)
                cfg = processes_dict.get(name, {})
                hep.histplot(vals, bins=bin_edges, histtype="step",
                             label=cfg.get("title", name),
                             color=cfg.get("color_mplhep", "black"), linewidth=2, ax=ax)
        # print(total_mc_vals)
        for k in mc_hists_withMinor.keys():
            y_maxes.append(mc_hists_withMinor[k].GetMaximum())
        # signals and data
        draw_signals(ax, signal_hists, processes_dict, bin_edges, divide_by_bin_width, wantSignal)

        if data_hist is not None and wantData:
            data_vals, data_errs = draw_data(ax, data_hist, bin_edges, divide_by_bin_width, wantData, blind_region)


    # -------------------------
    # Axes, scales, limits
    # -------------------------
    default_y_name = "$\\frac{\\mathrm{Events}}{\\mathrm{bin\\ width}}, \\left( \\frac{1}{\\mathrm{GeV}}\\right)$" if divide_by_bin_width else "Events"
    ax.set_ylabel(hist_cfg.get("y_title", default_y_name), fontsize=14)
    if not ratio_plot:
        x_label = hist_cfg.get("x_title", variable)
        for mu_idx in [1,2]:
            if f"mu{mu_idx}" in variable:
                x_label = x_label.format(mu_idx=mu_idx)
        ax.set_xlabel(x_label, fontsize=14)
    else:
        ax.get_xaxis().set_visible(False)
    ax.set_yscale("log" if wantLogY else "linear")
    ax.set_xscale("log" if wantLogX else "linear")

    # -------------------------
    # Plots limits: y_max
    # -------------------------

    y_max = np.max(y_maxes)
    if y_max is not None and np.isfinite(y_max) and y_max > 0:
        max_factor = hist_cfg.get("max_y_sf", 1.2) if not wantLogY else (10 ** (hist_cfg.get("max_y_sf", 1.0)))
        ax.set_ylim(top=y_max * max_factor)
        if wantLogY:
            y_min_log = max(0.01, y_max * 1e-4)
            ax.set_ylim(bottom=y_min_log)
            ax.set_ylim(bottom=0.00001)

    # -------------------------
    # Plots limits: x axis
    # -------------------------
    ax.set_xlim(bin_edges[0] * 0.99, bin_edges[-1] * 1.01)

    # -------------------------
    # Plot legend
    # -------------------------
    legend_cfg = page_cfg_dict.get("legend_mplhep", {})
    ax.legend(
        loc='upper right',
        facecolor=legend_cfg.get("fill_color", "white"),
        frameon=True,#bool(legend_cfg.get("border_size", 0) == 0),
        fontsize=legend_cfg.get("text_size", 0.02) * 60,
        framealpha=0.1,
        ncol=legend_cfg.get("ncols", 2),
        handleheight=1.5,
        labelspacing=0.1
    )


    if ratio_plot:
        if data_vals is not None and total_mc_vals is not None:

            x_label = hist_cfg.get("x_title", variable)
            for mu_idx in [1,2]:
                if f"mu{mu_idx}" in variable:
                    x_label = x_label.format(mu_idx=mu_idx)
            ax.set_xlabel(x_label, fontsize=14)
            draw_ratio(ax_ratio, bin_edges, data_vals, data_errs, total_mc_vals, total_mc_errs,
                       x_label=x_label, blind_region=blind_region)

    text_box_names = page_cfg_dict["page_setup"].get("text_boxes_mplhep", [])
    text_box_cfg = {name: page_cfg_dict.get(name, {}) for name in text_box_names}
    try:
        resolved_positions = resolve_text_positions(text_box_cfg)
    except NameError:
        resolved_positions = {}
    if not compare_mode and not compare_vars_mode:
        ax.text(0.4, 0.95, category, #','.join(el for el in list(set(histograms_dict.keys()))),
                transform=ax.transAxes,
                fontsize=11,
                verticalalignment='top',
                horizontalalignment='right')

    for name in text_box_names:
        cfg = text_box_cfg.get(name, {})
        pos = resolved_positions.get(name, cfg.get("pos", [0.02, 1.05]))
        year = period.split('_')[1]
        if year == "all": year = "2022-2023"
        if cfg.get("type") == "cms_mplhep":
            hep.cms.label(
                label=f"Preliminary", data=("data" in histograms_dict or "Data_Muon" in histograms_dict), ax=ax, loc=0,
                com=cfg.get("com", "13.6 TeV"),
                lumi=cfg.get("lumi", period_dict.get(period, "")),
                year=year,
                fontsize=cfg.get("text_size", 12)
            )
        else:
            text_content = cfg.get("text", "")
            text_content = text_content.format(category=category, channel=channel, variable=variable)
            ax.text(
                pos[0], pos[1], text_content, transform=ax.transAxes,
                fontsize=cfg.get("text_size", 10), ha="left", va="top"
            )

    plt.savefig(f"{filename_base}.pdf", bbox_inches="tight")
    print(f"Plot saved to {filename_base}.pdf")
    plt.savefig(f"{filename_base}.png", bbox_inches="tight")
    print(f"Plot saved to {filename_base}.png")
    plt.close()
