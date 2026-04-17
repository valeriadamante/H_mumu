import ROOT, re, array, math, numpy as np, bisect

folder_names = {
    "inclusive_etainclusive": "$p_T$ and $\\eta$ incl", "inclusive_BB": "$p_T$ incl, BB",
    "inclusive_BO": "$p_T$ incl, BO", "inclusive_BE": "$p_T$ incl, BE",
    "inclusive_OB": "$p_T$ incl, OB", "inclusive_OO": "$p_T$ incl, OO",
    "inclusive_OE": "$p_T$ incl, OE", "inclusive_EB": "$p_T$ incl, EB",
    "inclusive_EO": "$p_T$ incl, EO", "inclusive_EE": "$p_T$ incl, EE",
    "leading_mu_pt_upto26_etainclusive": " $p_T$ > 26 GeV, $\\eta$ incl",
    "leading_mu_pt_upto26_BB": "$p_T$ > 26 GeV, BB", "leading_mu_pt_upto26_BO": "$p_T$ > 26 GeV, BO",
    "leading_mu_pt_upto26_BE": "$p_T$ > 26 GeV, BE", "leading_mu_pt_upto26_OB": "$p_T$ > 26 GeV, OB",
    "leading_mu_pt_upto26_OO": "$p_T$ > 26 GeV, OO", "leading_mu_pt_upto26_OE": "$p_T$ > 26 GeV, OE",
    "leading_mu_pt_upto26_EB": "$p_T$ > 26 GeV, EB", "leading_mu_pt_upto26_EO": "$p_T$ > 26 GeV, EO",
    "leading_mu_pt_upto26_EE": "$p_T$ > 26 GeV, EE",
    "leading_mu_pt_26to45_etainclusive": "26 < $p_T$ < 45 GeV, $\\eta$ incl",
    "leading_mu_pt_26to45_BB": "26 < $p_T$ < 45 GeV, BB", "leading_mu_pt_26to45_BO": "26 < $p_T$ < 45 GeV, BO",
    "leading_mu_pt_26to45_BE": "26 < $p_T$ < 45 GeV, BE", "leading_mu_pt_26to45_OB": "26 < $p_T$ < 45 GeV, OB",
    "leading_mu_pt_26to45_OO": "26 < $p_T$ < 45 GeV, OO", "leading_mu_pt_26to45_OE": "26 < $p_T$ < 45 GeV, OE",
    "leading_mu_pt_26to45_EB": "26 < $p_T$ < 45 GeV, EB", "leading_mu_pt_26to45_EO": "26 < $p_T$ < 45 GeV, EO",
    "leading_mu_pt_26to45_EE": "26 < $p_T$ < 45 GeV, EE",
    "leading_mu_pt_upto45_etainclusive": " $p_T$ < 45 GeV,  $\\eta$ incl",
    "leading_mu_pt_upto45_BB": " $p_T$ < 45 GeV, BB", "leading_mu_pt_upto45_BO": " $p_T$ < 45 GeV, BO",
    "leading_mu_pt_upto45_BE": " $p_T$ < 45 GeV, BE", "leading_mu_pt_upto45_OB": " $p_T$ < 45 GeV, OB",
    "leading_mu_pt_upto45_OO": " $p_T$ < 45 GeV, OO", "leading_mu_pt_upto45_OE": " $p_T$ < 45 GeV, OE",
    "leading_mu_pt_upto45_EB": " $p_T$ < 45 GeV, EB", "leading_mu_pt_upto45_EO": " $p_T$ < 45 GeV, EO",
    "leading_mu_pt_upto45_EE": " $p_T$ < 45 GeV, EE",
    "leading_mu_pt_45to52_etainclusive": "45 < $p_T$ < 52 GeV, $\\eta$ incl",
    "leading_mu_pt_45to52_BB": "45 < $p_T$ < 52 GeV, BB", "leading_mu_pt_45to52_BO": "45 < $p_T$ < 52 GeV, BO",
    "leading_mu_pt_45to52_BE": "45 < $p_T$ < 52 GeV, BE", "leading_mu_pt_45to52_OB": "45 < $p_T$ < 52 GeV, OB",
    "leading_mu_pt_45to52_OO": "45 < $p_T$ < 52 GeV, OO", "leading_mu_pt_45to52_OE": "45 < $p_T$ < 52 GeV, OE",
    "leading_mu_pt_45to52_EB": "45 < $p_T$ < 52 GeV, EB", "leading_mu_pt_45to52_EO": "45 < $p_T$ < 52 GeV, EO",
    "leading_mu_pt_45to52_EE": "45 < $p_T$ < 52 GeV, EE",
    "leading_mu_pt_52to62_etainclusive": "52 < $p_T$ < 62 GeV, $\\eta$ incl",
    "leading_mu_pt_52to62_BB": "52 < $p_T$ < 62 GeV, BB", "leading_mu_pt_52to62_BO": "52 < $p_T$ < 62 GeV, BO",
    "leading_mu_pt_52to62_BE": "52 < $p_T$ < 62 GeV, BE", "leading_mu_pt_52to62_OB": "52 < $p_T$ < 62 GeV, OB",
    "leading_mu_pt_52to62_OO": "52 < $p_T$ < 62 GeV, OO", "leading_mu_pt_52to62_OE": "52 < $p_T$ < 62 GeV, OE",
    "leading_mu_pt_52to62_EB": "52 < $p_T$ < 62 GeV, EB", "leading_mu_pt_52to62_EO": "52 < $p_T$ < 62 GeV, EO",
    "leading_mu_pt_52to62_EE": "52 < $p_T$ < 62 GeV, EE",
    "leading_mu_pt_above62_etainclusive": "$p_T$ > 62 GeV, $\\eta$ incl",
    "leading_mu_pt_above62_BB": "$p_T$ > 62 GeV, BB", "leading_mu_pt_above62_BO": "$p_T$ > 62 GeV, BO",
    "leading_mu_pt_above62_BE": "$p_T$ > 62 GeV, BE", "leading_mu_pt_above62_OB": "$p_T$ > 62 GeV, OB",
    "leading_mu_pt_above62_OO": "$p_T$ > 62 GeV, OO", "leading_mu_pt_above62_OE": "$p_T$ > 62 GeV, OE",
    "leading_mu_pt_above62_EB": "$p_T$ > 62 GeV, EB", "leading_mu_pt_above62_EO": "$p_T$ > 62 GeV, EO",
    "leading_mu_pt_above62_EE": "$p_T$ > 62 GeV, EE",
}

period_dict = {"Run3_2022": "7.9804", "Run3_2022EE": "26.6717", "Run3_2023": "18.063",
               "Run3_2023BPix": "9.693", "Run3_all": "62.3"}


def findBinEntry(hist_cfg_dict, var_name):
    matches = [p for p in hist_cfg_dict.keys() if re.fullmatch(p, var_name)]
    if not matches: raise KeyError(f"No config for '{var_name}'")
    if len(matches) > 1: raise RuntimeError(f"Ambiguous: {matches}")
    return matches[0]


def get_hist_arrays(hist, divide_by_bin_width=False, scale=1.0):
    nbins = hist.GetNbinsX()
    bin_edges = np.array([hist.GetBinLowEdge(i) for i in range(1, nbins + 2)])
    bin_widths = np.array([hist.GetBinWidth(i) for i in range(1, nbins + 1)])
    vals = np.array([hist.GetBinContent(i + 1) for i in range(nbins)], dtype=float) * scale
    errs = np.array([hist.GetBinError(i + 1) for i in range(nbins)], dtype=float) * scale
    if divide_by_bin_width:
        vals = np.divide(vals, bin_widths, out=np.zeros_like(vals), where=bin_widths != 0)
        errs = np.divide(errs, bin_widths, out=np.zeros_like(errs), where=bin_widths != 0)
    return vals, errs, bin_edges, bin_widths


def calculate_bin_integrals(hist):
    nbins = hist.GetNbinsX()
    return np.array([hist.GetBinContent(i) * hist.GetBinWidth(i) for i in range(1, nbins + 1)])


def find_bin_edges(hist_list, n_new_bins, grid_size, min_bin_width, max_bin_width, xmin, xmax):
    first_hist = hist_list[0]
    nbins_old = first_hist.GetNbinsX()
    old_edges = np.array([first_hist.GetBinLowEdge(i) for i in range(1, nbins_old + 2)])
    total_bin_integrals = sum(calculate_bin_integrals(h) for h in hist_list)
    total_cumulative = np.concatenate([[0.0], np.cumsum(total_bin_integrals)])

    def snap(x): return round(x / grid_size) * grid_size
    def cum_at(x): return float(np.interp(x, old_edges, total_cumulative))
    def x_at_cum(target):
        idx = np.searchsorted(total_cumulative, target)
        if idx <= 0: return float(old_edges[0])
        if idx >= len(old_edges): return float(old_edges[-1])
        c_lo, c_hi = total_cumulative[idx - 1], total_cumulative[idx]
        x_lo, x_hi = old_edges[idx - 1], old_edges[idx]
        frac = (target - c_lo) / (c_hi - c_lo) if (c_hi - c_lo) > 0 else 0.5
        return x_lo + frac * (x_hi - x_lo)

    cum_xmin, step = cum_at(xmin), (cum_at(xmax) - cum_at(xmin)) / n_new_bins
    edges, lower_edge, k = [float(xmin)], float(xmin), 1

    while lower_edge < xmax:
        ue = snap(x_at_cum(cum_xmin + k * step))
        ue = max(ue, snap(lower_edge + min_bin_width))
        ue = min(ue, snap(lower_edge + max_bin_width))
        if ue >= xmax or ue <= lower_edge:
            edges.append(float(xmax)); break
        edges.append(ue); lower_edge = ue; k += 1

    out = [edges[0]]
    for i in range(1, len(edges)):
        w = edges[i] - out[-1]
        if w > max_bin_width:
            n_sub = int(np.ceil(w / max_bin_width))
            sw = max(snap(w / n_sub), min_bin_width)
            for _ in range(n_sub - 1): out.append(snap(out[-1] + sw))
        out.append(float(edges[i]))
    return np.array(out)


def resolve_text_positions(text_cfgs):
    # print(text_cfgs)
    resolved, resolving = {}, set()
    def resolve(name):
        if name in resolved: return resolved[name]
        if name in resolving: raise ValueError(f"Cyclic ref: '{name}'")
        cfg = text_cfgs[name]
        rel_pos, ref = cfg.get("pos", [0, 0]), cfg.get("ref")
        abs_pos = resolve(ref) + rel_pos if ref else rel_pos
        resolved[name] = abs_pos
        if name in resolving:
            resolving.remove(name)
        return abs_pos
    for name in text_cfgs: resolve(name)
    return resolved


def get_histograms_from_dir(directory, sample_type, hist_dict, pre_path=None):
    keys = [k.GetName() for k in directory.GetListOfKeys()]
    pre_path_list = pre_path.split('/') if pre_path else []
    if sample_type in keys:
        obj = directory.Get(sample_type)
        if obj.IsA().InheritsFrom(ROOT.TH1.Class()):
            obj.SetDirectory(0)
            path = directory.GetPath().split(':')[-1].strip('/')
            hist_dict.setdefault(path, {})
            if sample_type not in hist_dict[path]: hist_dict[path][sample_type] = obj
            else: hist_dict[path][sample_type].Add(obj)
    for key in keys:
        if pre_path_list and key not in pre_path_list: continue
        sub_dir = directory.Get(key)
        if sub_dir.IsA().InheritsFrom(ROOT.TDirectory.Class()):
            get_histograms_from_dir(sub_dir, sample_type, hist_dict)


def AdaptBinningToHistogram(hist, desired_binning):
    axis = hist.GetXaxis()
    original_edges = [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 2)]
    adapted = []
    for x in desired_binning:
        idx = bisect.bisect_left(original_edges, x)
        if idx == 0: closest = original_edges[0]
        elif idx == len(original_edges): closest = original_edges[-1]
        else:
            before, after = original_edges[idx - 1], original_edges[idx]
            closest = before if abs(x - before) < abs(x - after) else after
        adapted.append(closest)
    return sorted(set(adapted))


def FixNegativeContributions(histogram):
    orig_integral = histogram.Integral(0, histogram.GetNbinsX() + 1)
    if orig_integral < 0:
        print(f"Integral negative for {histogram.GetName()}")
        return False, "", ""
    for n in range(1, histogram.GetNbinsX() + 1):
        if histogram.GetBinContent(n) < 0:
            error = abs(histogram.GetBinContent(n))
            new_error = math.sqrt(error**2 + histogram.GetBinError(n)**2)
            histogram.SetBinContent(n, 0)
            histogram.SetBinError(n, new_error)
    if orig_integral > 0: histogram.Scale(1.0)
    return True, "", ""


def RebinHisto(hist_initial, new_binning, sample, wantOverflow=True, verbose=False):
    adapted = AdaptBinningToHistogram(hist_initial, new_binning)
    if len(adapted) < 2: raise RuntimeError("Adapted binning < 2 edges!")
    new_hist = hist_initial.Rebin(len(adapted) - 1, sample, array.array('d', adapted))
    if sample == 'data': new_hist.SetBinErrorOption(ROOT.TH1.kPoisson)
    if wantOverflow:
        n_final = new_hist.GetBinContent(new_hist.GetNbinsX())
        n_over = new_hist.GetBinContent(new_hist.GetNbinsX() + 1)
        new_hist.SetBinContent(new_hist.GetNbinsX(), n_final + n_over)
    FixNegativeContributions(new_hist)
    return new_hist


def findNewBins(hist_cfg_dict, var, **keys):
    cfg = hist_cfg_dict.get(var, {})
    if 'x_rebin' not in cfg: return cfg.get('x_bins', [])
    x_rebin = cfg['x_rebin']
    if isinstance(x_rebin, list): return x_rebin
    def recursive_search(d, remaining_keys):
        if isinstance(d, list): return d
        if not remaining_keys and isinstance(d, dict) and 'other' in d: return d['other']
        if not isinstance(d, dict): return None
        for k_name, k_value in remaining_keys.items():
            if k_value in d:
                found = recursive_search(d[k_value], {kk: vv for kk, vv in remaining_keys.items() if kk != k_name})
                if found is not None: return found
        return d.get('other') if isinstance(d, dict) else None
    return recursive_search(x_rebin, {k: v for k, v in keys.items() if v is not None}) or cfg.get('x_bins', [])


def getNewBins(bins):
    if isinstance(bins, list): return bins
    n_bins_str, bin_range = bins.split('|')
    start, stop = map(float, bin_range.split(':'))
    n_bins = int(n_bins_str)
    return [start + i * (stop - start) / n_bins for i in range(n_bins + 1)]


__all__ = [
    "findBinEntry", "get_hist_arrays", "calculate_bin_integrals", "find_bin_edges",
    "resolve_text_positions", "get_histograms_from_dir", "RebinHisto", "findNewBins",
    "getNewBins", "FixNegativeContributions", "AdaptBinningToHistogram",
    "folder_names", "period_dict",
]