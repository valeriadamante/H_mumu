# --- Helper functions for ROOT histogram manipulation (unchanged) ---
# These are crucial for handling ROOT objects within the Python environment.
import ROOT
import re
import array
import math
import numpy as np
import bisect

folder_names = {
 "inclusive_etainclusive" : "$p_T$ and $\\eta$ incl",
 "inclusive_BB" : "$p_T$ incl, BB",
 "inclusive_BO" : "$p_T$ incl, BO",
 "inclusive_BE" : "$p_T$ incl, BE",
 "inclusive_OB" : "$p_T$ incl, OB",
 "inclusive_OO" : "$p_T$ incl, OO",
 "inclusive_OE" : "$p_T$ incl, OE",
 "inclusive_EB" : "$p_T$ incl, EB",
 "inclusive_EO" : "$p_T$ incl, EO",
 "inclusive_EE" : "$p_T$ incl, EE",
 "leading_mu_pt_upto26_etainclusive" : " $p_T$ > 26 GeV, $\\eta$ incl",
 "leading_mu_pt_upto26_BB" : "$p_T$ > 26 GeV, BB",
 "leading_mu_pt_upto26_BO" : "$p_T$ > 26 GeV, BO",
 "leading_mu_pt_upto26_BE" : "$p_T$ > 26 GeV, BE",
 "leading_mu_pt_upto26_OB" : "$p_T$ > 26 GeV, OB",
 "leading_mu_pt_upto26_OO" : "$p_T$ > 26 GeV, OO",
 "leading_mu_pt_upto26_OE" : "$p_T$ > 26 GeV, OE",
 "leading_mu_pt_upto26_EB" : "$p_T$ > 26 GeV, EB",
 "leading_mu_pt_upto26_EO" : "$p_T$ > 26 GeV, EO",
 "leading_mu_pt_upto26_EE" : "$p_T$ > 26 GeV, EE",
 "leading_mu_pt_26to45_etainclusive" : "26 < $p_T$ < 45 GeV, $\\eta$ incl",
 "leading_mu_pt_26to45_BB" : " 26 < $p_T$ < 45 GeV, BB",
 "leading_mu_pt_26to45_BO" : " 26 < $p_T$ < 45 GeV, BO",
 "leading_mu_pt_26to45_BE" : " 26 < $p_T$ < 45 GeV, BE",
 "leading_mu_pt_26to45_OB" : " 26 < $p_T$ < 45 GeV, OB",
 "leading_mu_pt_26to45_OO" : " 26 < $p_T$ < 45 GeV, OO",
 "leading_mu_pt_26to45_OE" : " 26 < $p_T$ < 45 GeV, OE",
 "leading_mu_pt_26to45_EB" : " 26 < $p_T$ < 45 GeV, EB",
 "leading_mu_pt_26to45_EO" : " 26 < $p_T$ < 45 GeV, EO",
 "leading_mu_pt_26to45_EE" : " 26 < $p_T$ < 45 GeV, EE",
 "leading_mu_pt_upto45_etainclusive" : " $p_T$ < 45 GeV,  $\\eta$ incl",
 "leading_mu_pt_upto45_BB" : " $p_T$ < 45 GeV, BB",
 "leading_mu_pt_upto45_BO" : " $p_T$ < 45 GeV, BO",
 "leading_mu_pt_upto45_BE" : " $p_T$ < 45 GeV, BE",
 "leading_mu_pt_upto45_OB" : " $p_T$ < 45 GeV, OB",
 "leading_mu_pt_upto45_OO" : " $p_T$ < 45 GeV, OO",
 "leading_mu_pt_upto45_OE" : " $p_T$ < 45 GeV, OE",
 "leading_mu_pt_upto45_EB" : " $p_T$ < 45 GeV, EB",
 "leading_mu_pt_upto45_EO" : " $p_T$ < 45 GeV, EO",
 "leading_mu_pt_upto45_EE" : " $p_T$ < 45 GeV, EE",
 "leading_mu_pt_45to52_etainclusive" : " 45 < $p_T$ < 52 GeV, $\\eta$ incl",
 "leading_mu_pt_45to52_BB" : " 45 < $p_T$ < 52 GeV, BB",
 "leading_mu_pt_45to52_BO" : " 45 < $p_T$ < 52 GeV, BO",
 "leading_mu_pt_45to52_BE" : " 45 < $p_T$ < 52 GeV, BE",
 "leading_mu_pt_45to52_OB" : " 45 < $p_T$ < 52 GeV, OB",
 "leading_mu_pt_45to52_OO" : " 45 < $p_T$ < 52 GeV, OO",
 "leading_mu_pt_45to52_OE" : " 45 < $p_T$ < 52 GeV, OE",
 "leading_mu_pt_45to52_EB" : " 45 < $p_T$ < 52 GeV, EB",
 "leading_mu_pt_45to52_EO" : " 45 < $p_T$ < 52 GeV, EO",
 "leading_mu_pt_45to52_EE" : " 45 < $p_T$ < 52 GeV, EE",
 "leading_mu_pt_52to62_etainclusive" : " 52 < $p_T$ < 62 GeV, $\\eta$ incl",
 "leading_mu_pt_52to62_BB" : " 52 < $p_T$ < 62 GeV, BB",
 "leading_mu_pt_52to62_BO" : " 52 < $p_T$ < 62 GeV, BO",
 "leading_mu_pt_52to62_BE" : " 52 < $p_T$ < 62 GeV, BE",
 "leading_mu_pt_52to62_OB" : " 52 < $p_T$ < 62 GeV, OB",
 "leading_mu_pt_52to62_OO" : " 52 < $p_T$ < 62 GeV, OO",
 "leading_mu_pt_52to62_OE" : " 52 < $p_T$ < 62 GeV, OE",
 "leading_mu_pt_52to62_EB" : " 52 < $p_T$ < 62 GeV, EB",
 "leading_mu_pt_52to62_EO" : " 52 < $p_T$ < 62 GeV, EO",
 "leading_mu_pt_52to62_EE" : " 52 < $p_T$ < 62 GeV, EE",
 "leading_mu_pt_above62_etainclusive" : "$p_T$ > 62 GeV, $\\eta$ incl",
 "leading_mu_pt_above62_BB" : " $p_T$ > 62 GeV, BB",
 "leading_mu_pt_above62_BO" : " $p_T$ > 62 GeV, BO",
 "leading_mu_pt_above62_BE" : " $p_T$ > 62 GeV, BE",
 "leading_mu_pt_above62_OB" : " $p_T$ > 62 GeV, OB",
 "leading_mu_pt_above62_OO" : " $p_T$ > 62 GeV, OO",
 "leading_mu_pt_above62_OE" : " $p_T$ > 62 GeV, OE",
 "leading_mu_pt_above62_EB" : " $p_T$ > 62 GeV, EB",
 "leading_mu_pt_above62_EO" : " $p_T$ > 62 GeV, EO",
 "leading_mu_pt_above62_EE" : " $p_T$ > 62 GeV, EE",
}

period_dict = {
    "Run3_2022": "7.9804",
    "Run3_2022EE": "26.6717",
    "Run3_2023": "18.063",
    "Run3_2023BPix": "9.693",
    "Run3_all": "62.3", # sum of all periods, used for plotting
}

def findBinEntry(hist_cfg_dict, var_name):
    """
    Match variable name against regex-based histogram config entries.
    """

    matches = []

    for pattern in hist_cfg_dict.keys():
        if re.fullmatch(pattern, var_name):
            matches.append(pattern)

    if not matches:
        raise KeyError(f"No histogram config pattern matches variable '{var_name}'")

    if len(matches) > 1:
        raise RuntimeError(f"Ambiguous histogram config for '{var_name}': {matches}")

    return matches[0]

def compute_stat_unc(histograms):
    """
    Computes the statistical uncertainty for a sum of histograms.
    Returns the total histogram content and the total statistical error per bin.
    """
    if not histograms:
        return np.array([]), np.array([])

    bin_counts = [np.array([h.GetBinContent(i + 1) for i in range(h.GetNbinsX())]) for h in histograms.values()]
    bin_errors_sq = [np.array([h.GetBinError(i + 1)**2 for i in range(h.GetNbinsX())]) for h in histograms.values()]

    total_content = np.sum(bin_counts, axis=0)
    total_error = np.sqrt(np.sum(bin_errors_sq, axis=0))

    return total_content, total_error

def extract_config_for_sample(contrib, input_cfg):
    """
    Finds the configuration for a given contribution type.
    """
    for group in input_cfg:
        if "name" in group.keys() and contrib == group["name"]:
            return group
    return {}

def resolve_text_positions(text_cfgs):
    """
    Resolves relative positions of text boxes from configuration.
    """
    resolved = {}
    resolving = set()
    def resolve(name):
        if name in resolved:
            return resolved[name]
        if name in resolving:
            raise ValueError(f"Cyclic reference detected in textbox positions at '{name}'")
        if name not in text_cfgs:
            raise ValueError(f"Textbox '{name}' not found in config")

        resolving.add(name)
        cfg = text_cfgs[name]
        rel_pos = cfg.get("pos", [0, 0])
        ref = cfg.get("ref")

        if ref:
            ref_pos = resolve(ref)
            abs_pos = [ref_pos[0] + rel_pos[0], ref_pos[1] + rel_pos[1]]
        else:
            abs_pos = rel_pos

        resolved[name] = abs_pos
        resolving.remove(name)
        return abs_pos

    for name in text_cfgs:
        resolve(name)
    return resolved

def GetHistName(sample_name, sample_type, uncName, unc_scale,global_cfg_dict):
    sample_namehist = sample_type if sample_type in global_cfg_dict['sample_types_to_merge'] else sample_name
    onlyCentral = sample_name == 'data' or uncName == 'Central'
    histName = sample_namehist
    if not onlyCentral:
        histName = f"{sample_namehist}_{uncName}{unc_scale}"
    return histName


def FixNegativeContributions(histogram):
    correction_factor = 0.0

    ss_debug = ""
    ss_negative = ""

    original_Integral = histogram.Integral(0, histogram.GetNbinsX() + 1)
    ss_debug += "\nSubtracted hist for '{}'.\n".format(histogram.GetName())
    ss_debug += "Integral after bkg subtraction: {}.\n".format(original_Integral)
    if original_Integral < 0:
        print(ss_debug)
        print(
            "Integral after bkg subtraction is negative for histogram '{}'".format(
                histogram.GetName()
            )
        )
        return False, ss_debug, ss_negative

    for n in range(1, histogram.GetNbinsX() + 1):
        if histogram.GetBinContent(n) >= 0:
            continue
        prefix = (
            "WARNING"
            if histogram.GetBinContent(n) + histogram.GetBinError(n) >= 0
            else "ERROR"
        )

        ss_negative += (
            "{}: {} Bin {}, content = {}, error = {}, bin limits=[{},{}].\n".format(
                prefix,
                histogram.GetName(),
                n,
                histogram.GetBinContent(n),
                histogram.GetBinError(n),
                histogram.GetBinLowEdge(n),
                histogram.GetBinLowEdge(n + 1),
            )
        )

        error = correction_factor - histogram.GetBinContent(n)
        new_error = math.sqrt(
            math.pow(error, 2) + math.pow(histogram.GetBinError(n), 2)
        )
        histogram.SetBinContent(n, correction_factor)
        histogram.SetBinError(n, new_error)

    RenormalizeHistogram(histogram, original_Integral, True)
    return True, ss_debug, ss_negative


def RenormalizeHistogram(histogram, norm, include_overflows=True):
    integral = (
        histogram.Integral(0, histogram.GetNbinsX() + 1)
        if include_overflows
        else histogram.Integral()
    )
    if integral != 0:
        histogram.Scale(norm / integral)

def AdaptBinningToHistogram(hist, desired_binning):
    axis = hist.GetXaxis()

    # tutti gli edge originali
    original_edges = [
        axis.GetBinLowEdge(i)
        for i in range(1, axis.GetNbins() + 2)
    ]

    adapted_binning = []

    for x in desired_binning:
        idx = bisect.bisect_left(original_edges, x)

        if idx == 0:
            closest = original_edges[0]
        elif idx == len(original_edges):
            closest = original_edges[-1]
        else:
            before = original_edges[idx - 1]
            after  = original_edges[idx]
            closest = before if abs(x - before) < abs(x - after) else after

        adapted_binning.append(closest)

    # rimuovi duplicati e ordina
    adapted_binning = sorted(set(adapted_binning))

    return adapted_binning


def RebinHisto(hist_initial, new_binning, sample, wantOverflow=True, verbose=False):

    # print(f"rebinning histogram. Printing initial bin parameters")
    # print(f"initial number of bins = {hist_initial.GetNbinsX()}, "
    #       f"bin0 low edge = {hist_initial.GetXaxis().GetBinLowEdge(0)}, "
    #       f"bin(N+1) low edge = {hist_initial.GetXaxis().GetBinLowEdge(hist_initial.GetNbinsX()+1)}, "
    #       f"bin0 Up edge = {hist_initial.GetXaxis().GetBinUpEdge(0)}, "
    #       f"bin(N+1) Up edge = {hist_initial.GetXaxis().GetBinUpEdge(hist_initial.GetNbinsX()+1)}")

    # print("Requested binning:")
    # print(new_binning)

    # 🔥 adattamento automatico
    adapted_binning = AdaptBinningToHistogram(hist_initial, new_binning)

    # print("Adapted binning (ROOT-safe):")
    # print(adapted_binning)

    # sicurezza minima
    if len(adapted_binning) < 2:
        raise RuntimeError("Adapted binning has less than 2 edges!")

    new_binning_array = array.array('d', adapted_binning)

    new_hist = hist_initial.Rebin(len(adapted_binning)-1, sample, new_binning_array)

    # errori poisson per data
    if sample == 'data':
        new_hist.SetBinErrorOption(ROOT.TH1.kPoisson)

    # gestione overflow
    if wantOverflow:
        n_finalbin = new_hist.GetBinContent(new_hist.GetNbinsX())
        n_overflow = new_hist.GetBinContent(new_hist.GetNbinsX()+1)
        new_hist.SetBinContent(new_hist.GetNbinsX(), n_finalbin + n_overflow)

        err_finalbin = new_hist.GetBinError(new_hist.GetNbinsX())
        err_overflow = new_hist.GetBinError(new_hist.GetNbinsX()+1)
        new_hist.SetBinError(
            new_hist.GetNbinsX(),
            math.sqrt(err_finalbin*err_finalbin + err_overflow*err_overflow)
        )

    # # debug
    # if verbose:
    #     for nbin in range(0, new_hist.GetNbinsX()+1):
    #         print(f"nbin = {nbin}, content = {new_hist.GetBinContent(nbin)}, error {new_hist.GetBinError(nbin)}")

    # fix negativi
    fix_negative_contributions, debug_info, negative_bins_info = FixNegativeContributions(new_hist)

    return new_hist

def findNewBins(hist_cfg_dict, var, **keys):
    """
    Trova il binning per hist_cfg_dict[var] in modo ricorsivo e generico.

    Parametri:
    - hist_cfg_dict: dict con le configurazioni degli istogrammi
    - var: variabile di cui recuperare il binning
    - keys: coppie chiave=valore per channel, category, region, ecc.

    Restituisce:
    - Lista dei bin, cercando in x_rebin con match ricorsivo,
      o x_bins se non trovato.
    """
    cfg = hist_cfg_dict.get(var, {})

    # Caso base: se non esiste x_rebin, ritorna x_bins
    if 'x_rebin' not in cfg:
        return cfg.get('x_bins', [])

    x_rebin = cfg['x_rebin']

    # Caso base: se x_rebin è già una lista, ritorna direttamente
    if isinstance(x_rebin, list):
        return x_rebin

    # Ricerca ricorsiva: tenta tutte le combinazioni delle chiavi fornite
    def recursive_search(d, remaining_keys):
        # Se d è una lista, abbiamo trovato il binning
        if isinstance(d, list):
            return d
        # Se non ci sono più chiavi da provare e c'è 'other'
        if not remaining_keys and isinstance(d, dict) and 'other' in d:
            return d['other']
        if not isinstance(d, dict):
            return None

        # Prova per ogni chiave ancora disponibile se esiste nel dizionario
        for k_name, k_value in remaining_keys.items():
            if k_value in d:  # match trovato
                found = recursive_search(d[k_value], {kk: vv for kk, vv in remaining_keys.items() if kk != k_name})
                if found is not None:
                    return found

        # Se nessuna delle chiavi matcha ma esiste 'other', usalo
        if 'other' in d:
            return d['other']
        return None

    # Avvia la ricerca ricorsiva
    result = recursive_search(x_rebin, {k: v for k, v in keys.items() if v is not None})

    # Se non trovato, fallback a x_bins
    # print(result)
    return result if result is not None else cfg.get('x_bins', [])


def getNewBins(bins):
    if type(bins) == list:
        final_bins = bins
    else: # Format like "10|0:100"
        n_bins_str, bin_range = bins.split('|')
        start_str,stop_str = bin_range.split(':')
        n_bins = int(n_bins_str)
        start = float(start_str)
        stop = float(stop_str)
        bin_width = (stop - start)/n_bins
        final_bins = []
        for i in range(n_bins + 1):
            final_bins.append(start + i * bin_width)
    return final_bins

def compute_kde(data, bw_method=0.3, n_points=500, x_min=None, x_max=None):
    import numpy as np
    try:
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(data, bw_method=bw_method)
        if x_min is None: x_min = np.min(data)
        if x_max is None: x_max = np.max(data)
        xs = np.linspace(x_min, x_max, n_points)
        ys = kde(xs)
        return xs, ys
    except ImportError:
        # fallback: Gaussian smoothing manuale
        if x_min is None: x_min = np.min(data)
        if x_max is None: x_max = np.max(data)
        xs = np.linspace(x_min, x_max, n_points)
        bw = bw_method * np.std(data)
        ys = np.zeros_like(xs)
        for xi in data:
            ys += np.exp(-0.5*((xs - xi)/bw)**2)
        ys /= (len(data) * bw * np.sqrt(2*np.pi))
        return xs, ys


def get_histograms_from_dir(directory, sample_type, hist_dict, pre_path=None):
    keys = [k.GetName() for k in directory.GetListOfKeys()]
    # print(keys)
    pre_path_list = pre_path.split('/') if pre_path else []
    # print(sample_type)
    # print(sample_type in keys)
    if sample_type in keys:
        obj = directory.Get(sample_type)
        if obj.IsA().InheritsFrom(ROOT.TH1.Class()):
            obj.SetDirectory(0)

            path = directory.GetPath().split(':')[-1].strip('/')

            if path not in hist_dict:
                hist_dict[path] = {}

            if sample_type not in hist_dict[path]:
                hist_dict[path][sample_type] = obj
            else:
                hist_dict[path][sample_type].Add(obj)

    for key in keys:
        if pre_path_list and key not in pre_path_list:continue
        sub_dir = directory.Get(key)
        if sub_dir.IsA().InheritsFrom(ROOT.TDirectory.Class()):
            get_histograms_from_dir(sub_dir, sample_type,  hist_dict)



# -------------------------
# KDE utilities
# -------------------------
def _try_import_gaussian_kde():
    try:
        from scipy.stats import gaussian_kde  # type: ignore
        return gaussian_kde
    except Exception:
        return None


def kde_from_binned(vals, bin_centers, bw, n_points=400, x_min=None, x_max=None):
    """
    Build a KDE-like smooth curve from a binned histogram:
      - vals: array of counts (or density) per bin center
      - bin_centers: x positions for bins
      - bw: gaussian kernel sigma (same units as x)
      - n_points: output resolution
    Returns xs, ys where ys have same integral (sum*binwidth) as input vals*binwidth.
    """
    if len(vals) == 0:
        return np.array([]), np.array([])

    if x_min is None:
        x_min = bin_centers[0] - 0.5 * (bin_centers[1] - bin_centers[0])
    if x_max is None:
        x_max = bin_centers[-1] + 0.5 * (bin_centers[-1] - bin_centers[-2] if len(bin_centers) > 1 else 0.0)

    xs = np.linspace(x_min, x_max, n_points)
    # gaussian kernel evaluated vectorized
    # construct (n_vals x n_xs) differences
    # use broadcasting: (xs[None, :] - centers[:, None])
    sigma = float(bw)
    if sigma <= 0:
        # fallback: no smoothing -> step interpolation
        ys = np.interp(xs, bin_centers, vals, left=0, right=0)
        return xs, ys

    # compute kernel contributions
    diffs = (xs[None, :] - bin_centers[:, None]) / sigma
    kernel = np.exp(-0.5 * diffs ** 2) / (sigma * np.sqrt(2 * np.pi))
    # weight by bin content
    ys = np.dot(vals, kernel)  # shape (n_xs,)
    # normalize: ensure integral(ys dx) equals sum(vals * binwidth)
    # compute input area:
    # approximate binwidth from centers spacing:
    if len(bin_centers) > 1:
        bin_width = np.diff(np.concatenate(([bin_centers[0] - (bin_centers[1] - bin_centers[0]) / 2],
                                           0.5 * (bin_centers[1:] + bin_centers[:-1]),
                                           [bin_centers[-1] + (bin_centers[-1] - bin_centers[-2]) / 2])))[0]
        # simpler: use mean spacing
        mean_binw = np.mean(np.diff(bin_centers))
    else:
        mean_binw = 1.0
    input_area = np.sum(vals * mean_binw)
    out_area = np.trapz(ys, xs)
    if out_area > 0:
        ys *= (input_area / out_area)
    return xs, ys


def compute_kde_for_hist(hist, divide_by_bin_width=False, bw=None, n_points=500):
    """
    Compute KDE-like smooth curve for a ROOT TH1 (no raw values required).
    bw: If None -> heuristic = mean bin width * 1.0 (adjustable)
    """
    vals, _, bin_edges, bin_widths = get_hist_arrays(hist, divide_by_bin_width)
    if len(vals) == 0:
        return np.array([]), np.array([])

    # bin centers
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    # choose bw if not provided: typical choice ~ 1 * mean bin width
    if bw is None:
        bw = np.mean(bin_widths) * 1.0
    return kde_from_binned(vals, centers, bw=bw, n_points=n_points,
                           x_min=bin_edges[0], x_max=bin_edges[-1])



def compute_kde_from_hist(plot_vals, bin_edges):
    from scipy.stats import gaussian_kde
    x_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    kde = gaussian_kde(x_centers, weights=plot_vals)
    x_dense = np.linspace(bin_edges[0], bin_edges[-1], 500)
    y_dense = kde(x_dense)
    # integral_hist = np.sum(plot_vals)
    # integral_kde = np.trapz(y_dense, x_dense)

    # if integral_kde > 0:
    #     y_dense *= integral_hist / integral_kde
    return x_dense, y_dense
