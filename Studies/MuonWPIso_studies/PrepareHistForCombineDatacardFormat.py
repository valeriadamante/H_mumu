import ROOT
import os
import bisect
import array
from HistogramsForCombine import run

INPUT_FILE = "DNN_NNOutput.root"
OUTPUT_FILE = "DNN_NNOutput_forCombine.root"
eos_path = "/eos/user/v/vdamante/H_mumu/NewCustomProd_MuonWPStudies_Apr8/Hists_merged/Run3_all/DNN_NNOutput/"
CHANNEL = "muMu"
BASE_DIR = f"{CHANNEL}/Signal_Fit"
CATEGORIES = ['VBF_WPcombo1', 'VBF_WPcombo2', 'VBF_WPcombo3', 'VBF_WPcombo4', 'VBF_WPcurrentcombo']#, 'VBF_WPcombo1_trueMuons', 'VBF_WPcombo2_trueMuons', 'VBF_WPcombo3_trueMuons', 'VBF_WPcombo4_trueMuons', 'VBF_WPcurrentcombo_trueMuons', 'VBF_WPcombo1_NOtrueMuons', 'VBF_WPcombo2_NOtrueMuons', 'VBF_WPcombo3_NOtrueMuons', 'VBF_WPcombo4_NOtrueMuons', 'VBF_WPcurrentcombo_NOtrueMuons','VBF_JetVeto_WPcombo1', 'VBF_JetVeto_WPcombo2', 'VBF_JetVeto_WPcombo3', 'VBF_JetVeto_WPcurrentcombo']
REGIONS = ["Signal_Fit"]
# definisci signal
SIGNALS = [
    "GluGluHto2Mu",
    "VBFHto2Mu_M125_powheg"
]

def is_signal(proc):
    return proc in SIGNALS

def get_categories(tf):
    base = tf.Get(BASE_DIR)
    print(f"Categories found in {BASE_DIR}:{[key.GetName() for key in base.GetListOfKeys()]}")
    return [key.GetName() for key in base.GetListOfKeys()]

def get_processes(tf, category):
    d = tf.Get(f"{BASE_DIR}/{category}")
    return [key.GetName() for key in d.GetListOfKeys()]

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

def compute_quantile_binning(signal_hist, n_bins):
    """
    Compute bin edges from signal quantiles
    """

    probs = [i / float(n_bins) for i in range(n_bins + 1)]
    probs_array = array.array('d', probs)

    quantiles = array.array('d', [0.] * (n_bins + 1))

    signal_hist.GetQuantiles(n_bins + 1, quantiles, probs_array)

    # force exact min/max
    # for b in range(0, signal_hist.GetNbinsX()+1):
    #     print(f"bin number= {b}, low edge= {signal_hist.GetXaxis().GetBinLowEdge(b)}, up edge= {signal_hist.GetXaxis().GetBinUpEdge(b)}, Center= {signal_hist.GetXaxis().GetBinCenter(b)}")


    quantiles[0] = signal_hist.GetXaxis().GetBinLowEdge(0)
    quantiles[-1] = signal_hist.GetXaxis().GetBinUpEdge(signal_hist.GetNbinsX()+1)
    # print(quantiles)

    # protezione contro bin duplicati
    cleaned = [quantiles[0]]
    for q in quantiles[1:]:
        if q > cleaned[-1]:
            cleaned.append(q)

    if len(cleaned) < 2:
        raise RuntimeError("Quantile binning failed: not enough unique edges")
    cleaned_new = AdaptBinningToHistogram(signal_hist, cleaned)
    print(f"old binning with quantiles = {cleaned}")
    print(f"new binning with quantiles = {cleaned_new}")
    return cleaned_new

def rebin_hist(hist, bin_edges):
    """
    Variable bin rebinning
    """
    nbins = len(bin_edges) - 1

    rebinned = hist.Rebin(
        nbins,
        hist.GetName() + "_rebin",
        array.array('d', bin_edges)
    )

    rebinned.SetDirectory(0)
    return rebinned


def make_output():
    tf_in = ROOT.TFile.Open(f"{eos_path}/{INPUT_FILE}", "READ")
    tf_out = ROOT.TFile(f"{eos_path}/{OUTPUT_FILE}", "RECREATE")
    categories = CATEGORIES
    if CATEGORIES is None:
        categories = get_categories(tf_in)

    all_processes = set()

    for cat in categories:
        tf_out.mkdir(cat)
        tf_out.cd(cat)

        d = tf_in.Get(f"{BASE_DIR}/{cat}")
        procs = get_processes(tf_in, cat)

        h_data = None
        # bins = compute_quantile_binning(tf_in.Get(f"{BASE_DIR}/VBF_WPcurrentcombo").Get("VBFHto2Mu_M125_powheg"), n_bins=12)
        bins = [0.0, 0.2199999988079071, 0.38999998569488525, 0.5099999904632568, 0.6100000143051147, 0.699999988079071, 0.7699999809265137, 0.8199999928474426, 0.8700000047683716, 0.9100000262260437, 0.9399999976158142, 0.9599999785423279, 0.9900000095367432, 1.0] # N = 12

        # bins = [0.0, 0.25999999046325684, 0.4399999976158142, 0.5799999833106995, 0.6800000071525574, 0.7699999809265137, 0.8299999833106995, 0.8899999856948853, 0.9300000071525574, 0.9599999785423279, 0.9900000095367432, 1.0] # N = 10
        # print(f"Using bin edges: {bins}")
        for proc in procs:
            h = d.Get(proc)
            h_clone = rebin_hist(h.Clone(proc), bins)
            h_clone.SetDirectory(0)

            tf_out.cd(cat)
            h_clone.Write(proc)
            if proc != "Data_Muon":
                all_processes.add(proc)
            else:
                if h_data is None:
                    h_data = rebin_hist(h.Clone("data_obs"),bins)


        tf_out.cd(cat)
        h_data.Write("data_obs")

    tf_out.Close()
    tf_in.Close()
    print(f"{eos_path}/{OUTPUT_FILE} created.")
    return categories, sorted(list(all_processes))


def write_datacard(cat, processes):
    with open(f"/afs/cern.ch/user/v/vdamante/H_mumu/Studies/MuonWPIso_studies/NewDatacards_Apr8/datacard_{cat}.txt", "w") as f:

        # header
        f.write("imax * number of channels\n")
        f.write("jmax * number of backgrounds\n")
        f.write("kmax * number of nuisance parameters\n")
        f.write("------------\n")

        f.write(f"shapes * * {eos_path}/{OUTPUT_FILE} $CHANNEL/$PROCESS\n")
        f.write("------------\n")

        # bins
        f.write("bin " + cat + "\n")
        f.write("observation " + " -1 " + "\n")
        f.write("------------\n")

        # costruisci righe
        bin_line = []
        proc_name_line = []
        proc_id_line = []
        rate_line = []

        proc_id_map = {}
        sig_id = -1
        bkg_id = 1

        for proc in processes:
            if is_signal(proc):
                proc_id_map[proc] = sig_id
                sig_id -= 1
            else:
                proc_id_map[proc] = bkg_id
                bkg_id += 1
            bin_line.append(cat)
            proc_name_line.append(proc)
            proc_id_line.append(str(proc_id_map[proc]))
            rate_line.append("-1")

        f.write("bin " + " ".join(bin_line) + "\n")
        f.write("process " + " ".join(proc_name_line) + "\n")
        f.write("process " + " ".join(proc_id_line) + "\n")
        f.write("rate " + " ".join(rate_line) + "\n")

        f.write("------------\n")
        f.write("# autoMCStats\n")

        # esempio sistematiche (placeholder)
        # f.write("lumi lnN " + " ".join(["1.025"] * len(rate_line)) + "\n")
        f.write("* autoMCStats 10")

        print(f"/afs/cern.ch/user/v/vdamante/H_mumu/Studies/MuonWPIso_studies/NewDatacards_Apr7/datacard_{cat}.txt created")


if __name__ == "__main__":


    # run(f"{eos_path}/{INPUT_FILE}", REGIONS, CATEGORIES,outfile_name_base=f"{eos_path}/DNN_NNOutput_forCombine")
    cats, procs = make_output()
    for cat in cats:
        write_datacard(cat, procs)
