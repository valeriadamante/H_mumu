#!/usr/bin/env python3
import argparse
import os
import sys
import yaml
import importlib

import ROOT

if __name__ == "__main__":
    sys.path.append(os.environ["ANALYSIS_PATH"])
import FLAF.Common.Utilities as Utilities
from FLAF.Common.Setup import Setup

from Analysis.plotting_tools.drawer_functions import plot_histogram_from_config


def parse_args():
    parser = argparse.ArgumentParser(description="Plot histograms from ROOT files")
    parser.add_argument('--outFile', required=True, help="Output file name (e.g., my_plot)")
    parser.add_argument('--inFile', required=True, help="Input ROOT file")
    parser.add_argument('--period', required=True, default='Run3_2022', help="Run era")
    parser.add_argument('--var', required=True, default='m_mumu', help="Variable to plot")
    parser.add_argument('--channel', default='muMu', help="Analysis channel")
    parser.add_argument('--region', default='OS_Iso', help="Region")
    parser.add_argument('--category', default='inclusive', help="Category")
    parser.add_argument('--subregion', default=None, help="Subregion")
    parser.add_argument('--wantLogY', action='store_true', help="Log Y axis")
    parser.add_argument('--wantLogX', action='store_true', help="Log X axis")
    parser.add_argument('--wantData', action='store_true', help="Include data")
    parser.add_argument('--wantRatio', action='store_true', help="Include ratio")
    parser.add_argument('--wantSignal', action='store_true', help="Include signals")
    parser.add_argument('--contribution', default='all', help="Specific contributions (comma-separated)")
    parser.add_argument('--pre_path', default=None, help="Custom histogram path")
    parser.add_argument('--rebin', action='store_true', help="Enable rebinning")
    parser.add_argument('--compare_list', default=None, help="Comma-separated regions to compare")
    parser.add_argument('--ref_region', default=None, help="Reference region for ratio")
    parser.add_argument('--stack', action='store_true', help="Stack MC contributions")
    return parser.parse_args()


def main():
    args = parse_args()

    setup = Setup.getGlobal(
        os.environ["ANALYSIS_PATH"],
        args.period,
        None,
        custom_model_selection="WPStudiesModelPlusData"
    )

    analysis = importlib.import_module(f"{setup.global_params['analysis_import']}")

    phys_model_dict = {
        'signals': setup.phys_model.processes("signals"),
        'backgrounds': setup.phys_model.processes("backgrounds"),
        'data': setup.phys_model.processes("data"),
    }

    processes_dict = {
        proc_name: {
            "color_mplhep": proc_cfg.get("color_mplhep", "gray"),
            "scale": proc_cfg.get("scale", 1.0),
        }
        for proc_name, proc_cfg in setup.parent_processes.items()
    }

    page_cfg = os.path.join(os.environ["ANALYSIS_PATH"], 'config', 'plot', "cms_stacked.yaml")
    with open(page_cfg, 'r') as f: page_cfg_dict = yaml.safe_load(f)

    page_cfg_custom_path = os.path.join(os.environ["ANALYSIS_PATH"], 'config', 'plot', f"{args.period}.yaml")
    with open(page_cfg_custom_path, 'r') as f: page_cfg_custom_dict = yaml.safe_load(f)

    from Analysis.plotting_tools.HelpersForHistograms import findBinEntry, findNewBins, getNewBins, RebinHisto, get_histograms_from_dir

    all_contributions = phys_model_dict['backgrounds']
    if args.wantSignal: all_contributions += phys_model_dict['signals']
    if args.wantData: all_contributions += phys_model_dict['data']
    if args.contribution != 'all': all_contributions = args.contribution.split(",")

    var_entry = findBinEntry(setup.hists, args.var)

    new_bins = None
    if args.rebin and "x_rebin" in setup.hists[var_entry]:
        bins_to_compute = findNewBins(setup.hists, var_entry, channel=args.channel, category=args.category, region=args.region)
        new_bins = getNewBins(bins_to_compute)

    inFile_root = ROOT.TFile.Open(args.inFile, "READ")
    hists_to_plot = {}
    for sample_type in all_contributions:
        get_histograms_from_dir(inFile_root, sample_type, hists_to_plot)

    hists_to_plot_binned = {}
    for path, hist_dict in hists_to_plot.items():
        hists_to_plot_binned.setdefault(path, {})
        for hist_key, hist in hist_dict.items():
            if new_bins:
                hists_to_plot_binned[path][hist_key] = RebinHisto(hist, new_bins, hist_key, wantOverflow=False)
            else:
                hists_to_plot_binned[path][hist_key] = hist

    pre_path = f"{args.channel}/{args.region}/{args.category}"
    if args.pre_path: pre_path = args.pre_path
    if args.subregion: pre_path += f"/{args.subregion}"

    if args.compare_list:
        regions_to_compare = args.compare_list.split(",")
        comparison_paths = [f"{pre_path}/{reg}" for reg in regions_to_compare]
        if all(path in hists_to_plot_binned for path in comparison_paths):
            comparison_data = {reg: hists_to_plot_binned[f"{pre_path}/{reg}"] for reg in regions_to_compare}
            plot_histogram_from_config(
                variable=args.var, histograms_dict=comparison_data, phys_model_dict=phys_model_dict,
                processes_dict=processes_dict, axes_cfg_dict=setup.hists, page_cfg_dict=page_cfg_dict,
                page_cfg_custom_dict=page_cfg_custom_dict, filename_base=f"{args.outFile}_regions_comparison",
                period=args.period, stacked=args.stack, compare_mode=True, wantLogX=args.wantLogX,
                wantLogY=args.wantLogY, wantData=args.wantData, wantSignal=args.wantSignal, wantRatio=args.wantRatio,
                category=args.category, channel=args.channel, ref_region=args.ref_region or regions_to_compare[0]
            )
        else:

            print("Regions not found in histogram dict.")
    elif pre_path in hists_to_plot_binned:
        main_region_data = {contrib: hists_to_plot_binned[pre_path][contrib]
                          for contrib in all_contributions if contrib in hists_to_plot_binned[pre_path]}
        if main_region_data:
            plot_histogram_from_config(
                variable=args.var, histograms_dict=main_region_data, phys_model_dict=phys_model_dict,
                processes_dict=processes_dict, axes_cfg_dict=setup.hists, page_cfg_dict=page_cfg_dict,
                page_cfg_custom_dict=page_cfg_custom_dict, filename_base=args.outFile,
                period=args.period, stacked=args.stack, wantLogX=args.wantLogX, wantLogY=args.wantLogY,
                wantData=args.wantData, wantSignal=args.wantSignal, wantRatio=args.wantRatio,
                category=args.category, channel=args.channel
            )
        else:
            print("No contributions found")
    else:
        print(f"Path {pre_path} not found in file.")

    inFile_root.Close()


if __name__ == "__main__":
    main()