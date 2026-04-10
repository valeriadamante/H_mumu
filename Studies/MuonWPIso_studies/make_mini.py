
import os
import sys
import ROOT
import glob
import argparse

sys.path.append(os.environ["ANALYSIS_PATH"])
from Studies.MuonWPIso_studies.utils import *

ROOT.EnableThreadSafety()
ROOT.EnableImplicitMT()

### setup everything

parser = argparse.ArgumentParser()
parser.add_argument("--year", required=False, default="2022")
parser.add_argument("--samples", required=False, default="")
args = parser.parse_args()
period = f"Run3_{args.year}"
setup = Setup.getGlobal(os.environ["ANALYSIS_PATH"], period,"")

### input files
for process_name,samples_list in GetAnaTuplePathDict(period, args.samples).items():
    columns = []
    print(process_name, expand_filelist(samples_list))
    inputDataFrame = ROOT.RDataFrame("Events", Utilities.ListToVector(expand_filelist(samples_list)))
    print(f"dataframe successfully opened from {samples_list}")
    Utilities.InitializeCorrections(setup, find_dataset_name(setup.datasets, process_name), stage="HistTuple")
    print("corrections initialized")
    corrections=Corrections.getGlobal()
    print("corrections get global")
    bTagAlgo = analysis.Taggers_branchesNames[
        setup.global_params.get("bTagAlgo", "particleNet")
    ]
    print(f"btagalgo = {bTagAlgo}")
    bTagWPDict = corrections.btag.getWPValues()
    print(f"bTagWPDict = {bTagWPDict}")
    miniTuple_path = miniTuple_path.format(period=period)
    print(f"processing {process_name} and using as minituple path = {miniTuple_path}")
    os.makedirs(miniTuple_path, exist_ok=True)
    output_file = f"{miniTuple_path}/mini_{process_name}.root"
    print(f"output_file is {output_file}")
    isVBF_String = 'return true' if 'VBFHto2Mu' in process_name else 'return false'
    print(f"is VBF? {isVBF_String}")
    isggH_String = 'return true' if 'GluGluHto2Mu' in process_name else 'return false'
    print(f"is ggH? {isggH_String}")
    isTT_String = 'return true' if 'TT' in process_name else 'return false'
    print(f"is TT? {isTT_String}")
    isDY_String = 'return true' if 'DY' in process_name else 'return false'
    print(f"is DY? {isDY_String}")
    isW_String = 'return true' if 'W' in process_name else 'return false'
    print(f"is W? {isW_String}")
    inputDataFrame = inputDataFrame.Define(
        "isVBF",
        isVBF_String,
    )
    columns.append("isVBF")
    inputDataFrame = inputDataFrame.Define(
        "isggH",
        isggH_String,
    )
    columns.append("isggH")
    inputDataFrame = inputDataFrame.Define(
        "isDY",
        isDY_String,
    )
    columns.append("isDY")
    inputDataFrame = inputDataFrame.Define(
        "isTT",
        isTT_String,
    )
    columns.append("isTT")
    inputDataFrame = inputDataFrame.Define(
        "isW",
        isW_String,
    )
    columns.append("isW")
    print(f"going to define all WP stuff")
    if os.path.exists(output_file):
        os.remove(output_file)
    df_singleMuonOnly,df_singleAndDiMuonTrg,columns = GetAllDfStuff(inputDataFrame,columns,bTagAlgo,bTagWPDict)
    df_singleMuonOnly.Snapshot("singleMu",output_file,columns)
    print(f"finished snapshot singleMu")
    options_update = ROOT.RDF.RSnapshotOptions()
    options_update.fMode = "UPDATE" # Opens the existing file in append mode
    df_singleAndDiMuonTrg.Snapshot("singleMuOrDiMuon",output_file,columns, options_update)
    print(f"finished snapshot singleMuOrDiMuon")
    print(f"Mini-ntupla salvata in {output_file}")
