import argparse
import os
import tomllib
from glob import glob

import numpy as np
import onnxruntime as ort
from model_generation.kfold import KFolder
from model_generation.pandas_loader import PandasLoader


def get_arguments():
    """
    Builds an argument parser to get CLI arguments for the config file and dataset directory.
    """
    parser = argparse.ArgumentParser(
        prog="ONNX Inference Script",
        description="A script to run ONNX inference on a dataset using k-fold models",
    )
    parser.add_argument(
        "--directory",
        required=True,
        help="the results directory to infer on. should have a directories for each k fold model",
    )
    parser.add_argument(
        "--datafile",
        required=True,
        help="the input pandas df pkl file to use for inference",
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    # Set correct multiprocessing (needed for DataLoader parallelism)
    # multiprocessing.set_start_method("spawn", force=True)
    # Read the CLI arguments
    args = get_arguments()
    os.chdir(args.directory)

    # Read in config and datasets from args
    print("Reading config...")
    c = glob("config*.toml")[0]
    with open(c, "rb") as f:
        config = tomllib.load(f)

    # Load and split
    print("Reading in dataframe...")
    pl = PandasLoader(args.datafile, **config["dataset"])
    df = pl.load_to_dataframe()

    # Begin k-fold training loop
    print("Performing k-fold split...")
    folds = config["splitting"]["k"]
    kfolder = KFolder(k=folds, fold_idx_only=True)
    all_predictions = np.zeros(len(df))
    for k, fold_idx in enumerate(kfolder.split(df)):
        print("* ON FOLD:", k)
        test_df = df.loc[fold_idx].copy()
        os.chdir(f"{k}_fold")

        # prepare data
        x = test_df[config["dataset"]["data_columns"]].values

        # Load and run model
        filename = glob("*.onnx")[0]
        sess = ort.InferenceSession(filename)
        input_name = sess.get_inputs()[0].name
        label_name = sess.get_outputs()[0].name
        predictions = sess.run([label_name], {input_name: x})[0]
        print(predictions.shape)
        predictions = predictions.reshape(
            predictions.shape[0],
        )
        all_predictions[fold_idx] = predictions

        # Back to the top and do it all again
        os.chdir(args.directory)

    # Save the final inference plots
    print("K-fold complete! Saving final plots...")
    df["NN_Output"] = all_predictions
    df.to_pickle("ONNX_EVALUATED.pkl")
