from __future__ import annotations

import os
import sys
import tomllib

import numpy as np
import onnxruntime as ort
import psutil
import yaml

import Analysis.H_mumu as analysis
import FLAF.Common.Utilities as Utilities
from Analysis.JetRelatedFunctions import VBFNetJetCollectionDef
from Corrections.Corrections import Corrections
from FLAF.Common.Utilities import DeclareHeader


class PostVBFNetDNNProducer:
    def __init__(self, cfg, payload_name, period):
        print("VBFNet Producer init!")
        # cfg is H_mumu/configs/global.yaml
        self.period = period
        self.cfg = cfg
        self.global_cfg = self._load_global_config()
        self.payload_name = payload_name
        self._load_framework()
        self.vbf_parity, self.vbf_input_features = self._load_network_config(
            "VBFNet_models"
        )
        self.parity, self.input_features = self._load_network_config(
            "PostVBFNetDNN_models"
        )
        # Columns for tmp file
        self.vars_to_save = set(
            self.vbf_input_features + self.input_features
        ).difference(set(["VBFNet_Output"]))
        # Final columns
        self.cols_to_save = [
            f"{self.payload_name}_{col}" for col in self.cfg["columns"]
        ]
        self.vbf_models = self._load_models("VBFNet_models", self.vbf_parity)
        self.models = self._load_models("PostVBFNetDNN_models", self.parity)

    ### Init helpers ###

    def _load_framework(self):
        sys.path.append(os.environ["ANALYSIS_PATH"])
        for header in [
            "FLAF/include/Utilities.h",
            "include/Helper.h",
            "include/HmumuCore.h",
            "FLAF/include/AnalysisTools.h",
            "FLAF/include/AnalysisMath.h",
            "FLAF/include/HistHelper.h",
        ]:
            DeclareHeader(os.environ["ANALYSIS_PATH"] + "/" + header)

    def _load_global_config(self):
        filepath = os.path.join(os.environ["ANALYSIS_PATH"], "config", "global.yaml")
        with open(filepath, "r") as f:
            global_config = yaml.safe_load(f)
        period_based_config_filepath = os.path.join(
            os.environ["ANALYSIS_PATH"], "config", self.period, "global.yaml"
        )
        with open(period_based_config_filepath, "r") as f:
            global_config.update(yaml.safe_load(f))
        return global_config

    def _load_models(self, model_path, parity):
        """
        Load in the trained ONNX models
        """
        directory = os.path.join(os.environ["ANALYSIS_PATH"], "Analysis", model_path)
        models = []
        for i in range(parity):
            # The model itself
            filename = os.path.join(directory, f"trained_model_{i}.onnx")
            model = ort.InferenceSession(filename)
            models.append(model)
        return models

    def _load_network_config(self, model_path):
        """
        Read in the config used to train the network and input features used
        """
        directory = os.path.join(os.environ["ANALYSIS_PATH"], "Analysis", model_path)
        filepath = os.path.join(directory, "config.toml")
        with open(filepath, "rb") as f:
            config = tomllib.load(f)
        parity = config["splitting"]["k"]
        input_features = config["dataset"]["data_columns"]
        return parity, input_features

    ### Functions for running the inference ###

    def prepare_dfw(self, dfw, dataset_name):
        print("*********** Running prepare_dfw...")
        corrections = Corrections.getGlobal()
        dfw = analysis.DataFrameBuilderForHistograms(
            dfw.df, self.global_cfg, self.period, corrections
        )
        dfw = analysis.PrepareDFBuilder(dfw)
        dfw.df = VBFNetJetCollectionDef(dfw.df)
        return dfw

    def ApplyVBFNet(self, branches):
        print("*********** Running ApplyVBFNet...")
        nEvents = len(branches)
        event_number = np.array(getattr(branches, "FullEventId"))
        input_array = np.array(
            [
                getattr(branches, feature_name)
                for feature_name in self.vbf_input_features
            ]
        ).transpose()
        all_predictions = np.zeros([nEvents, self.vbf_parity])
        for parityIdx, sess in enumerate(self.vbf_models):
            input_name = sess.get_inputs()[0].name
            label_name = sess.get_outputs()[0].name
            predictions = sess.run([label_name], {input_name: input_array})[0]
            mask = (event_number % self.parity) != parityIdx
            predictions[mask] = 0
            all_predictions[:, parityIdx] = predictions.reshape(predictions.shape[0])
        final_predictions = np.sum(all_predictions, axis=1)
        # Last save the branches
        branches["VBFNet_Output"] = final_predictions.transpose().astype(np.float32)
        print("Finishing call, memory?")
        process = psutil.Process(os.getpid())
        mem_mb = process.memory_info().rss / 1024 / 1024
        print(f"Current memory usage: {mem_mb:.2f} MB")
        return branches

    def ApplyPostVBFNetDNN(self, branches):
        print("*********** Running ApplyPostVBFNetDNN...")
        nEvents = len(branches)
        event_number = np.array(getattr(branches, "FullEventId"))
        input_array = np.array(
            [getattr(branches, feature_name) for feature_name in self.input_features]
        ).transpose()
        all_predictions = np.zeros([nEvents, self.parity])
        for parityIdx, sess in enumerate(self.models):
            input_name = sess.get_inputs()[0].name
            label_name = sess.get_outputs()[0].name
            predictions = sess.run([label_name], {input_name: input_array})[0]
            mask = (event_number % self.parity) != parityIdx
            predictions[mask] = 0
            all_predictions[:, parityIdx] = predictions.reshape(predictions.shape[0])
        final_predictions = np.sum(all_predictions, axis=1)
        # Last save the branches
        branches["NNOutput"] = final_predictions.transpose().astype(np.float32)
        print("Finishing call, memory?")
        process = psutil.Process(os.getpid())
        mem_mb = process.memory_info().rss / 1024 / 1024
        print(f"Current memory usage: {mem_mb:.2f} MB")
        return branches

    def run(self, array):
        print("*********** Running VBFNet producer")
        array = self.ApplyVBFNet(array)
        array = self.ApplyPostVBFNetDNN(array)
        # Delete not-needed branches
        for col in array.fields:
            if col not in self.cfg["columns"]:
                if col != "FullEventId":
                    del array[col]
        # Rename the branches
        for col in self.cfg["columns"]:
            if col in array.fields:
                array[f"{self.payload_name}_{col}"] = array[f"{col}"]
                del array[f"{col}"]
            else:
                print(f"Expected column {col} not found in your payload array!")
        return array
