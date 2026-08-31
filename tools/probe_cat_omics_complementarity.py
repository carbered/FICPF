#!/usr/bin/env python3
"""Probe whether CAT omics add ranking information to saved VLSA predictions."""

from pathlib import Path
import argparse
import json

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from compare_paired_survival_predictions import (
    build_meta,
    c_index,
    load_config,
    load_predictions,
    predicted_event_times,
)
from probe_cat_omics_signal import fit_cox, signature_features


ROOT = Path(__file__).resolve().parents[1]


def rank_unit(values):
    values = np.asarray(values, dtype=np.float64)
    if len(values) < 2:
        return np.zeros_like(values)
    ranks = pd.Series(values).rank(method="average").to_numpy(dtype=np.float64)
    return (ranks - 1.0) / (len(values) - 1.0)


def train_empirical_cdf(train_values, evaluation_values):
    train_values = np.sort(np.asarray(train_values, dtype=np.float64))
    evaluation_values = np.asarray(evaluation_values, dtype=np.float64)
    if len(train_values) == 0:
        raise ValueError("Training values are empty.")
    return np.searchsorted(train_values, evaluation_values, side="right") / len(
        train_values
    )


def load_ids(split_path, column):
    split = pd.read_csv(split_path, dtype=str)
    return split[column].dropna().astype(str).tolist()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default="tcga_brca")
    parser.add_argument("--split", required=True)
    parser.add_argument("--baseline-predictions", required=True)
    parser.add_argument("--baseline-train-predictions", required=True)
    parser.add_argument("--pca-components", type=int, default=64)
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--lr", type=float, default=0.02)
    parser.add_argument("--l2", type=float, default=0.10)
    parser.add_argument("--omic-weight", type=float, default=0.20)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    split_path = Path(args.split)
    if not split_path.is_absolute():
        split_path = ROOT / split_path
    prediction_path = Path(args.baseline_predictions)
    if not prediction_path.is_absolute():
        prediction_path = ROOT / prediction_path
    train_prediction_path = Path(args.baseline_train_predictions)
    if not train_prediction_path.is_absolute():
        train_prediction_path = ROOT / train_prediction_path

    cohort = args.dataset.replace("tcga_", "")
    cat_root = ROOT / "external/CAT-with-TPL-master"
    cat = pd.read_csv(
        cat_root / f"dataset_csv/tcga_{cohort}_all_clean.csv.zip",
        low_memory=False,
    ).drop_duplicates("case_id").set_index("case_id")
    signatures = pd.read_csv(cat_root / "dataset_csv/signatures.csv")
    features, group_sizes = signature_features(cat.columns, signatures)
    if not features:
        raise RuntimeError(f"No signature-linked RNA/CNV features found for {args.dataset}.")

    train_ids = load_ids(split_path, "train")
    validation_ids = load_ids(split_path, "validation")
    required = set(train_ids) | set(validation_ids)
    if not required.issubset(cat.index):
        missing = sorted(required - set(cat.index))
        raise RuntimeError(f"CAT table is missing {len(missing)} split patients: {missing[:5]}")

    train_raw = cat.loc[train_ids, features].to_numpy(dtype=np.float64)
    validation_raw = cat.loc[validation_ids, features].to_numpy(dtype=np.float64)
    scaler = StandardScaler()
    train_scaled = scaler.fit_transform(train_raw)
    validation_scaled = scaler.transform(validation_raw)
    components = min(args.pca_components, len(train_ids) - 1, len(features))
    pca = PCA(n_components=components, svd_solver="randomized", random_state=args.seed)
    train_x = pca.fit_transform(train_scaled)
    validation_x = pca.transform(validation_scaled)

    train_t = cat.loc[train_ids, "survival_months"].to_numpy(dtype=np.float64)
    train_e = 1.0 - cat.loc[train_ids, "censorship"].to_numpy(dtype=np.float64)
    validation_t = cat.loc[validation_ids, "survival_months"].to_numpy(dtype=np.float64)
    validation_e = 1.0 - cat.loc[validation_ids, "censorship"].to_numpy(dtype=np.float64)
    weight, bias, final_loss = fit_cox(
        train_x, train_t, train_e, args.steps, args.lr, args.l2
    )
    omic_risk = validation_x @ weight + bias

    baseline, baseline_survival, _ = load_predictions(prediction_path)
    if set(baseline["patient_id"]) != set(validation_ids):
        raise RuntimeError("Saved baseline predictions do not exactly match validation IDs.")
    baseline = baseline.set_index("patient_id").loc[validation_ids]
    if not np.allclose(baseline["e"].to_numpy(dtype=float), validation_e):
        raise RuntimeError("Baseline event labels do not match the CAT omics table.")
    config, _ = load_config(prediction_path)
    meta = build_meta(config)
    actual = meta.get_patient_data(pids=validation_ids, ret_columns=["t", "e"])
    train_meta = meta.get_patient_data(split="train", ret_columns=["t", "e"])
    event_times = actual.t.to_numpy(dtype=np.float64)
    event_indicators = actual.e.to_numpy(dtype=np.float64)
    if not np.allclose(event_indicators, validation_e):
        raise RuntimeError("SurvivalEVAL and CAT event labels do not match.")
    baseline_times_all = predicted_event_times(baseline_survival, meta.time_coordinates)
    baseline_times = baseline_times_all[baseline["row_index"].to_numpy(dtype=int)]
    baseline_train, baseline_train_survival, _ = load_predictions(train_prediction_path)
    if set(baseline_train["patient_id"]) != set(train_ids):
        raise RuntimeError("Saved baseline training predictions do not match train IDs.")
    baseline_train = baseline_train.set_index("patient_id").loc[train_ids]
    baseline_train_times_all = predicted_event_times(
        baseline_train_survival, meta.time_coordinates
    )
    baseline_train_times = baseline_train_times_all[
        baseline_train["row_index"].to_numpy(dtype=int)
    ]
    train_omic_risk = train_x @ weight + bias
    omic_predicted_times = -omic_risk
    train_omic_predicted_times = -train_omic_risk
    alpha = float(args.omic_weight)
    fixed_blend = (
        (1.0 - alpha)
        * train_empirical_cdf(baseline_train_times, baseline_times)
        + alpha
        * train_empirical_cdf(train_omic_predicted_times, omic_predicted_times)
    )

    evaluation_args = (
        event_times,
        event_indicators,
        train_meta.t.to_numpy(dtype=np.float64),
        train_meta.e.to_numpy(dtype=np.float64),
    )
    baseline_c = c_index(baseline_times, *evaluation_args)
    omic_c = c_index(omic_predicted_times, *evaluation_args)
    blend_c = c_index(fixed_blend, *evaluation_args)
    output = {
        "dataset": args.dataset,
        "seed": args.seed,
        "train_n": len(train_ids),
        "train_events": int(train_e.sum()),
        "validation_n": len(validation_ids),
        "validation_events": int(validation_e.sum()),
        "metric_contract": "SurvivalEVAL_mean_continuous_comparable_all_ties",
        "selected_omics": len(features),
        "signature_group_sizes": group_sizes,
        "pca_components": components,
        "pca_variance": float(pca.explained_variance_ratio_.sum()),
        "final_cox_loss": final_loss,
        "baseline_c_index": baseline_c,
        "omic_c_index": omic_c,
        "fixed_omic_weight": alpha,
        "fixed_blend_c_index": blend_c,
        "fixed_blend_delta": blend_c - baseline_c,
        "train_only_preprocessing": True,
        "validation_weight_search": False,
        "calibration_fit": "training_empirical_cdf_only",
    }
    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = ROOT / output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))
    print(f"output={output_path}")


if __name__ == "__main__":
    main()
