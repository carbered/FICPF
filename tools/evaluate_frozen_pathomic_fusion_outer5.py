#!/usr/bin/env python3
"""Frozen outer-5 evaluation of inductively calibrated pathology-omics fusion."""

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
from probe_cat_omics_complementarity import train_empirical_cdf
from probe_cat_omics_signal import fit_cox, signature_features


ROOT = Path(__file__).resolve().parents[1]


def resolve_path(path):
    path = Path(path)
    return path if path.is_absolute() else ROOT / path


def relocate_saved_data_paths(config):
    config = dict(config)
    for key in ("data_split_path", "path_table"):
        saved = str(config[key]).replace("\\", "/")
        marker = "data_split/"
        marker_index = saved.find(marker)
        if marker_index < 0:
            continue
        local_path = ROOT / saved[marker_index:]
        if not local_path.exists():
            raise FileNotFoundError(
                f"Cannot relocate saved {key} to local project: {local_path}"
            )
        config[key] = str(local_path)
    return config


def one_path(root, pattern, description, required_suffix=None):
    paths = sorted(root.glob(pattern))
    if required_suffix is not None:
        paths = [path for path in paths if path.name.endswith(required_suffix)]
    if len(paths) != 1:
        raise RuntimeError(
            f"Expected one {description}, found {len(paths)}: {[str(p) for p in paths]}"
        )
    return paths[0]


def align_predicted_times(path, patient_ids, meta):
    frame, survival, _ = load_predictions(path)
    if set(frame["patient_id"]) != set(patient_ids):
        raise RuntimeError(f"Prediction patient IDs do not match {path}.")
    frame = frame.set_index("patient_id").loc[patient_ids]
    all_times = predicted_event_times(survival, meta.time_coordinates)
    return all_times[frame["row_index"].to_numpy(dtype=int)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default="tcga_blca")
    parser.add_argument("--baseline-root", default="result/blca-baseline-e10")
    parser.add_argument("--alpha", type=float, default=0.20)
    parser.add_argument("--pca-components", type=int, default=64)
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--lr", type=float, default=0.02)
    parser.add_argument("--l2", type=float, default=0.10)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument(
        "--output-dir", default="result/blca-frozen-pathomic-fusion-outer5-20260820"
    )
    args = parser.parse_args()

    cohort = args.dataset.replace("tcga_", "")
    baseline_root = resolve_path(args.baseline_root)
    cat_root = ROOT / "external/CAT-with-TPL-master"
    cat = pd.read_csv(
        cat_root / f"dataset_csv/tcga_{cohort}_all_clean.csv.zip",
        low_memory=False,
    ).drop_duplicates("case_id").set_index("case_id")
    signatures = pd.read_csv(cat_root / "dataset_csv/signatures.csv")
    features, group_sizes = signature_features(cat.columns, signatures)
    raw = cat[features].to_numpy(dtype=np.float64)
    raw_by_patient = {pid: raw[idx] for idx, pid in enumerate(cat.index)}

    fold_rows = []
    fold_data = []
    all_prediction_rows = []
    for fold in range(5):
        result_dir = one_path(
            baseline_root,
            f"{cohort}-VLSA-CONCH-fold_{fold}-*",
            f"{cohort} fold{fold} baseline result directory",
            required_suffix="coop_rank",
        )
        train_prediction_path = result_dir / "bag_clf_train_last_pred_train.csv"
        test_prediction_path = result_dir / "bag_clf_train_last_pred_test.csv"
        config, _ = load_config(test_prediction_path)
        config = relocate_saved_data_paths(config)
        meta = build_meta(config)
        split = pd.read_csv(ROOT / config["data_split_path"], dtype=str)
        train_ids = split["train"].dropna().tolist()
        test_column = "test" if "test" in split.columns else "val"
        test_ids = split[test_column].dropna().tolist()
        required = set(train_ids) | set(test_ids)
        if not required.issubset(raw_by_patient):
            raise RuntimeError(f"CAT omics are missing fold{fold} patients.")

        train_raw = np.stack([raw_by_patient[pid] for pid in train_ids])
        test_raw = np.stack([raw_by_patient[pid] for pid in test_ids])
        scaler = StandardScaler()
        train_scaled = scaler.fit_transform(train_raw)
        test_scaled = scaler.transform(test_raw)
        components = min(args.pca_components, len(train_ids) - 1, len(features))
        pca = PCA(
            n_components=components,
            svd_solver="randomized",
            random_state=args.seed,
        )
        train_x = pca.fit_transform(train_scaled)
        test_x = pca.transform(test_scaled)
        train_t = cat.loc[train_ids, "survival_months"].to_numpy(dtype=float)
        train_e = 1.0 - cat.loc[train_ids, "censorship"].to_numpy(dtype=float)
        weight, bias, final_loss = fit_cox(
            train_x, train_t, train_e, args.steps, args.lr, args.l2
        )
        train_omic_time = -(train_x @ weight + bias)
        test_omic_time = -(test_x @ weight + bias)

        train_path_time = align_predicted_times(
            train_prediction_path, train_ids, meta
        )
        test_path_time = align_predicted_times(test_prediction_path, test_ids, meta)
        path_calibrated = train_empirical_cdf(train_path_time, test_path_time)
        omic_calibrated = train_empirical_cdf(train_omic_time, test_omic_time)
        fused_time = (
            (1.0 - args.alpha) * path_calibrated
            + args.alpha * omic_calibrated
        )

        actual = meta.get_patient_data(pids=test_ids, ret_columns=["t", "e"])
        train_actual = meta.get_patient_data(split="train", ret_columns=["t", "e"])
        evaluation_args = (
            actual.t.to_numpy(dtype=float),
            actual.e.to_numpy(dtype=float),
            train_actual.t.to_numpy(dtype=float),
            train_actual.e.to_numpy(dtype=float),
        )
        baseline_c = c_index(test_path_time, *evaluation_args)
        omic_c = c_index(test_omic_time, *evaluation_args)
        fused_c = c_index(fused_time, *evaluation_args)
        fold_rows.append({
            "fold": fold,
            "n": len(test_ids),
            "events": int(actual.e.sum()),
            "baseline_c_index": baseline_c,
            "omics_c_index": omic_c,
            "fused_c_index": fused_c,
            "delta": fused_c - baseline_c,
            "pca_variance": float(pca.explained_variance_ratio_.sum()),
            "final_cox_loss": final_loss,
        })
        fold_data.append((test_path_time, fused_time, evaluation_args))
        all_prediction_rows.extend({
            "fold": fold,
            "patient_id": pid,
            "event_time": evaluation_args[0][idx],
            "event": evaluation_args[1][idx],
            "pathology_predicted_time": test_path_time[idx],
            "omics_predicted_time": test_omic_time[idx],
            "pathology_calibrated": path_calibrated[idx],
            "omics_calibrated": omic_calibrated[idx],
            "fused_predicted_time": fused_time[idx],
        } for idx, pid in enumerate(test_ids))

    fold_metrics = pd.DataFrame(fold_rows)
    rng = np.random.RandomState(args.seed)
    bootstrap_deltas = []
    for _ in range(args.bootstrap):
        sampled_fold_deltas = []
        for baseline_time, fused_time, evaluation_args in fold_data:
            size = len(baseline_time)
            sampled = rng.randint(0, size, size=size)
            sampled_event = evaluation_args[1][sampled]
            if not np.any(sampled_event):
                continue
            sampled_args = (
                evaluation_args[0][sampled],
                sampled_event,
                evaluation_args[2],
                evaluation_args[3],
            )
            sampled_fold_deltas.append(
                c_index(fused_time[sampled], *sampled_args)
                - c_index(baseline_time[sampled], *sampled_args)
            )
        if len(sampled_fold_deltas) == 5:
            bootstrap_deltas.append(float(np.mean(sampled_fold_deltas)))
    bootstrap_deltas = np.asarray(bootstrap_deltas)
    summary = {
        "dataset": args.dataset,
        "method": "frozen_pathology_plus_independent_pca_cox_omics",
        "alpha": args.alpha,
        "calibration_fit": "outer_training_empirical_cdf_only",
        "validation_or_test_weight_search": False,
        "selected_omics": len(features),
        "signature_group_sizes": group_sizes,
        "fold_mean_baseline_c_index": float(
            fold_metrics["baseline_c_index"].mean()
        ),
        "fold_mean_omics_c_index": float(fold_metrics["omics_c_index"].mean()),
        "fold_mean_fused_c_index": float(fold_metrics["fused_c_index"].mean()),
        "fold_mean_delta": float(fold_metrics["delta"].mean()),
        "positive_folds": int((fold_metrics["delta"] > 0).sum()),
        "bootstrap_samples": len(bootstrap_deltas),
        "delta_ci025": float(np.quantile(bootstrap_deltas, 0.025)),
        "delta_ci975": float(np.quantile(bootstrap_deltas, 0.975)),
        "probability_delta_gt_zero": float((bootstrap_deltas > 0).mean()),
    }
    output_dir = ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    fold_metrics.to_csv(output_dir / "fold_metrics.csv", index=False)
    pd.DataFrame(all_prediction_rows).to_csv(
        output_dir / "patient_predictions.csv", index=False
    )
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print("===== Frozen pathomic fusion outer-5 fold metrics =====")
    print(fold_metrics.to_string(index=False))
    print("===== Frozen pathomic fusion outer-5 summary =====")
    print(json.dumps(summary, indent=2))
    print(f"output_dir={output_dir}")


if __name__ == "__main__":
    main()
