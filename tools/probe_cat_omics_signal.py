#!/usr/bin/env python3
"""Leakage-safe BLCA omics signal probe using only an outer-training partition."""

from pathlib import Path
import argparse
import json

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.preprocessing import StandardScaler
import torch


ROOT = Path(__file__).resolve().parents[1]


def comparable_c_index(times, events, risks):
    concordant = 0.0
    comparable = 0
    for i in range(len(times)):
        if events[i] < 0.5:
            continue
        later = np.flatnonzero(times > times[i])
        for j in later:
            comparable += 1
            if risks[i] > risks[j]:
                concordant += 1.0
            elif risks[i] == risks[j]:
                concordant += 0.5
    if comparable == 0:
        raise ValueError("No comparable validation pairs.")
    return concordant / comparable, comparable


def cox_loss(risk, time, event, weight, l2):
    order = torch.argsort(time, descending=True)
    ordered_risk = risk[order]
    ordered_event = event[order]
    log_denominator = torch.logcumsumexp(ordered_risk, dim=0)
    partial = ordered_risk - log_denominator
    events = ordered_event.sum().clamp_min(1.0)
    return -(partial * ordered_event).sum() / events + l2 * weight.square().mean()


def fit_cox(train_x, train_t, train_e, steps, lr, l2):
    x = torch.as_tensor(train_x, dtype=torch.float64)
    t = torch.as_tensor(train_t, dtype=torch.float64)
    e = torch.as_tensor(train_e, dtype=torch.float64)
    weight = torch.zeros(x.shape[1], dtype=torch.float64, requires_grad=True)
    bias = torch.zeros((), dtype=torch.float64, requires_grad=True)
    optimizer = torch.optim.Adam([weight, bias], lr=lr)
    for _ in range(steps):
        optimizer.zero_grad()
        risk = x @ weight + bias
        loss = cox_loss(risk, t, e, weight, l2)
        loss.backward()
        optimizer.step()
    return weight.detach().numpy(), float(bias.detach()), float(loss.detach())


def signature_features(columns, signatures):
    available = set(columns)
    selected = set()
    group_sizes = {}
    for group in signatures.columns:
        genes = signatures[group].dropna().astype(str).unique()
        expected = {
            f"{gene}_{suffix}" for gene in genes for suffix in ("rnaseq", "cnv")
        }
        matched = expected & available
        selected.update(matched)
        group_sizes[group] = len(matched)
    return sorted(selected), group_sizes


def inner_split(patient_frame, seed, validation_fraction):
    times = patient_frame["survival_months"].to_numpy(dtype=float)
    events = (1.0 - patient_frame["censorship"].to_numpy(dtype=float)).astype(int)
    time_bin = pd.qcut(times, q=4, labels=False, duplicates="drop").astype(int)
    strata = np.asarray([f"{event}_{bin_id}" for event, bin_id in zip(events, time_bin)])
    splitter = StratifiedShuffleSplit(
        n_splits=1, test_size=validation_fraction, random_state=seed
    )
    return next(splitter.split(np.zeros(len(patient_frame)), strata))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outer-fold", type=int, default=0)
    parser.add_argument("--inner-seeds", type=int, nargs="+", default=[3404, 3405])
    parser.add_argument("--validation-fraction", type=float, default=0.20)
    parser.add_argument("--pca-components", type=int, default=64)
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--lr", type=float, default=0.02)
    parser.add_argument("--l2", type=float, default=0.10)
    parser.add_argument(
        "--output-dir", default="result/cat-omics-signal-probe-blca-20260820",
    )
    args = parser.parse_args()

    cat_root = ROOT / "external/CAT-with-TPL-master"
    cat = pd.read_csv(
        cat_root / "dataset_csv/tcga_blca_all_clean.csv.zip", low_memory=False
    ).drop_duplicates("case_id").set_index("case_id")
    signatures = pd.read_csv(cat_root / "dataset_csv/signatures.csv")
    features, group_sizes = signature_features(cat.columns, signatures)
    if not features:
        raise RuntimeError("No signature-linked RNA/CNV features found.")

    outer = pd.read_csv(
        ROOT / f"data_split/5foldcv/tcga_blca/splits_{args.outer_fold}.csv",
        dtype=str,
    )
    outer_train = outer["train"].dropna().astype(str).tolist()
    if not set(outer_train).issubset(cat.index):
        raise RuntimeError("Outer-training patients are missing from the CAT omics table.")
    patients = cat.loc[outer_train].copy()
    raw_x = patients[features].to_numpy(dtype=np.float64)
    if not np.isfinite(raw_x).all():
        raise RuntimeError("Omics matrix contains non-finite values.")

    rows = []
    split_audit = {}
    for seed in args.inner_seeds:
        train_idx, val_idx = inner_split(patients, seed, args.validation_fraction)
        scaler = StandardScaler()
        train_scaled = scaler.fit_transform(raw_x[train_idx])
        val_scaled = scaler.transform(raw_x[val_idx])
        components = min(
            args.pca_components, train_scaled.shape[0] - 1, train_scaled.shape[1]
        )
        pca = PCA(n_components=components, svd_solver="randomized", random_state=seed)
        train_x = pca.fit_transform(train_scaled)
        val_x = pca.transform(val_scaled)

        train_t = patients.iloc[train_idx]["survival_months"].to_numpy(dtype=float)
        train_e = 1.0 - patients.iloc[train_idx]["censorship"].to_numpy(dtype=float)
        val_t = patients.iloc[val_idx]["survival_months"].to_numpy(dtype=float)
        val_e = 1.0 - patients.iloc[val_idx]["censorship"].to_numpy(dtype=float)
        weight, bias, final_loss = fit_cox(
            train_x, train_t, train_e, args.steps, args.lr, args.l2
        )
        train_risk = train_x @ weight + bias
        val_risk = val_x @ weight + bias
        train_c, train_pairs = comparable_c_index(train_t, train_e, train_risk)
        val_c, val_pairs = comparable_c_index(val_t, val_e, val_risk)
        rows.append(
            {
                "inner_seed": seed,
                "train_n": len(train_idx),
                "train_events": int(train_e.sum()),
                "validation_n": len(val_idx),
                "validation_events": int(val_e.sum()),
                "selected_omics": len(features),
                "pca_components": components,
                "pca_variance": float(pca.explained_variance_ratio_.sum()),
                "train_c_index": train_c,
                "validation_c_index": val_c,
                "train_comparable_pairs": train_pairs,
                "validation_comparable_pairs": val_pairs,
                "final_cox_loss": final_loss,
            }
        )
        split_audit[str(seed)] = {
            "train_patient_ids": patients.index[train_idx].tolist(),
            "validation_patient_ids": patients.index[val_idx].tolist(),
        }

    result = pd.DataFrame(rows)
    summary = {
        "outer_fold": args.outer_fold,
        "outer_training_patients_only": True,
        "outer_validation_patients_used": False,
        "selected_omics": len(features),
        "signature_group_sizes": group_sizes,
        "mean_validation_c_index": float(result["validation_c_index"].mean()),
        "min_validation_c_index": float(result["validation_c_index"].min()),
        "both_above_random": bool((result["validation_c_index"] > 0.5).all()),
    }
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = ROOT / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_dir / "probe_metrics.csv", index=False)
    (output_dir / "protocol.json").write_text(
        json.dumps({"summary": summary, "splits": split_audit}, indent=2),
        encoding="utf-8",
    )
    print(result.to_string(index=False))
    print(json.dumps(summary, indent=2))
    print(f"output_dir={output_dir}")


if __name__ == "__main__":
    main()
