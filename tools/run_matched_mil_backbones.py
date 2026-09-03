#!/usr/bin/env python3
"""Prepare, audit, run, and summarize matched ABMIL/TransMIL baselines.

The experiment reuses the survival-analysis implementations and official SA
hyperparameters already shipped with VLSA. Only the WSI aggregation backbone
is varied. CONCH feature bags, outer folds, labels, and evaluation code are
shared with the matched VLSA experiment.
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
COHORTS = ("blca", "brca", "gbmlgg", "luad", "ucec")
FOLDS = tuple(range(5))
MODELS = ("ABMIL", "TransMIL")
FEATURE_DIRS = {
    "blca": Path("/root/autodl-tmp/tcga_blca/feats-CONCH-vl-proj/pt_files"),
    "brca": Path("/root/autodl-tmp/tcga_brca/feats-CONCH-vl-proj/pt_files"),
    "gbmlgg": Path("/root/autodl-tmp/tcga_gbmlgg/feats-CONCH-vl-proj/pt_files"),
    "luad": Path("/root/autodl-tmp/feats-CONCH-vl-proj/pt_files"),
    "ucec": Path("/root/autodl-tmp/tcga_ucec/feats-CONCH-vl-proj/pt_files"),
}
LOCKED_VLSA_MEANS = {
    "blca": 0.6179738476146908,
    "brca": 0.6676134405982314,
    "gbmlgg": 0.8004703042754624,
    "luad": 0.6356355680407797,
    "ucec": 0.7572158796793654,
}
METRIC_PATTERN = re.compile(
    r"^test/pred_c_index\s+-->\s+([0-9eE+.-]+)\s*$", re.MULTILINE
)


def source_config_path(cohort):
    return ROOT / "config" / "IFMLE" / f"tcga_{cohort}" / "cfg_sa_base_conch.yaml"


def generated_config_path(model, cohort, fold):
    return (
        ROOT
        / "config"
        / "IFMLE"
        / "matched_mil_backbones_20260902"
        / model.lower()
        / f"cfg_{model.lower()}_conch_{cohort}_fold{fold}_seed42.yaml"
    )


def result_directory(model, cohort, fold):
    return (
        ROOT
        / "result"
        / "matched-mil-backbones"
        / "seed_42"
        / model.lower()
        / f"{cohort}-{model}-CONCH-fold_{fold}"
    )


def table_path(cohort):
    return (
        ROOT
        / "data_split"
        / "5foldcv"
        / f"tcga_{cohort}"
        / f"mahmoodlab_tcga_{cohort}_survival.csv"
    )


def split_path(cohort, fold):
    return (
        ROOT
        / "data_split"
        / "5foldcv"
        / f"tcga_{cohort}"
        / f"splits_{fold}.csv"
    )


def load_yaml(path):
    with Path(path).open("r", encoding="utf-8") as stream:
        return yaml.safe_load(stream)


def validate_source_config(config, cohort):
    expected = {
        "task": "sa",
        "arch": "DeepMIL",
        "deepmil_network": "ABMIL",
        "deepmil_use_feat_proj": False,
        "deepmil_pooling": "attention",
        "loss_type": "SurvIFMLE",
        "evaluator": "NLL-IF",
        "net_dims": "512-256-0",
        "net_output_converter": "softmax",
        "opt_name": "adam",
        "opt_lr": 0.0002,
        "opt_weight_decay": 0.00001,
        "seed": 42,
        "ckpt_for_eval": "last",
        "es": False,
        "lrs": False,
        "test": False,
    }
    expected_epochs = 5 if cohort == "gbmlgg" else 10
    expected_bp = 16 if cohort == "gbmlgg" else 32
    expected.update({"epochs": expected_epochs, "bp_every_batch": expected_bp})
    mismatches = {
        key: {"expected": value, "actual": config.get(key)}
        for key, value in expected.items()
        if config.get(key) != value
    }
    if mismatches:
        raise AssertionError(
            f"Official SA config drift for {cohort}: {json.dumps(mismatches)}"
        )


def build_config(source, model, cohort, fold):
    config = dict(source)
    config.update(
        {
            "cuda_id": 0,
            "seed": 42,
            "wandb_dir": ".",
            "dataset_name": f"tcga_{cohort}",
            "path_patch": FEATURE_DIRS[cohort].as_posix(),
            "path_coord": None,
            "path_table": (
                f"./data_split/5foldcv/tcga_{cohort}/"
                f"mahmoodlab_tcga_{cohort}_survival.csv"
            ),
            "data_split_path": (
                f"./data_split/5foldcv/tcga_{cohort}/splits_{fold}.csv"
            ),
            "data_split_seed": fold,
            "num_shot": -1,
            "deepmil_network": model,
            "save_path": "./" + result_directory(model, cohort, fold)
            .relative_to(ROOT)
            .as_posix(),
        }
    )
    return config


def audit_feature_coverage(cohort, verify_tensor=False):
    feature_dir = FEATURE_DIRS[cohort]
    if not feature_dir.is_dir():
        raise FileNotFoundError(f"Feature directory not found: {feature_dir}")
    features = sorted(feature_dir.glob("*.pt"))
    if not features:
        raise RuntimeError(f"No .pt feature files found in {feature_dir}")
    feature_stems = {path.stem for path in features}

    table = pd.read_csv(table_path(cohort))
    required_columns = {"pathology_id", "patient_id", "e", "t"}
    missing_columns = required_columns - set(table.columns)
    if missing_columns:
        raise RuntimeError(
            f"Missing survival-table columns for {cohort}: {sorted(missing_columns)}"
        )

    split_patients = set()
    fold_counts = []
    for fold in FOLDS:
        frame = pd.read_csv(split_path(cohort, fold))
        columns = [column for column in frame.columns if not str(column).startswith("Unnamed")]
        if len(columns) < 2:
            raise RuntimeError(f"Expected train/test columns in {split_path(cohort, fold)}")
        train_ids = set(frame[columns[0]].dropna().astype(str))
        test_ids = set(frame[columns[1]].dropna().astype(str))
        overlap = train_ids & test_ids
        if overlap:
            raise RuntimeError(
                f"Train/test overlap in {cohort} fold{fold}: {sorted(overlap)[:5]}"
            )
        split_patients.update(train_ids)
        split_patients.update(test_ids)
        fold_counts.append(
            {"fold": fold, "train_patients": len(train_ids), "test_patients": len(test_ids)}
        )

    table_patients = set(table["patient_id"].astype(str))
    patients_missing_from_table = sorted(split_patients - table_patients)
    if patients_missing_from_table:
        raise RuntimeError(
            f"{cohort}: {len(patients_missing_from_table)} split patients missing from table"
        )

    relevant = table[table["patient_id"].astype(str).isin(split_patients)].copy()
    relevant["feature_available"] = relevant["pathology_id"].astype(str).isin(feature_stems)
    available_by_patient = relevant.groupby("patient_id")["feature_available"].any()
    patients_without_features = sorted(
        patient for patient in split_patients if not bool(available_by_patient.get(patient, False))
    )
    if patients_without_features:
        raise RuntimeError(
            f"{cohort}: {len(patients_without_features)} patients have no usable feature bag; "
            f"first={patients_without_features[:5]}"
        )

    sample_shape = None
    sample_dtype = None
    if verify_tensor:
        import torch

        tensor = torch.load(str(features[0]), map_location="cpu")
        if not torch.is_tensor(tensor) or tensor.ndim != 2 or tensor.shape[1] != 512:
            raise RuntimeError(
                f"Unexpected feature object in {features[0]}: "
                f"type={type(tensor)}, shape={getattr(tensor, 'shape', None)}"
            )
        if not torch.isfinite(tensor).all():
            raise RuntimeError(f"Non-finite values found in {features[0]}")
        sample_shape = list(tensor.shape)
        sample_dtype = str(tensor.dtype)

    return {
        "cohort": cohort,
        "feature_directory": feature_dir.as_posix(),
        "feature_files": len(features),
        "outer5_unique_patients": len(split_patients),
        "registered_slides": len(relevant),
        "available_registered_slides": int(relevant["feature_available"].sum()),
        "missing_registered_slides": int((~relevant["feature_available"]).sum()),
        "patients_without_features": patients_without_features,
        "folds": fold_counts,
        "sample_shape": sample_shape,
        "sample_dtype": sample_dtype,
    }


def prepare(models, cohorts, skip_feature_audit=False, verify_tensors=False):
    rows = []
    for cohort in cohorts:
        source_path = source_config_path(cohort)
        source = load_yaml(source_path)
        validate_source_config(source, cohort)
        for model in models:
            for fold in FOLDS:
                config = build_config(source, model, cohort, fold)
                output = generated_config_path(model, cohort, fold)
                output.parent.mkdir(parents=True, exist_ok=True)
                with output.open("w", encoding="utf-8") as stream:
                    yaml.safe_dump(config, stream, sort_keys=False)
                if load_yaml(output) != config:
                    raise RuntimeError(f"YAML round-trip failed: {output}")
                rows.append(
                    {
                        "model": model,
                        "cohort": cohort,
                        "fold": fold,
                        "config": output.relative_to(ROOT).as_posix(),
                        "result": result_directory(model, cohort, fold)
                        .relative_to(ROOT)
                        .as_posix(),
                    }
                )

    feature_audit = []
    if not skip_feature_audit:
        for cohort in cohorts:
            feature_audit.append(audit_feature_coverage(cohort, verify_tensors))

    audit_path = ROOT / "protocol" / "matched_mil_backbones_20260902_audit.json"
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(
        json.dumps(
            {
                "protocol": "official_vlsa_sa_backbones_matched_outer5_conch",
                "training_seed": 42,
                "checkpoint": "raw_last",
                "models": list(models),
                "cohorts": list(cohorts),
                "folds": list(FOLDS),
                "generated_configs": rows,
                "feature_audit_skipped": skip_feature_audit,
                "feature_audit": feature_audit,
                "notes": [
                    "ABMIL and TransMIL use the VLSA repository SA implementation.",
                    "Model-specific training follows cfg_sa_base_conch.yaml.",
                    "CONCH bags, patient splits, endpoints, seed, and evaluator are matched.",
                    "GBMLGG uses 5 epochs and bp_every_batch=16; other cohorts use 10 and 32.",
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print("Matched MIL backbone configs prepared.")
    print(f"generated_configs={len(rows)}")
    print(f"feature_audit_skipped={skip_feature_audit}")
    print(f"audit={audit_path}")
    return audit_path


def raise_file_limit():
    try:
        import resource

        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        target = min(max(soft, 65535), hard)
        resource.setrlimit(resource.RLIMIT_NOFILE, (target, hard))
        return target
    except (ImportError, OSError, ValueError):
        return None


def append_status(path, row):
    exists = path.exists()
    with path.open("a", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(row.keys()))
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def stream_command(command, log_path, environment):
    with log_path.open("w", encoding="utf-8", errors="replace") as log_stream:
        process = subprocess.Popen(
            command,
            cwd=str(ROOT),
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log_stream.write(line)
            log_stream.flush()
        return process.wait()


def run(args):
    prepare(
        args.models,
        args.datasets,
        skip_feature_audit=args.skip_feature_audit,
        verify_tensors=args.verify_tensors,
    )
    file_limit = raise_file_limit()
    log_root = ROOT / "logs" / "matched_mil_backbones"
    log_root.mkdir(parents=True, exist_ok=True)
    status_path = log_root / "run_status.csv"
    environment = os.environ.copy()
    environment.update(
        {
            "CUDA_VISIBLE_DEVICES": str(args.cuda_id),
            "MKL_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
            "PYTHONUNBUFFERED": "1",
            "PYTORCH_CUDA_ALLOC_CONF": "max_split_size_mb:128",
            "WANDB_MODE": "offline",
        }
    )
    tasks = [
        (model, cohort, fold)
        for model in args.models
        for cohort in args.datasets
        for fold in args.folds
    ]
    print("Matched ABMIL/TransMIL outer-fold runner.")
    print(f"root={ROOT}")
    print(f"tasks={len(tasks)}")
    print(f"models={args.models}")
    print(f"datasets={args.datasets}")
    print(f"folds={args.folds}")
    print(f"file_descriptor_limit={file_limit}")

    for index, (model, cohort, fold) in enumerate(tasks, start=1):
        config = generated_config_path(model, cohort, fold)
        metrics = result_directory(model, cohort, fold) / "train_metrics-last.txt"
        log = log_root / f"{model.lower()}_{cohort}_fold{fold}_seed42.log"
        command = [
            args.python,
            "-u",
            "main.py",
            "--config",
            config.relative_to(ROOT).as_posix(),
            "--handler",
            "SA",
        ]
        if metrics.exists() and not args.rerun:
            print(f"[{index}/{len(tasks)}] SKIP {model} {cohort} fold{fold}: metrics exist")
            status = "skipped_existing"
            return_code = 0
        elif args.dry_run:
            print(f"[{index}/{len(tasks)}] DRY-RUN " + " ".join(command))
            continue
        else:
            print(f"[{index}/{len(tasks)}] RUN " + " ".join(command))
            return_code = stream_command(command, log, environment)
            status = "completed" if return_code == 0 and metrics.exists() else "failed"
        append_status(
            status_path,
            {
                "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                "model": model,
                "cohort": cohort,
                "fold": fold,
                "status": status,
                "return_code": return_code,
                "config": config.relative_to(ROOT).as_posix(),
                "metrics": metrics.relative_to(ROOT).as_posix(),
                "log": log.relative_to(ROOT).as_posix(),
            },
        )
        if status == "failed":
            raise RuntimeError(
                f"Run failed: model={model}, cohort={cohort}, fold={fold}, log={log}"
            )
    print("All requested matched MIL backbone tasks are complete or already present.")


def read_c_index(path):
    text = path.read_text(encoding="utf-8")
    match = METRIC_PATTERN.search(text)
    if match is None:
        raise RuntimeError(f"test/pred_c_index not found in {path}")
    return float(match.group(1))


def summarize(models, cohorts, allow_incomplete=False):
    rows = []
    missing = []
    for model in models:
        for cohort in cohorts:
            for fold in FOLDS:
                metrics = result_directory(model, cohort, fold) / "train_metrics-last.txt"
                if not metrics.exists():
                    missing.append(metrics.relative_to(ROOT).as_posix())
                    continue
                rows.append(
                    {
                        "model": model,
                        "dataset": cohort.upper(),
                        "fold": fold,
                        "test_c_index": read_c_index(metrics),
                        "metrics_path": metrics.relative_to(ROOT).as_posix(),
                    }
                )
    if missing and not allow_incomplete:
        raise RuntimeError(
            f"Missing {len(missing)} metrics files; first paths:\n  " + "\n  ".join(missing[:10])
        )
    folds = pd.DataFrame(rows)
    if folds.empty:
        raise RuntimeError("No completed baseline metrics found.")
    dataset_summary = (
        folds.groupby(["model", "dataset"], as_index=False)
        .agg(
            completed_folds=("fold", "count"),
            mean_c_index=("test_c_index", "mean"),
            std_c_index=("test_c_index", lambda x: float(np.std(x, ddof=0))),
        )
        .sort_values(["model", "dataset"])
    )
    dataset_summary["matched_vlsa_mean"] = dataset_summary["dataset"].str.lower().map(
        LOCKED_VLSA_MEANS
    )
    dataset_summary["delta_vs_matched_vlsa"] = (
        dataset_summary["mean_c_index"] - dataset_summary["matched_vlsa_mean"]
    )
    complete = dataset_summary[dataset_summary["completed_folds"] == 5]
    model_summary = (
        complete.groupby("model", as_index=False)
        .agg(
            completed_datasets=("dataset", "count"),
            macro_mean_c_index=("mean_c_index", "mean"),
        )
        .sort_values("model")
    )
    model_summary["matched_vlsa_macro_mean"] = float(np.mean(list(LOCKED_VLSA_MEANS.values())))
    model_summary["delta_vs_matched_vlsa_macro"] = (
        model_summary["macro_mean_c_index"] - model_summary["matched_vlsa_macro_mean"]
    )

    output = ROOT / "result" / "matched-mil-backbones" / "summary"
    output.mkdir(parents=True, exist_ok=True)
    folds.to_csv(output / "fold_metrics.csv", index=False)
    dataset_summary.to_csv(output / "dataset_summary.csv", index=False)
    model_summary.to_csv(output / "model_summary.csv", index=False)
    (output / "summary_audit.json").write_text(
        json.dumps(
            {
                "models": list(models),
                "cohorts": list(cohorts),
                "missing_metrics": missing,
                "formal_complete": not missing,
                "matched_vlsa_macro_mean": float(np.mean(list(LOCKED_VLSA_MEANS.values()))),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print("===== Matched MIL dataset summary =====")
    print(dataset_summary.to_string(index=False))
    print("===== Matched MIL model summary =====")
    print(model_summary.to_string(index=False))
    print(f"missing_metrics={len(missing)}")
    print(f"output={output}")


def add_common_selection_arguments(parser):
    parser.add_argument("--models", nargs="+", choices=MODELS, default=list(MODELS))
    parser.add_argument("--datasets", nargs="+", choices=COHORTS, default=list(COHORTS))


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare")
    add_common_selection_arguments(prepare_parser)
    prepare_parser.add_argument("--skip-feature-audit", action="store_true")
    prepare_parser.add_argument("--verify-tensors", action="store_true")

    run_parser = subparsers.add_parser("run")
    add_common_selection_arguments(run_parser)
    run_parser.add_argument("--folds", nargs="+", type=int, choices=FOLDS, default=list(FOLDS))
    run_parser.add_argument("--python", default=sys.executable)
    run_parser.add_argument("--cuda-id", default="0")
    run_parser.add_argument("--rerun", action="store_true")
    run_parser.add_argument("--dry-run", action="store_true")
    run_parser.add_argument("--skip-feature-audit", action="store_true")
    run_parser.add_argument("--verify-tensors", action="store_true")

    summary_parser = subparsers.add_parser("summarize")
    add_common_selection_arguments(summary_parser)
    summary_parser.add_argument("--allow-incomplete", action="store_true")

    args = parser.parse_args()
    if args.command == "prepare":
        prepare(
            args.models,
            args.datasets,
            skip_feature_audit=args.skip_feature_audit,
            verify_tensors=args.verify_tensors,
        )
    elif args.command == "run":
        run(args)
    else:
        summarize(args.models, args.datasets, args.allow_incomplete)


if __name__ == "__main__":
    main()
