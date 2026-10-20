"""Run repository DSMIL with the existing matched SA inputs and evaluator.

No dependencies are installed, no MOTCat process is resumed, and previous
experiments are never overwritten. Start with one fold to measure runtime.
"""

import argparse
import ast
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

from run_vlsa_ilra_baseline import (
    COHORTS, CORE_FILES, load_previous_runner, sha, write_frozen,
)


def check_motcat_inactive(expected_pid, proc_root=Path("/proc")):
    if not proc_root.is_dir():
        raise RuntimeError("Linux /proc is unavailable; cannot verify MOTCat state.")
    scripts = {"run_matched_motcat_cloud.py", "run_matched_motcat_cloud_vlsa38.py"}
    stopped, active = [], []
    # Scan every PID: the old job may have exited or restarted under a new PID.
    for proc in proc_root.iterdir():
        if not proc.name.isdigit() or not proc.is_dir():
            continue
        try:
            args = [arg.decode(errors="replace") for arg in
                    (proc / "cmdline").read_bytes().split(b"\0") if arg]
            if "run" not in args or not any(Path(arg).name in scripts for arg in args):
                continue
            state = (proc / "stat").read_text().rsplit(")", 1)[1].split()[0]
        except FileNotFoundError:
            continue
        if state in ("T", "t"):
            stopped.append(int(proc.name))
        elif state not in ("Z", "X", "x"):
            active.append((int(proc.name), state))
    if active:
        raise RuntimeError("MOTCat training is still active; do not start DSMIL: " + str(active))
    if stopped:
        print("Confirmed paused MOTCat PID(s):", sorted(stopped), flush=True)
    else:
        print("No live MOTCat training process found; no MOTCat job will be resumed.", flush=True)
    if expected_pid not in stopped:
        print("Previous MOTCat PID", expected_pid,
              "is no longer a paused MOTCat job; no signal was sent.", flush=True)


def config_path(root, cohort, fold):
    return root / "config/IFMLE/additional_vlsa_baselines/dsmil" / (
        "cfg_dsmil_conch_{}_fold{}_seed42.yaml".format(cohort, fold))


def result_path(root, cohort, fold):
    return root / "result/additional-vlsa-baselines/seed_42/dsmil" / (
        "{}-DSMIL-CONCH-fold_{}".format(cohort, fold))


def build_config(base, root, cohort, fold):
    source = base.load_yaml(base.source_config_path(cohort))
    base.validate_source_config(source, cohort)
    config = base.build_config(source, "DSMIL", cohort, fold)
    config.update({
        "save_path": "./" + result_path(root, cohort, fold).relative_to(root).as_posix(),
        "save_prediction": True,
    })
    return config


def prepare(base, root, cohorts):
    import yaml

    files = CORE_FILES + ["model/layers.py", "tools/run_matched_mil_backbones.py"]
    hashes = {name: sha(root / name) for name in files}
    hashes["deployment/driver"] = sha(Path(__file__))
    hashes["deployment/ilra_helpers"] = sha(Path(__file__).with_name("run_vlsa_ilra_baseline.py"))
    tree = ast.parse((root / "model/deepmil.py").read_text(encoding="utf-8"))
    if not any(isinstance(node, ast.ClassDef) and node.name == "DSMIL" for node in tree.body):
        raise RuntimeError("This repository has no DSMIL implementation.")
    rows, features = [], []
    for cohort in cohorts:
        features.append(base.audit_feature_coverage(cohort, verify_tensor=True))
        for fold in range(5):
            config = build_config(base, root, cohort, fold)
            path = config_path(root, cohort, fold)
            write_frozen(path, yaml.safe_dump(config, sort_keys=False))
            rows.append({"dataset": cohort, "fold": fold, "epochs": config["epochs"],
                         "config": path.relative_to(root).as_posix(),
                         "config_sha256": sha(path),
                         "split_sha256": sha(base.split_path(cohort, fold)),
                         "labels_sha256": sha(base.table_path(cohort))})
    audit = {
        "model": "repository DSMIL adapted to projected CONCH features and SA loss",
        "seed": 42, "checkpoint": "raw_last", "source_hashes": hashes,
        "training_budget": "GBMLGG 5 epochs; other cohorts 10, matching existing MIL baselines.",
        "scope": "Adapted backbone comparison, not an original DSMIL classification reproduction.",
        "selection": "Fixed configuration, no outer-test-based selection.",
        "label_processing": "Existing SA label conversion retained unchanged; no new claim of train-only time-bin fitting.",
        "features": features, "configs": rows,
    }
    path = root / "protocol" / ("additional_vlsa_dsmil_" + "_".join(cohorts) + "_audit.json")
    write_frozen(path, json.dumps(audit, indent=2))
    print("Prepared", len(rows), "frozen DSMIL configs. Audit:", path, flush=True)
    return hashes


def check_model(root, python):
    code = """import torch
from model.utils import load_model
torch.set_num_threads(2)
assert torch.cuda.is_available(), 'CUDA is unavailable'
torch.manual_seed(42)
net = load_model('DeepMIL', [512, 256, 4], network='DSMIL',
                 use_feat_proj=False, drop_rate=0.25).cuda()
for patches in (128, 2048):
    net.zero_grad()
    output = net(torch.randn(1, patches, 512, device='cuda'))
    assert tuple(output.shape) == (1, 4), tuple(output.shape)
    loss = torch.nn.functional.cross_entropy(output, torch.tensor([0], device='cuda'))
    loss.backward()
    assert torch.isfinite(loss).item()
    assert all(torch.isfinite(p.grad).all().item() for p in net.parameters() if p.grad is not None)
torch.cuda.synchronize()
print('DSMIL synthetic CUDA forward/backward: PASS. Not an experimental result.')
"""
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = "0"
    subprocess.run([python, "-c", code], cwd=str(root), env=env, check=True)


def read_completed(base, root, cohort, fold, protocol=None):
    folder = result_path(root, cohort, fold)
    files = [folder / name for name in (
        "train_metrics-last.txt", "bag_clf_train_last_pred_test.csv",
        "additional_baseline_protocol.json", "dsmil_completed.json")]
    if not all(path.is_file() for path in files):
        raise RuntimeError("No verified completion record: " + str(folder))
    record = json.loads(files[2].read_text(encoding="utf-8"))
    if protocol is not None and record != protocol:
        raise RuntimeError("Existing run protocol differs: " + str(folder))
    completion = json.loads(files[3].read_text(encoding="utf-8"))
    for path in files[:3]:
        if sha(path) != completion["sha256"][path.name]:
            raise RuntimeError("Completed file changed: " + str(path))
    score = base.read_c_index(files[0])
    if not math.isfinite(score) or not 0 <= score <= 1:
        raise ValueError("Invalid C-index: " + str(folder))
    return score, completion["wall_seconds"]


def run(base, args, root, hashes):
    log_root = root / "logs/additional_vlsa_baselines/dsmil"
    log_root.mkdir(parents=True, exist_ok=True)
    base.raise_file_limit()
    env = os.environ.copy()
    env.update({"CUDA_VISIBLE_DEVICES": "0", "MKL_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
                "PYTHONUNBUFFERED": "1", "WANDB_MODE": "offline",
                "PYTORCH_CUDA_ALLOC_CONF": "max_split_size_mb:128"})
    tasks = [(cohort, fold) for cohort in args.datasets for fold in args.folds]
    for index, (cohort, fold) in enumerate(tasks, 1):
        folder = result_path(root, cohort, fold)
        config = config_path(root, cohort, fold)
        protocol = {"model": "DSMIL", "dataset": cohort, "fold": fold, "seed": 42,
                    "config_sha256": sha(config), "source_hashes": hashes,
                    "split_sha256": sha(base.split_path(cohort, fold)),
                    "labels_sha256": sha(base.table_path(cohort)), "checkpoint": "raw_last"}
        if (folder / "dsmil_completed.json").is_file():
            score, seconds = read_completed(base, root, cohort, fold, protocol)
            print("SKIP verified completed", cohort, fold, score, seconds, flush=True)
            continue
        if folder.exists() and any(folder.iterdir()):
            raise RuntimeError("Incomplete output preserved; do not overwrite: " + str(folder))
        log = log_root / ("{}_fold{}_seed42.log".format(cohort, fold))
        if log.exists():
            raise RuntimeError("Existing log preserved: " + str(log))
        folder.mkdir(parents=True, exist_ok=True)
        record = folder / "additional_baseline_protocol.json"
        record.write_text(json.dumps(protocol, indent=2), encoding="utf-8")
        command = [args.python, "-u", "main.py", "--config", config.relative_to(root).as_posix(),
                   "--handler", "SA"]
        print("[{}/{}] RUN DSMIL {} fold {}".format(index, len(tasks), cohort, fold), flush=True)
        started = time.monotonic()
        code = base.stream_command(command, log, env)
        seconds = time.monotonic() - started
        metrics = folder / "train_metrics-last.txt"
        predictions = folder / "bag_clf_train_last_pred_test.csv"
        if code or not metrics.is_file() or not predictions.is_file():
            raise RuntimeError("DSMIL incomplete; inspect " + str(log))
        score = base.read_c_index(metrics)
        if not math.isfinite(score) or not 0 <= score <= 1:
            raise ValueError("Invalid C-index; no completion record written")
        completion = {"wall_seconds": seconds, "c_index": score,
                      "sha256": {path.name: sha(path) for path in (metrics, predictions, record)}}
        write_frozen(folder / "dsmil_completed.json", json.dumps(completion, indent=2))
        print("COMPLETE DSMIL", cohort, fold, "C-index", score,
              "wall_minutes", round(seconds / 60, 2), flush=True)
    print("All requested DSMIL folds completed or verified existing.", flush=True)


def summarize(base, root, cohorts, folds):
    import numpy as np
    import pandas as pd

    rows = []
    for cohort in cohorts:
        for fold in folds:
            score, seconds = read_completed(base, root, cohort, fold)
            rows.append({"model": "DSMIL-CONCH", "dataset": cohort.upper(), "fold": fold,
                         "c_index": score, "wall_seconds": seconds})
    frame = pd.DataFrame(rows)
    summary = frame.groupby(["model", "dataset"], as_index=False).agg(
        completed_folds=("fold", "count"), mean_c_index=("c_index", "mean"),
        std_c_index=("c_index", lambda x: float(np.std(x, ddof=0))),
        total_wall_seconds=("wall_seconds", "sum"))
    output = root / "result/additional-vlsa-baselines/summary"
    output.mkdir(parents=True, exist_ok=True)
    suffix = "_".join(cohorts) + "_folds_" + "_".join(map(str, folds))
    frame.to_csv(output / ("dsmil_fold_metrics_" + suffix + ".csv"), index=False)
    summary.to_csv(output / ("dsmil_dataset_summary_" + suffix + ".csv"), index=False)
    print(summary.to_string(index=False))
    if set(cohorts) == set(COHORTS) and set(folds) == set(range(5)):
        print("Five-cohort macro C-index:", summary.mean_c_index.mean())
    else:
        print("Partial results only; no five-cohort macro C-index.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["check", "run", "summarize"])
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--datasets", nargs="+", choices=COHORTS, default=COHORTS)
    parser.add_argument("--folds", nargs="+", type=int, choices=range(5), default=list(range(5)))
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--motcat-pid", type=int, required=True,
                        help="Previous MOTCat PID; all current MOTCat jobs are also checked.")
    args = parser.parse_args()
    if len(set(args.datasets)) != len(args.datasets) or len(set(args.folds)) != len(args.folds):
        parser.error("Repeated datasets/folds are not allowed")
    root = args.repo_root.resolve()
    check_motcat_inactive(args.motcat_pid)
    base = load_previous_runner(root)
    if args.command == "summarize":
        summarize(base, root, args.datasets, args.folds)
        return
    hashes = prepare(base, root, args.datasets)
    check_model(root, args.python)
    if args.command == "run":
        run(base, args, root, hashes)


if __name__ == "__main__":
    main()
