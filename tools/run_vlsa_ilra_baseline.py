"""Run the repository ILRA backbone under the existing matched SA protocol.

This entry point reuses the previous runner, configs, main.py and evaluator.
It does not install dependencies or overwrite previous baseline experiments.
"""

import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys


COHORTS = ["blca", "brca", "gbmlgg", "luad", "ucec"]
CORE_FILES = ["main.py", "model/deepmil.py", "model/utils.py", "runner/sa_handler.py",
              "runner/base_handler.py", "dataset/PatchWSI.py", "dataset/label_converter.py",
              "loss/loss_surv.py", "loss/loss_surv_ext.py", "eval/evaluator_surv.py"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_frozen(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_text(encoding="utf-8") != content:
            raise RuntimeError("Existing frozen file differs; preserve it instead of overwriting: " + str(path))
        return
    path.write_text(content, encoding="utf-8")


def check_paused(pid):
    if pid is None:
        return
    proc = Path("/proc") / str(pid)
    command = (proc / "cmdline").read_bytes().split(b"\0")
    if not any(Path(arg.decode()).name == "run_matched_motcat_cloud_vlsa38.py" for arg in command if arg):
        raise RuntimeError("The requested paused PID is not the expected MOTCat process.")
    state = (proc / "stat").read_text().rsplit(")", 1)[1].split()[0]
    if state not in ("T", "t"):
        raise RuntimeError("MOTCat is not paused; do not start another GPU training job yet.")
    print("Confirmed paused MOTCat PID", pid, flush=True)


def load_previous_runner(root):
    path = root / "tools/run_matched_mil_backbones.py"
    spec = importlib.util.spec_from_file_location("existing_matched_runner", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.ROOT.resolve() != root:
        raise RuntimeError("Previous runner resolves to another repository.")
    return module


def config_path(root, cohort, fold):
    return root / "config/IFMLE/additional_vlsa_baselines/ilra" / ("cfg_ilra_conch_" + cohort + "_fold" + str(fold) + "_seed42.yaml")


def result_path(root, cohort, fold):
    return root / "result/additional-vlsa-baselines/seed_42/ilra" / (cohort + "-ILRA-CONCH-fold_" + str(fold))


def build_config(base, root, cohort, fold):
    source = base.load_yaml(base.source_config_path(cohort))
    base.validate_source_config(source, cohort)
    config = base.build_config(source, "ILRA", cohort, fold)
    config.update({"save_path": "./" + result_path(root, cohort, fold).relative_to(root).as_posix(),
                   "save_prediction": True, "deepmil_num_layers": 2, "deepmil_num_heads": 8,
                   "deepmil_topk": 1, "deepmil_ln": False})
    return config


def prepare(base, root, cohorts):
    import yaml
    files = CORE_FILES + ["tools/run_matched_mil_backbones.py"]
    hashes = {name: sha(root / name) for name in files}
    # Confirm architecture presence before preparing any expensive run.
    import ast
    tree = ast.parse((root / "model/deepmil.py").read_text(encoding="utf-8"))
    if not any(isinstance(node, ast.ClassDef) and node.name == "ILRA" for node in tree.body):
        raise RuntimeError("This cloud repository has no ILRA implementation.")
    feature_audits = []
    rows = []
    for cohort in cohorts:
        feature_audits.append(base.audit_feature_coverage(cohort, verify_tensor=False))
        for fold in range(5):
            config = build_config(base, root, cohort, fold)
            path = config_path(root, cohort, fold)
            write_frozen(path, yaml.safe_dump(config, sort_keys=False))
            rows.append({"dataset": cohort, "fold": fold, "epochs": config["epochs"],
                         "config": path.relative_to(root).as_posix(), "config_sha256": sha(path),
                         "split_sha256": sha(base.split_path(cohort, fold)),
                         "labels_sha256": sha(base.table_path(cohort))})
    context = {"model": "repository ILRA adapted to projected CONCH features",
               "seed": 42, "checkpoint": "raw_last", "source_hashes": hashes,
               "scope": "Matched SA backbone comparison; not an exact reproduction of the original VLSA visual-only feature setting.",
               "training_budget": "GBMLGG 5 epochs; other cohorts 10 epochs, as in existing matched ABMIL/TransMIL runs.",
               "selection": "Fixed settings, no outer-test selection.",
               "features": feature_audits, "configs": rows}
    # Different cohort selections receive separate audit files.
    audit = root / "protocol" / ("additional_vlsa_ilra_" + "_".join(cohorts) + "_audit.json")
    write_frozen(audit, json.dumps(context, indent=2))
    print("Prepared", len(rows), "frozen ILRA configs. Audit:", audit, flush=True)
    return hashes


def check_model(root, python):
    code = """import torch
from model.utils import load_model
torch.set_num_threads(2)
assert torch.cuda.is_available(), 'CUDA is unavailable'
torch.manual_seed(42)
net = load_model('DeepMIL', [512, 256, 4], network='ILRA', num_layers=2, num_heads=8, topk=1, ln=False).cuda()
output = net(torch.randn(1, 128, 512, device='cuda'))
assert tuple(output.shape) == (1, 4), tuple(output.shape)
loss = torch.nn.functional.cross_entropy(output, torch.tensor([0], device='cuda'))
loss.backward()
assert torch.isfinite(loss).item()
print('ILRA synthetic CUDA forward/backward: PASS. This is not a formal experimental result.')
"""
    environment = os.environ.copy()
    environment["CUDA_VISIBLE_DEVICES"] = "0"
    subprocess.run([python, "-c", code], cwd=str(root), env=environment, check=True)


def run(base, args, root, hashes):
    log_root = root / "logs/additional_vlsa_baselines/ilra"
    log_root.mkdir(parents=True, exist_ok=True)
    base.raise_file_limit()
    environment = os.environ.copy()
    environment.update({"CUDA_VISIBLE_DEVICES": "0", "MKL_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
                        "PYTHONUNBUFFERED": "1", "WANDB_MODE": "offline",
                        "PYTORCH_CUDA_ALLOC_CONF": "max_split_size_mb:128"})
    tasks = [(cohort, fold) for cohort in args.datasets for fold in args.folds]
    for index, (cohort, fold) in enumerate(tasks, 1):
        folder = result_path(root, cohort, fold)
        config = config_path(root, cohort, fold)
        protocol = {"model": "ILRA", "dataset": cohort, "fold": fold, "seed": 42,
                    "config_sha256": sha(config), "source_hashes": hashes,
                    "split_sha256": sha(base.split_path(cohort, fold)),
                    "labels_sha256": sha(base.table_path(cohort)),
                    "runner_sha256": sha(Path(__file__)),
                    "checkpoint": "raw_last"}
        record = folder / "additional_baseline_protocol.json"
        metrics = folder / "train_metrics-last.txt"
        predictions = folder / "bag_clf_train_last_pred_test.csv"
        if metrics.is_file() and predictions.is_file() and record.is_file():
            if json.loads(record.read_text()) != protocol:
                raise RuntimeError("Completed run protocol differs; preserve existing results: " + str(folder))
            score = base.read_c_index(metrics)
            if not math.isfinite(score) or not 0 <= score <= 1:
                raise ValueError("Invalid existing C-index")
            print("SKIP verified completed", cohort, fold, flush=True)
            continue
        if folder.exists() and any(folder.iterdir()):
            raise RuntimeError("Incomplete output exists; do not overwrite it: " + str(folder))
        folder.mkdir(parents=True, exist_ok=True)
        record.write_text(json.dumps(protocol, indent=2), encoding="utf-8")
        command = [args.python, "-u", "main.py", "--config", config.relative_to(root).as_posix(), "--handler", "SA"]
        print("[" + str(index) + "/" + str(len(tasks)) + "] RUN ILRA", cohort, "fold", fold, flush=True)
        log = log_root / (cohort + "_fold" + str(fold) + "_seed42.log")
        if log.exists():
            raise RuntimeError("Existing training log must be preserved: " + str(log))
        code = base.stream_command(command, log, environment)
        if code or not metrics.is_file() or not predictions.is_file():
            raise RuntimeError("ILRA run incomplete; inspect " + str(log))
        score = base.read_c_index(metrics)
        if not math.isfinite(score) or not 0 <= score <= 1:
            raise ValueError("Invalid completed C-index")
        print("COMPLETE ILRA", cohort, fold, "C-index", score, flush=True)
    print("All requested ILRA folds completed or verified existing.", flush=True)


def summarize(base, root, cohorts):
    import numpy as np
    import pandas as pd
    rows = []
    for cohort in cohorts:
        for fold in range(5):
            folder = result_path(root, cohort, fold)
            metrics = folder / "train_metrics-last.txt"
            if not metrics.is_file() or not (folder / "bag_clf_train_last_pred_test.csv").is_file():
                raise RuntimeError("Not all requested folds are complete; missing " + str(folder))
            score = base.read_c_index(metrics)
            if not math.isfinite(score) or not 0 <= score <= 1:
                raise ValueError("Invalid C-index in " + str(metrics))
            rows.append({"model": "ILRA-CONCH", "dataset": cohort.upper(), "fold": fold, "c_index": score})
    frame = pd.DataFrame(rows)
    summary = frame.groupby(["model", "dataset"], as_index=False).agg(
        completed_folds=("fold", "count"), mean_c_index=("c_index", "mean"),
        std_c_index=("c_index", lambda x: float(np.std(x, ddof=0))))
    output = root / "result/additional-vlsa-baselines/summary"
    output.mkdir(parents=True, exist_ok=True)
    # Deterministic derived files may be refreshed without touching source predictions.
    frame.to_csv(output / "ilra_fold_metrics.csv", index=False)
    summary.to_csv(output / "ilra_dataset_summary.csv", index=False)
    print(summary.to_string(index=False))
    if len(cohorts) == 5:
        print("Five-cohort macro C-index:", summary.mean_c_index.mean())
    print("Derived summaries:", output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["check", "run", "summarize"])
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--datasets", nargs="+", choices=COHORTS, default=COHORTS)
    parser.add_argument("--folds", nargs="+", type=int, choices=range(5), default=list(range(5)))
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--motcat-pid", type=int)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    check_paused(args.motcat_pid)
    base = load_previous_runner(root)
    if args.command == "summarize":
        summarize(base, root, args.datasets)
        return
    hashes = prepare(base, root, args.datasets)
    check_model(root, args.python)
    if args.command == "run":
        run(base, args, root, hashes)


if __name__ == "__main__":
    main()
