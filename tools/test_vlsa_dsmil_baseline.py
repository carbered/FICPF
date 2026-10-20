"""Local safety/configuration tests; no cloud data or formal training is used."""

import ast
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import run_vlsa_dsmil_baseline as driver


class DSMILDriverTests(unittest.TestCase):
    def fake_process(self, root, pid, state="T", script="run_matched_motcat_cloud_vlsa38.py"):
        proc = root / str(pid)
        proc.mkdir()
        (proc / "cmdline").write_bytes(
            ("python\0tools/" + script + "\0run\0--datasets\0blca\0").encode())
        (proc / "stat").write_text("{} (python) {} 1 2 3".format(pid, state))

    def test_absent_old_pid_is_allowed_after_scan(self):
        with TemporaryDirectory() as tmp, redirect_stdout(io.StringIO()):
            driver.check_motcat_inactive(3021, Path(tmp))

    def test_paused_job_with_new_pid_is_allowed(self):
        with TemporaryDirectory() as tmp, redirect_stdout(io.StringIO()):
            root = Path(tmp)
            self.fake_process(root, 9000)
            driver.check_motcat_inactive(3021, root)

    def test_active_job_with_new_pid_is_blocked(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.fake_process(root, 9000, "R")
            with self.assertRaisesRegex(RuntimeError, "still active"):
                driver.check_motcat_inactive(3021, root)

    def test_reused_old_pid_does_not_hide_other_active_job(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.fake_process(root, 3021, "R", "unrelated.py")
            self.fake_process(root, 9000, "S", "run_matched_motcat_cloud.py")
            with self.assertRaisesRegex(RuntimeError, "still active"):
                driver.check_motcat_inactive(3021, root)

    def test_paused_old_pid_does_not_hide_other_active_job(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.fake_process(root, 3021)
            self.fake_process(root, 9000, "D")
            with self.assertRaisesRegex(RuntimeError, "still active"):
                driver.check_motcat_inactive(3021, root)

    def test_missing_proc_is_blocked(self):
        with TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(RuntimeError, "unavailable"):
                driver.check_motcat_inactive(3021, Path(tmp) / "absent")

    def test_config_only_changes_backbone_and_output(self):
        root = Path("/example")
        validated = []
        source = {"epochs": 10, "loss_type": "SurvIFMLE", "deepmil_use_feat_proj": False}
        expected = dict(source, deepmil_network="DSMIL", save_path="old")
        base = SimpleNamespace(
            source_config_path=lambda cohort: cohort,
            load_yaml=lambda path: source,
            validate_source_config=lambda config, cohort: validated.append(cohort),
            build_config=lambda config, model, cohort, fold: dict(expected),
        )
        config = driver.build_config(base, root, "blca", 0)
        self.assertEqual(validated, ["blca"])
        self.assertEqual(config["deepmil_network"], "DSMIL")
        self.assertEqual(config["loss_type"], "SurvIFMLE")
        self.assertEqual(config["epochs"], 10)
        self.assertFalse(config["deepmil_use_feat_proj"])
        self.assertTrue(config["save_prediction"])
        self.assertIn("/dsmil/blca-DSMIL-CONCH-fold_0", config["save_path"])
        self.assertEqual(expected["save_path"], "old")

    def make_completion(self, root, score=0.6):
        folder = driver.result_path(root, "blca", 0)
        folder.mkdir(parents=True)
        contents = {"train_metrics-last.txt": str(score),
                    "bag_clf_train_last_pred_test.csv": "patient,prediction\na,1\n",
                    "additional_baseline_protocol.json": '{"model": "DSMIL"}'}
        for name, content in contents.items():
            (folder / name).write_text(content, encoding="utf-8")
        completion = {"wall_seconds": 60, "sha256": {
            name: driver.sha(folder / name) for name in contents}}
        (folder / "dsmil_completed.json").write_text(json.dumps(completion), encoding="utf-8")
        return folder

    def test_completed_hash_and_protocol_checks(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            folder = self.make_completion(root)
            base = SimpleNamespace(read_c_index=lambda p: float(p.read_text()))
            self.assertEqual(driver.read_completed(base, root, "blca", 0, {"model": "DSMIL"}), (0.6, 60))
            with self.assertRaisesRegex(RuntimeError, "protocol differs"):
                driver.read_completed(base, root, "blca", 0, {"model": "ILRA"})
            (folder / "bag_clf_train_last_pred_test.csv").write_text("changed")
            with self.assertRaisesRegex(RuntimeError, "file changed"):
                driver.read_completed(base, root, "blca", 0)

    def test_incomplete_and_invalid_results_are_rejected(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            base = SimpleNamespace(read_c_index=lambda p: float(p.read_text()))
            with self.assertRaisesRegex(RuntimeError, "completion record"):
                driver.read_completed(base, root, "blca", 0)
            self.make_completion(root, float("nan"))
            with self.assertRaisesRegex(ValueError, "Invalid C-index"):
                driver.read_completed(base, root, "blca", 0)

    def test_incomplete_run_is_not_overwritten(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = driver.config_path(root, "blca", 0)
            config.parent.mkdir(parents=True)
            config.write_text("config")
            folder = driver.result_path(root, "blca", 0)
            folder.mkdir(parents=True)
            protected = folder / "existing_result.txt"
            protected.write_text("keep")
            base = SimpleNamespace(raise_file_limit=lambda: None,
                                   split_path=lambda *a: config, table_path=lambda *a: config)
            args = SimpleNamespace(datasets=["blca"], folds=[0], python="python")
            with self.assertRaisesRegex(RuntimeError, "Incomplete output preserved"):
                driver.run(base, args, root, {})
            self.assertEqual(protected.read_text(), "keep")

    def test_partial_summary_is_not_macro_result(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.make_completion(root)
            base = SimpleNamespace(read_c_index=lambda p: float(p.read_text()))
            output = io.StringIO()
            with redirect_stdout(output):
                driver.summarize(base, root, ["blca"], [0])
            self.assertIn("Partial results only", output.getvalue())
            self.assertNotIn("Five-cohort macro C-index:", output.getvalue())

    @unittest.skipUnless((Path(__file__).resolve().parents[1] / "model/deepmil.py").is_file(),
                         "Requires upstream VLSA model code, intentionally not bundled in the analysis release.")
    def test_existing_dsmil_cpu_forward_backward(self):
        import torch
        import torch.nn as nn
        import torch.nn.functional as F

        root = Path(__file__).resolve().parents[1]
        tree = ast.parse((root / "model/deepmil.py").read_text(encoding="utf-8"))
        # Isolate unchanged DSMIL classes from unrelated optional model imports.
        names = {"FCLayer", "BClassifier", "DSMIL"}
        nodes = [n for n in tree.body if isinstance(n, ast.ClassDef) and n.name in names]
        self.assertEqual({n.name for n in nodes}, names)
        namespace = {"torch": torch, "nn": nn, "F": F}
        exec(compile(ast.Module(body=nodes, type_ignores=[]), "repository_dsmil", "exec"), namespace)
        torch.set_num_threads(2)
        for classes in (4, 29):
            net = namespace["DSMIL"](512, 256, classes, use_feat_proj=False)
            for patches in (1, 128, 2048):
                net.zero_grad()
                output = net(torch.randn(1, patches, 512))
                self.assertEqual(tuple(output.shape), (1, classes))
                output.log_softmax(dim=1)[0, 0].neg().backward()
                self.assertTrue(torch.isfinite(output).all().item())
                self.assertTrue(all(torch.isfinite(p.grad).all().item()
                                    for p in net.parameters() if p.grad is not None))


if __name__ == "__main__":
    unittest.main()
