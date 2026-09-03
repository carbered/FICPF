# Matched ABMIL and TransMIL baselines

This experiment uses the ABMIL and TransMIL implementations already included
in VLSA. It consumes precomputed `N x 512` CONCH feature bags and does not need
raw WSIs or additional pretrained weights.

Run a one-fold smoke test first:

```bash
cd /tmp/VLSA_matched_multiseed_20260828
PY=/root/miniconda3/envs/vlsa_env/bin/python

$PY tools/run_matched_mil_backbones.py prepare --verify-tensors

WANDB_MODE=offline $PY tools/run_matched_mil_backbones.py run \
  --models ABMIL TransMIL --datasets brca --folds 0
```

Run all 50 tasks:

```bash
mkdir -p logs
LOG=logs/matched_mil_backbones_$(date +%Y%m%d_%H%M%S).log

WANDB_MODE=offline nohup $PY -u tools/run_matched_mil_backbones.py run \
  --models ABMIL TransMIL \
  --datasets blca brca gbmlgg luad ucec \
  --folds 0 1 2 3 4 \
  > "$LOG" 2>&1 &

echo "PID=$!"
echo "LOG=$LOG"
```

Summarize completed folds:

```bash
$PY tools/run_matched_mil_backbones.py summarize
```

The runner skips tasks with an existing `train_metrics-last.txt`, so the full
command is safe to resume after interruption.
