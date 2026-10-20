# FICPF

Frozen Inductively Calibrated Pathomic Fusion (FICPF) is a leakage-conscious
late-fusion framework for pathology-omics survival analysis. It preserves a
trained VLSA pathology predictor using CONCH patch features, fits a compact RNA-seq/CNV PCA-Cox
branch inside each outer training fold, maps both scores to training-fold
percentiles, and combines them with a fixed pathology/omics weight of 0.80/0.20.

![FICPF method overview](figures/ficpf_model_names/FICPF_method_overview.png)

## Locked results

The frozen five-cohort outer five-fold analysis reports:

| Measure | Value |
|---|---:|
| VLSA macro C-index | 0.695782 |
| FICPF macro C-index | 0.705017 |
| Paired macro improvement | +0.009235 |
| 95% bootstrap CI | [0.0016, 0.0169] |
| Cohorts | BLCA, BRCA, GBMLGG, LUAD, UCEC |
| Pathology training seeds | 42, 3404, 3405 |
| Positive seed-by-cohort point estimates | 15/15 |

The main result remains the locked seed-42 analysis. The additional pathology
training seeds are a stability analysis and do not replace the main table.

## Current comparisons

The eight-method analysis includes five pathology methods, two multimodal
baselines, and FICPF. All use precomputed CONCH pathology features. Display
labels omit the encoder suffix; raw result identifiers are unchanged.

| Method | Five-cohort macro C-index |
|---|---:|
| ABMIL | 0.6215 |
| TransMIL | 0.6447 |
| ILRA | 0.6364 |
| DSMIL | 0.6591 |
| VLSA | 0.6958 |
| MCAT* | 0.6495 |
| CAT-with-TPL* | 0.6614 |
| FICPF | 0.7050 |

*The retained MCAT/CAT-with-TPL adaptations use cohort-wide event-time bins.
Their binning protocol differs from the training-fold-only bins of the
pathology baselines. See the retained configs and audit records.

Per-cohort values are in
[`eight_method_dataset_summary.csv`](result/ficpf-eight-method-comparison/eight_method_dataset_summary.csv).
The [current comparison figure](figures/ficpf_model_names/FICPF_eight_method_comparison.pdf)
uses the same locked results.

Six methods have retained survival curves. Their full metrics are in
[`macro_curve_metrics.csv`](result/ficpf-six-method-curve-comparison/macro_curve_metrics.csv).
FICPF has lower macro IBS (0.2001 versus 0.2296) and KM calibration error
(0.0178 versus 0.0283) than VLSA, but higher hinge MAE (29.89 versus 25.23
months). D-calibration counts are descriptive counts over 25 held-out folds,
not pooled patient-level calibration tests. The scalar-risk outputs retained
for MCAT and CAT-with-TPL are not used to invent survival curves.

## Dataset and resource audits

The locked 2,831 patients contribute 3,529 available slides and 15,615,143
patches. See [`dataset_statistics.csv`](result/ficpf-dataset-resource-audit/dataset_statistics.csv)
for per-cohort slide, patch, omics-variable and time-bin statistics. Missing
slide feature files are excluded consistently with the original loader.

The [resource benchmark](result/ficpf-dataset-resource-audit/pathology_benchmark.json)
uses a representative BLCA fold-0 input (6,956 patches, 512 features per patch,
3,062 omics variables), cached patient-independent text features, RTX 5070 Ti
Laptop GPU, and single-thread i9-14900HX CPU processing. The retained pathology
checkpoint reproduces its saved curve within 2.4e-8 maximum absolute error.
The 100-repeat median feature-to-curve-and-mean-time latencies are 25.61 ms for
VLSA and 62.17 ms for FICPF. These include CPU postprocessing and numerical
mean-time integration, not WSI encoding, file I/O or input transfer, and are
not cohort-average timings. FICPF is not claimed to be faster at inference.
Dense-linear MACs are 85.744 M versus 85.940 M; FICPF adds 65 Cox coefficients
including the intercept and no deep-network training parameters. Data-fitted
PCA/scaler/calibration state is not counted as Cox coefficients.

## Repository contents

- `tools/`: fitting, evaluation, ablation, robustness, statistics, and figure scripts.
- `data_split/`: fixed patient-level outer five-fold splits and survival tables.
- `result/`: locked patient-level predictions and derived statistical summaries.
- `figures/`: publication figures in PNG, PDF, and SVG formats.
- `docs/`: historical drafts, evidence notes, captions, and references. The current
  Chinese Word manuscript is private and is not included in this update.
- `eval/`, `dataset/`, `utils/`: evaluation and data utilities needed by the scripts.
- `manifest/`: environment metadata and a per-file SHA256 manifest.

Raw WSIs, CONCH weights, WSI patch features, VLSA checkpoints, and large
preprocessed omics tables are deliberately excluded. See [DATA.md](DATA.md).

## Installation

Python 3.10 or later is recommended.

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
```

## Verify the frozen package

This check does not require the omitted omics tables:

```bash
python tools/verify_ficpf_release.py
```

It verifies the SHA256 manifest, package scope, headline results, three-seed
stability summary, current result aggregations, and retained prediction hashes.
The public-file budget is 50 MiB, excluding Git internals; individual files
remain below 25 MiB.

## Reproduce analyses

Summary analyses based only on stored patient-level predictions can be rerun
directly. For example:

```bash
python tools/summarize_frozen_pathomic_fusion_five_datasets.py
```

Check the four retained MIL methods against the original continuous-label
C-index evaluator without retraining:

```bash
python tools/verify_current_results.py --recompute-c-index
```

Recompute the six-method curve comparison without the large omics tables:

```bash
python tools/evaluate_additional_baseline_curves.py \
  --output-dir result/reproduced_curve_comparison
```

The 100 retained test-curve CSVs, model configs, exact archived label converter,
fixed patient splits, and FICPF time grids are included. Full recomputation is
slower than summary verification because it integrates individual curves.

Regenerate the current eight-method figure:

```bash
python tools/plot_ficpf_six_method_comparison.py --include-dsmil \
  --canonical-names --output-stem figures/reproduced_eight_method_comparison
```

The method schematic can be regenerated with
`tools/refine_ficpf_main_figure.py --tinted --canonical-names --output-stem figures/reproduced_method`; its small public
TCGA example images and provenance are included. The optional local layout-audit
tool is not required to export the figures. Source audits and resource timing
that consume model weights or patch features still require the omitted upstream
inputs; see [DATA.md](DATA.md). No model is retrained by the verification commands.
The optional archive-audit utilities require the original cloud export archives
and a Python 3.12+ runtime; they are not the data-free verification entry points.

The full ablation, external-baseline comparison, PCA-Cox re-fitting, and omics
corruption experiments require the omitted preprocessed tables. Place them as
described in [DATA.md](DATA.md), then run, for example:

```bash
python tools/evaluate_ficpf_pca_sensitivity.py \
  --omics-table-dir external/CAT-with-TPL-master/dataset_csv \
  --output-dir result/reproduced_pca
```

Re-training the VLSA pathology backbone is an upstream workflow and additionally
requires the original VLSA code, CONCH model weights, and WSI patch features.

## Citation and license

The manuscript citation will be added after publication. Methodological sources
are listed in `docs/references.bib`. The repository preserves the license shipped
with the source workspace; see `LICENSE` and [NOTICE.md](NOTICE.md).

A Chinese quick-start is available in [README_中文.md](README_中文.md).
