# FICPF

Frozen Inductively Calibrated Pathomic Fusion (FICPF) is a leakage-conscious
late-fusion framework for pathology-omics survival analysis. It preserves a
trained VLSA-CONCH pathology predictor, fits a compact RNA-seq/CNV PCA-Cox
branch inside each outer training fold, maps both scores to training-fold
percentiles, and combines them with a fixed pathology/omics weight of 0.80/0.20.

![FICPF method overview](figures/ficpf_20260821/ficpf_method_overview.png)

## Locked results

The frozen five-cohort outer five-fold analysis reports:

| Measure | Value |
|---|---:|
| Matched VLSA macro C-index | 0.695782 |
| FICPF macro C-index | 0.705017 |
| Paired macro improvement | +0.009235 |
| 95% bootstrap CI | [0.0016, 0.0169] |
| Cohorts | BLCA, BRCA, GBMLGG, LUAD, UCEC |
| Pathology training seeds | 42, 3404, 3405 |
| Positive seed-by-cohort point estimates | 15/15 |

The main result remains the locked seed-42 analysis. The additional pathology
training seeds are a stability analysis and do not replace the main table.

## Repository contents

- `tools/`: fitting, evaluation, ablation, robustness, statistics, and figure scripts.
- `data_split/`: fixed patient-level outer five-fold splits and survival tables.
- `result/`: locked patient-level predictions and derived statistical summaries.
- `figures/`: publication figures in PNG, PDF, and SVG formats.
- `docs/`: English/Chinese manuscript drafts, evidence notes, captions, and references.
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
stability summary, and the GitHub package size limit.

## Reproduce analyses

Summary analyses based only on stored patient-level predictions can be rerun
directly. For example:

```bash
python tools/summarize_frozen_pathomic_fusion_five_datasets.py
```

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
