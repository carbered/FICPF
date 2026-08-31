# Data setup

## Included data

This repository includes fixed patient-level survival tables and outer five-fold
splits under `data_split/`, the CAT-with-TPL gene-signature definition, and the
locked patient-level predictions required to audit the reported results.

## Omitted large inputs

The preprocessed RNA-seq/CNV tables are omitted to keep the GitHub package below
25 MB and to avoid silently redistributing inherited TCGA derivatives. Obtain the
data under the applicable NCI Genomic Data Commons policy and reproduce the
CAT-with-TPL preprocessing. Place the complete tables here:

```text
external/CAT-with-TPL-master/dataset_csv/
  signatures.csv
  tcga_blca_all_clean.csv.zip
  tcga_brca_all_clean.csv.zip
  tcga_gbmlgg_all_clean.csv.zip
  tcga_luad_all_clean.csv.zip
  tcga_ucec_all_clean.csv.zip
```

Expected SHA256 values for the exact complete tables used in the locked analysis:

| File | SHA256 |
|---|---|
| `tcga_blca_all_clean.csv.zip` | `386e0adf343eb8b947385d175133a53b75e90a152ae2e0ceb1c180f82c8037e8` |
| `tcga_brca_all_clean.csv.zip` | `ebd2d31abb871ae1aec59606a9c4a5f7280fe238ce440ab376c70cac8bf6e974` |
| `tcga_gbmlgg_all_clean.csv.zip` | `6b4fc7a938dbf49bcf536f17e24eb52a8fc0721ac9243ab6fd671194d157d345` |
| `tcga_luad_all_clean.csv.zip` | `e90a1a28c79cab829737c56c2306cd7912bbf7b4485828a507d6b82033a174d1` |
| `tcga_ucec_all_clean.csv.zip` | `215bf025477bcf8ff6c7865f8c7a9a1a26ee0619eb36aba64618bd6f3e9d8e1b` |

Some unified external-baseline and ablation scripts also consume the
signature-filtered CAT-with-TPL tables under `dataset_csv_sig/`. Their hashes are
recorded in `manifest/omitted_data_sha256.csv`.

## Upstream pathology inputs

Raw WSIs, CONCH weights, WSI patch embeddings, and VLSA checkpoints are not
included. They are not needed to verify the frozen patient-level analysis, but
they are required to regenerate pathology predictions from the beginning.

## Integrity check

After adding any omitted table, verify it with:

```bash
sha256sum external/CAT-with-TPL-master/dataset_csv/*.zip
```

On Windows PowerShell, use `Get-FileHash <path> -Algorithm SHA256`.
