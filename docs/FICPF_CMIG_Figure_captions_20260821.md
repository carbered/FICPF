# Figure captions

**Graphical abstract. Frozen Inductively Calibrated Pathomic Fusion.** Frozen VLSA pathology predictions and a compact RNA/CNV PCA-Cox branch are aligned by training-fold empirical cumulative distribution functions and combined using a fixed 0.80/0.20 weight. The resulting calibrated survival score reaches a macro C-index of 0.705 across five TCGA cohorts, a gain of 0.009 over matched pathology prediction.

**Figure 1. Overview of FICPF.** A matched VLSA-CONCH model provides frozen pathology survival curves and mean predicted survival times. Signature-linked RNA-seq/CNV features are processed by training-fold standardization, PCA64, and an L2-regularized Cox model. The two time-oriented scores are independently mapped through training-fold empirical cumulative distribution functions and combined using a fixed pathology/molecular weight of 0.80/0.20. Missing molecular inputs trigger pathology fallback. Complete survival curves are reconstructed by a monotone training-fold Cox-Breslow projection.

**Figure 2. Unified outer five-fold comparison across five cohorts.** The left panel compares published VLSA, matched VLSA, MCAT-CONCH, CAT-TPL-CONCH, and FICPF. The right panel shows paired patient-level bootstrap differences between FICPF and matched VLSA with 95% confidence intervals.

**Figure 3. Component ablation and fixed-parameter sensitivity.** Macro C-index is shown for principal method components, molecular fusion weights, and requested PCA dimensions. Diamonds indicate the locked molecular weight of 0.20 and PCA64 configuration.

**Figure 4. Multi-metric comparison of survival-curve quality.** Matched VLSA and FICPF are compared across five cohorts using the integrated Brier score (IBS), Kaplan-Meier calibration error, hinge mean absolute error (MAE), and the number of outer folds passing D-calibration. Arrows indicate the preferred direction. FICPF improves probability error and population-level calibration but worsens absolute-time MAE.

**Figure 5. Test-time robustness to molecular corruption.** C-index gains over pathology are shown under random feature dropout, Gaussian noise, and patient-level molecular missingness. Black lines denote the unweighted macro average across the five cohorts.

**Figure 6. Fixed-threshold out-of-fold risk stratification.** (A-E) Kaplan-Meier curves for high- and low-risk patients in the five TCGA cohorts. Each patient is scored by the outer-fold model that did not use that patient for fitting. High risk is defined prospectively as a calibrated fused percentile below 0.50; no survival outcome is used to optimize the cutoff. Censoring marks are shown on the curves. (F) Unadjusted and clinically adjusted Cox hazard ratios with 95% confidence intervals. Adjustment uses age, identifiable sex, and sufficiently represented OncoTree histologic subtype; adjusted estimates use complete cases.
