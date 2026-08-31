# Results Draft: Frozen Inductively Calibrated Pathomic Fusion

This section follows the locked experimental record and is intended for direct adaptation into a *Computerized Medical Imaging and Graphics* manuscript. Values from the primary run and values from later common-runtime sensitivity refits are identified separately.

## 1. Primary outer five-fold performance

FICPF achieved the highest mean C-index among the compared methods on all five cohorts. Relative to the matched VLSA baseline, its mean C-index increased from 0.6180 to 0.6334 on BLCA, from 0.6676 to 0.6679 on BRCA, from 0.8005 to 0.8182 on GBMLGG, from 0.6356 to 0.6438 on LUAD, and from 0.7572 to 0.7618 on UCEC. The corresponding macro average increased from 0.6958 to 0.7050, a gain of 0.0092. FICPF also exceeded the published VLSA reference on every cohort, with a macro gain of 0.0096.

| Cohort | Published VLSA | Matched VLSA | MCAT-CONCH | CAT-TPL-CONCH | FICPF | FICPF vs matched |
|---|---:|---:|---:|---:|---:|---:|
| BLCA | 0.6176 | 0.6180 | 0.6214 | 0.6263 | **0.6334** | +0.0154 |
| BRCA | 0.6652 | 0.6676 | 0.5834 | 0.6084 | **0.6679** | +0.0003 |
| GBMLGG | 0.8002 | 0.8005 | 0.7918 | 0.8109 | **0.8182** | +0.0177 |
| LUAD | 0.6370 | 0.6356 | 0.6077 | 0.6067 | **0.6438** | +0.0082 |
| UCEC | 0.7571 | 0.7572 | 0.6432 | 0.6549 | **0.7618** | +0.0046 |
| Macro | 0.6954 | 0.6958 | 0.6495 | 0.6614 | **0.7050** | +0.0092 |

FICPF improved 18 of 25 outer folds. The paired macro difference relative to matched VLSA was 0.0092 (95% bootstrap CI, 0.0016 to 0.0169; (P[\Delta>0]=0.9928)). At cohort level, the confidence interval excluded zero for BLCA (gain 0.0154; 95% CI, 0.0007 to 0.0299) and GBMLGG (gain 0.0177; 95% CI, 0.0096 to 0.0265). LUAD improved in all five folds, but its cohort-level interval crossed zero (gain 0.0082; 95% CI, -0.0069 to 0.0237). BRCA and UCEC showed small mean gains with intervals crossing zero. These results support a consistent macro-level improvement but do not establish a statistically resolved benefit in every individual cohort.

## 2. Comparison with multimodal baselines

Under the unified patient splits and metric implementation, MCAT-CONCH and CAT-TPL-CONCH obtained macro C-indices of 0.6495 and 0.6614, respectively. FICPF improved these values by 0.0555 (95% CI, 0.0278 to 0.0829) and 0.0436 (95% CI, 0.0174 to 0.0710). The largest differences occurred on BRCA and UCEC, where end-to-end multimodal attention was less stable than the pathology baseline. CAT-TPL-CONCH remained competitive on GBMLGG, but FICPF was higher by 0.0073. These comparisons indicate that conservative calibrated fusion was more robust than forcing molecular information through a jointly optimized deep multimodal representation in the evaluated low-event cohorts.

## 3. Component ablation

The molecular-only branch obtained a macro C-index of 0.6273 and was 0.0685 below pathology-only prediction. Directly combining the uncalibrated pathology and molecular scores with the same 0.80/0.20 weights yielded only a 0.0014 macro gain over pathology (95% CI, -0.0003 to 0.0035). After training-fold ECDF alignment, locked FICPF improved the uncalibrated fusion by 0.0078 (95% CI, 0.0011 to 0.0152), demonstrating that score-distribution alignment accounted for a substantial part of the final improvement.

Equally weighted ECDF fusion achieved strong results on BLCA and GBMLGG but reduced performance on BRCA and UCEC. Its macro gain over pathology was 0.0031 (95% CI, -0.0153 to 0.0220), compared with 0.0092 for the conservative 0.80/0.20 configuration. Thus, the molecular branch was useful as complementary evidence even though it was not a strong standalone predictor, and limiting its weight reduced cross-cohort negative transfer.

| Variant | Macro C-index | Difference from pathology |
|---|---:|---:|
| Pathology only | 0.6958 | 0 |
| Molecular only | 0.6273 | -0.0685 |
| Raw fixed fusion | 0.6972 | +0.0014 |
| Equal ECDF fusion | 0.6988 | +0.0031 |
| FICPF | **0.7050** | **+0.0092** |

## 4. Patient-specific molecular correspondence

Randomly permuting the molecular predictions among test patients reduced C-index relative to the true FICPF pairing by 0.0211 on BLCA, 0.0286 on GBMLGG, 0.0123 on LUAD, and 0.0175 on UCEC. The probabilities that a random permutation equaled or exceeded the correctly paired result were 0.001, <0.001, 0.015, and 0.031, respectively. The BRCA permutation probability was 0.164, consistent with the negligible BRCA gain in the primary comparison. The negative control therefore supports the interpretation that the improvement was driven by patient-specific molecular complementarity rather than by adding a generic smoothing distribution to pathology scores.

## 5. Fixed-parameter sensitivity

Across fixed molecular weights, macro C-index was 0.6957, 0.7007, 0.7050, 0.7059, and 0.6988 for (\alpha=0,0.1,0.2,0.3,) and 0.5, respectively. The interval from 0.10 to 0.30 formed a stable plateau. Although (\alpha=0.30) was 0.0008 higher than the locked primary weight, this difference was observed after opening the outer tests and was not used to redefine the method. Equal weighting caused marked decreases on BRCA and UCEC.

In a common numerical environment, the requested PCA dimensions 16, 32, 64, and 128 produced macro C-indices of 0.6998, 0.7004, 0.7053, and 0.6984. Cohort-specific optima differed, whereas PCA64 provided the strongest common setting. The common-runtime PCA64 value differed from the locked primary macro result by only 0.0003 and was retained as a sensitivity analysis rather than substituted into the primary table.

## 6. Survival-curve quality

The monotone Cox-Breslow projection preserved scalar FICPF concordance exactly. Across cohorts, macro IBS decreased from 0.2296 for matched VLSA to 0.2001 for FICPF, an improvement of 0.0295, and IBS improved on all five cohorts. The number of folds passing D-calibration increased from 14/25 to 25/25. Mean Kaplan-Meier calibration error decreased from 0.0283 to 0.0178.

In contrast, hinge MAE increased from 25.23 to 29.89 months, a deterioration of 4.66 months, and no cohort showed an MAE improvement. FICPF should therefore be interpreted primarily as a ranking and distribution-calibration method. The current curve projection improved population-level probability calibration but did not improve individual absolute survival-time estimation.

## 7. Robustness to molecular corruption and missingness

The common-runtime clean FICPF refit obtained a macro C-index of 0.7053 compared with 0.6957 for pathology. With 50% random feature dropout, Gaussian noise of standard deviation 0.50 in standardized space, and 50% patient-level molecular missingness, macro C-index remained 0.6983, 0.7019, and 0.6987, respectively. These values retained gains of 0.0026, 0.0063, and 0.0031 over pathology. At 100% patient-level missingness, the explicit fallback reproduced pathology exactly.

The robustness effect was not uniform across cohorts. BLCA, GBMLGG, and LUAD retained positive gains under substantial feature corruption. BRCA fluctuated around its already small clean gain. UCEC, which had only 123 signature-linked variables, fell below pathology by 0.0141 with 50% feature dropout and by 0.0116 with noise standard deviation 0.50. Feature-level robustness is therefore contingent on sufficient molecular coverage, whereas patient-level fallback remains deterministic and safe.

## 8. Numerical stability and computational cost

Changing only the randomized PCA seed among 0, 1, 42, 3404, and 3405 produced macro C-indices between 0.7039 and 0.7053. The across-seed standard deviation was 0.0006 and the total range was 0.0014, substantially smaller than the primary gain over matched VLSA. This experiment characterizes the numerical stability of the frozen PCA-Cox branch and does not replace full multi-seed retraining of the VLSA pathology model.

FICPF added no trainable deep parameters and only 65 Cox parameters per fold. Across all 25 folds, scaler/PCA/Cox/ECDF fitting required 21.1 CPU seconds in total, with a mean of 0.85 seconds per fold. The mean persistent state was 1.26 MB in float64, and incremental CPU inference averaged 0.0149 ms per patient. These values exclude the pre-existing VLSA computation and molecular-data acquisition.

## 9. Overall interpretation

The principal result is a reproducible macro C-index increase with a paired confidence interval excluding zero, achieved without modifying the VLSA backbone or selecting test-specific parameters. The strongest cohort-level evidence was observed on BLCA and GBMLGG, while LUAD showed consistent but statistically uncertain gains and BRCA contributed almost no improvement. The ablations identify inductive score calibration and conservative weighting as necessary components. The method improves concordance, IBS, and calibration, but its benefit does not extend to absolute-time MAE and it remains sensitive to severe feature-level molecular degradation in low-coverage cohorts.
