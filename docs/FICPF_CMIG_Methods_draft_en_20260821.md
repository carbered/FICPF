# Methods Draft: Frozen Inductively Calibrated Pathomic Fusion

This text is written as a manuscript-ready draft for *Computerized Medical Imaging and Graphics*. Citation labels such as `[VLSA]`, `[CONCH]`, and `[Cox]` should be replaced by the final bibliography keys during manuscript assembly.

## 1. Problem formulation and evaluation protocol

We considered right-censored survival prediction from paired whole-slide pathology and molecular measurements. For patient (i), the available data were a bag of precomputed pathology patch embeddings (X_i^p), a molecular feature vector (x_i^o), an observed follow-up time (t_i), and an event indicator (e_i\in\{0,1\}). Experiments used five TCGA cohorts: BLCA, BRCA, GBMLGG, LUAD, and UCEC.

All model fitting and calibration were performed independently inside each outer training fold. The corresponding outer test fold was used once for evaluation and did not contribute labels, predictions, distribution statistics, feature preprocessing, fusion weights, checkpoint selection, or hyperparameter selection. The fixed patient-level five-fold partitions yielded 373/169, 956/130, 569/187, 453/158, and 480/75 out-of-fold patients/events for BLCA, BRCA, GBMLGG, LUAD, and UCEC, respectively.

## 2. Frozen pathology survival branch

The pathology branch was the matched VLSA-CONCH reproduction `[VLSA, CONCH]`. Precomputed 512-dimensional CONCH patch embeddings were aggregated by the VLFAN vision-language attention module using ordinal survival prompts. For each patient, the model produced a discrete survival curve

\[
S_i^p(t_j)=P(T_i>t_j\mid X_i^p),
\]

on the cohort-specific time grid. A scalar time-oriented score (T_i^p) was obtained as the mean predicted survival time under linear interpolation of (S_i^p). Larger (T_i^p) therefore indicated longer predicted survival.

FICPF did not modify or jointly retrain the pathology network. The matched raw last checkpoint and its saved outer-training and outer-test survival predictions were treated as fixed inputs. This design isolates the contribution of inductive multimodal calibration from changes in pathology optimization, early stopping, or checkpoint selection.

## 3. Compact molecular survival branch

The molecular input comprised RNA-seq expression and copy-number variation features linked to six curated gene groups: tumor suppressor genes, oncogenes, protein kinases, cell differentiation markers, transcription factors, and cytokines and growth factors. Only columns with the suffix `_rnaseq` or `_cnv` and a gene listed in one of these groups were retained. Depending on cohort coverage, this produced 3062, 3189, 2832, 3096, and 123 molecular variables for BLCA, BRCA, GBMLGG, LUAD, and UCEC, respectively.

Within outer fold (k), a StandardScaler was fitted to the molecular features of the training patients only. Principal component analysis with randomized SVD then reduced the standardized variables to

\[
d_k=\min(64,n_k^{train}-1,p_k)
\]

components, where (p_k) is the number of available molecular variables. The randomized PCA seed was fixed to 42 for the primary analysis.

An L2-regularized Cox proportional-hazards readout was fitted to the training-fold components (z_i\in\mathbb{R}^{d_k}). With linear predictor (\eta_i=w^\top z_i+b), the optimized objective was

\[
\mathcal{L}_{Cox}
=-\frac{1}{\sum_i e_i}
\sum_{i:e_i=1}
\left[\eta_i-\log\sum_{j:t_j\geq t_i}\exp(\eta_j)\right]
+\lambda\frac{\lVert w\rVert_2^2}{d_k}.
\]

We used (\lambda=0.10), Adam with learning rate 0.02, and 1200 optimization steps. The molecular time-oriented score was defined as (T_i^o=-\eta_i), so that both branches used the same direction: a larger score represented longer predicted survival.

## 4. Inductive cross-modal calibration

Directly averaging (T_i^p) and (T_i^o) is ill-posed because the two branches have unrelated score ranges and marginal distributions. FICPF aligns each modality to its own outer-training distribution. For modality (m\in\{p,o\}), the empirical cumulative distribution function was

\[
\widehat F_m^{train}(z)
=\frac{1}{n_k^{train}}\sum_{j\in\mathcal D_k^{train}}
\mathbb{1}(T_j^m\leq z).
\]

The calibrated patient scores were

\[
q_i^p=\widehat F_p^{train}(T_i^p),\qquad
q_i^o=\widehat F_o^{train}(T_i^o).
\]

Both empirical CDFs were estimated without held-out labels. The transformation is monotone and introduces no pairwise inversions within a modality, although finite training sets can map nearby test scores to the same percentile and therefore create ties.

## 5. Conservative pathomic fusion and missing-modality fallback

The final FICPF score was a fixed convex combination

\[
q_i=(1-\alpha)q_i^p+\alpha q_i^o,
\qquad \alpha=0.20.
\]

The same value was used for every cohort and outer fold. It was locked before final outer-test evaluation and was not adapted to cohort identity, validation performance, or test performance. The deliberately conservative molecular weight limits negative transfer when the molecular branch is weak while retaining patient-specific complementary information.

When molecular measurements were unavailable for a patient, inference reverted explicitly to the calibrated pathology score, (q_i=q_i^p). Thus, missing molecular data did not require imputation by a learned multimodal network and could not disable pathology-only inference.

## 6. Monotone survival-curve reconstruction

The scalar FICPF score directly defines the patient ranking used for C-index evaluation. To obtain complete survival curves for calibration and Brier-score analyses, we fitted a monotone Cox-Breslow projection using outer-training patients only. First, percentiles were clipped to ([0.5/n_k^{train},1-0.5/n_k^{train}]) and converted to a risk-oriented feature

\[
u_i=-\log\frac{q_i}{1-q_i}.
\]

A positive scalar (\gamma\in[0.05,5]) was estimated by minimizing the Cox partial likelihood on the outer training fold. A Breslow baseline survival function (S_0(t)) was then estimated and the patient curve was defined as

\[
S_i^{FICPF}(t)=S_0(t)^{\exp(\gamma u_i)}.
\]

The first and final points were set to one and zero, respectively, to match the finite discrete-time support of the VLSA output grid. Because (\gamma>0) and the logit mapping is monotone, this projection preserved the FICPF ordering exactly; the maximum observed C-index difference between scalar and curve evaluation was zero.

## 7. Comparators and evaluation metrics

The primary comparator was the matched VLSA-CONCH reproduction evaluated on the same patients, folds, and raw last-checkpoint policy. Published VLSA values were included as unpaired references. MCAT-CONCH and CAT-TPL-CONCH were reproduced as multimodal external comparators using the same outer splits, CONCH pathology features, and available molecular tables. Their fixed protocol used seed 1, 20 epochs, learning rate (10^{-4}), gradient accumulation of 16 patients, weighted sampling, no early stopping, and the final checkpoint.

The primary endpoint was Harrell-style concordance under the `SurvivalEVAL_mean_continuous_comparable_all_ties` contract: predicted mean survival time, comparable pairs, and all ties assigned 0.5. Secondary curve metrics were the inverse-probability-of-censoring-weighted integrated Brier score (IBS), hinge mean absolute error (MAE), D-calibration, and Kaplan-Meier calibration error.

For paired uncertainty estimates, patients were resampled with replacement within each outer fold. A fold-level C-index difference was calculated for each resample and averaged across the five folds. We report the 2.5th and 97.5th percentiles of 5000 bootstrap replicates and the empirical probability that the difference was positive. Published VLSA values were not used in paired bootstrap tests because patient-level predictions were unavailable.

## 8. Prespecified ablation, sensitivity, and robustness analyses

The component analysis compared pathology-only prediction, molecular-only prediction, uncalibrated fixed fusion, equally weighted ECDF fusion, and locked FICPF. A patient correspondence control randomly permuted test-patient molecular scores 1000 times without consulting survival labels. Fusion-weight sensitivity used (\alpha\in\{0,0.1,0.2,0.3,0.5\}), and PCA sensitivity used 16, 32, 64, and 128 requested components. These analyses were descriptive and did not replace the locked (\alpha=0.20), PCA64 primary configuration after outer-test observation.

Test-time molecular robustness was evaluated with random feature dropout and training-mean imputation, Gaussian noise in the training-standardized feature space, and patient-level molecular missingness with pathology fallback. Each stochastic condition was repeated 20 times, and the molecular model and all training-fold statistics remained fixed.

## 9. Incremental computational cost

FICPF added no trainable deep-network parameters. With PCA64, the Cox readout contained 64 weights and one bias per fold. The minimal persistent state comprised the scaler statistics, PCA mean and projection, Cox coefficients, and sorted outer-training scores required by the two empirical CDFs. Reported incremental timing excluded the pre-existing VLSA forward pass and the acquisition of pathology and molecular data.
