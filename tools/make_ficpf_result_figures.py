#!/usr/bin/env python3
"""Create publication-ready FICPF result, ablation, and robustness figures."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "ficpf_20260821"
COLORS = {
    "paper": "#202124",
    "matched": "#8A8F98",
    "mcat": "#0072B2",
    "cat": "#009E73",
    "ficpf": "#D55E00",
    "omics": "#CC79A7",
    "raw": "#56B4E9",
    "equal": "#E69F00",
}
DATASET_COLORS = {
    "BLCA": "#0072B2",
    "BRCA": "#7A7F87",
    "GBMLGG": "#009E73",
    "LUAD": "#E69F00",
    "UCEC": "#D55E00",
    "MACRO": "#202124",
}


def setup_style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 7.5,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.3,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def style_axis(ax, grid_axis="y"):
    ax.grid(axis=grid_axis, color="#DADCE0", linewidth=0.55, alpha=0.9)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def panel_label(ax, label):
    ax.text(
        -0.09,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        fontweight="bold",
    )


def save_all(fig, stem):
    OUT.mkdir(parents=True, exist_ok=True)
    paths = []
    for extension in ("png", "pdf", "svg"):
        path = OUT / f"{stem}.{extension}"
        kwargs = {"bbox_inches": "tight", "pad_inches": 0.06}
        if extension == "png":
            kwargs["dpi"] = 450
        fig.savefig(path, **kwargs)
        paths.append(str(path))
    plt.close(fig)
    return paths


def make_main_figure():
    summary_path = (
        ROOT
        / "result"
        / "external_multimodal_comparison_20260821"
        / "unified"
        / "dataset_summary_unified.csv"
    )
    bootstrap_path = summary_path.parent / "paired_bootstrap_unified.csv"
    summary = pd.read_csv(summary_path).set_index("dataset")
    order = ["BLCA", "BRCA", "GBMLGG", "LUAD", "UCEC", "MACRO"]
    summary = summary.loc[order]
    methods = [
        ("paper_vlsa", "Published VLSA", COLORS["paper"], "x"),
        ("matched_vlsa", "Matched VLSA", COLORS["matched"], "o"),
        ("mcat_conch", "MCAT-CONCH", COLORS["mcat"], "s"),
        ("cat_tpl_conch", "CAT-TPL-CONCH", COLORS["cat"], "^"),
        ("ficpf", "FICPF", COLORS["ficpf"], "D"),
    ]

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(10.4, 3.45),
        gridspec_kw={"width_ratios": [1.65, 1.0]},
    )
    ax = axes[0]
    x = np.arange(len(order), dtype=float)
    offsets = np.linspace(-0.24, 0.24, len(methods))
    for offset, (column, label, color, marker) in zip(offsets, methods):
        values = summary[column].to_numpy(dtype=float)
        scatter_kwargs = {
            "x": x + offset,
            "y": values,
            "label": label,
            "color": color,
            "marker": marker,
            "s": 34 if column != "ficpf" else 42,
            "linewidth": 0.9,
            "zorder": 3,
        }
        if marker != "x":
            scatter_kwargs["edgecolor"] = color
        ax.scatter(**scatter_kwargs)
    ax.axvline(4.5, color="#B8BABC", linewidth=0.7, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(["Macro" if item == "MACRO" else item for item in order])
    ax.set_ylabel("C-index")
    ax.set_ylim(0.55, 0.84)
    ax.set_title("Unified outer 5-fold comparison")
    style_axis(ax, "y")
    ax.legend(
        frameon=False,
        ncol=5,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.12),
        columnspacing=1.1,
        handletextpad=0.4,
    )
    panel_label(ax, "A")

    ax = axes[1]
    boot = pd.read_csv(bootstrap_path)
    boot = boot[boot["comparison"] == "ficpf_vs_matched"].set_index("dataset").loc[order]
    y = np.arange(len(order), dtype=float)
    delta = boot["observed_delta"].to_numpy(dtype=float)
    lower = boot["delta_ci025"].to_numpy(dtype=float)
    upper = boot["delta_ci975"].to_numpy(dtype=float)
    significant = (lower > 0.0) | (upper < 0.0)
    for index, (yi, value, lo, hi, is_significant) in enumerate(
        zip(y, delta, lower, upper, significant)
    ):
        color = COLORS["ficpf"] if index < len(order) - 1 else COLORS["paper"]
        ax.errorbar(
            value,
            yi,
            xerr=np.array([[value - lo], [hi - value]]),
            fmt="o",
            markersize=5.4,
            markerfacecolor=color if is_significant else "white",
            markeredgecolor=color,
            markeredgewidth=1.1,
            ecolor=color,
            elinewidth=1.2,
            capsize=2.4,
            zorder=3,
        )
        ax.text(hi + 0.0013, yi, f"{value:+.3f}", va="center", fontsize=7.5)
    ax.axvline(0.0, color="#4A4D50", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(["Macro" if item == "MACRO" else item for item in order])
    ax.invert_yaxis()
    ax.set_xlabel("C-index gain vs matched VLSA (95% CI)")
    ax.set_title("Paired bootstrap effect")
    ax.set_xlim(-0.025, 0.038)
    style_axis(ax, "x")
    panel_label(ax, "B")

    fig.tight_layout(w_pad=2.2)
    fig.subplots_adjust(top=0.78)
    return save_all(fig, "ficpf_main_comparison")


def make_ablation_figure():
    ablation_dir = ROOT / "result" / "ficpf-ablation-sensitivity-20260821"
    component = pd.read_csv(ablation_dir / "dataset_summary_ablation.csv")
    alpha = pd.read_csv(ablation_dir / "dataset_summary_alpha_sensitivity.csv")
    pca = pd.read_csv(ROOT / "result" / "ficpf-pca-sensitivity-20260821" / "dataset_summary.csv")

    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.15))
    ax = axes[0]
    component_order = [
        ("pathology_only", "Pathology", COLORS["matched"]),
        ("omics_only", "Omics", COLORS["omics"]),
        ("raw_fixed_fusion", "Raw fusion", COLORS["raw"]),
        ("ecdf_equal_fusion", "ECDF equal", COLORS["equal"]),
        ("ficpf_locked", "FICPF", COLORS["ficpf"]),
    ]
    macro = component[component["dataset"] == "MACRO"].set_index("method")
    y = np.arange(len(component_order), dtype=float)
    values = np.array([macro.loc[key, "mean_c_index"] for key, _, _ in component_order])
    colors = [color for _, _, color in component_order]
    labels = [label for _, label, _ in component_order]
    ax.hlines(y, 0.61, values, color=colors, linewidth=2.0)
    ax.scatter(values, y, color=colors, s=42, zorder=3, edgecolor="#333333", linewidth=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlim(0.61, 0.712)
    ax.set_xlabel("Macro C-index")
    ax.set_title("Component ablation")
    for yi, value in zip(y, values):
        ax.text(value + 0.0015, yi, f"{value:.3f}", va="center", fontsize=7.4)
    style_axis(ax, "x")
    panel_label(ax, "A")

    ax = axes[1]
    macro_alpha = alpha[alpha["dataset"] == "MACRO"].sort_values("alpha")
    ax.plot(
        macro_alpha["alpha"],
        macro_alpha["mean_c_index"],
        color=COLORS["ficpf"],
        marker="o",
        markersize=4.8,
    )
    locked = macro_alpha[np.isclose(macro_alpha["alpha"], 0.20)].iloc[0]
    ax.scatter(
        [locked["alpha"]],
        [locked["mean_c_index"]],
        marker="D",
        s=52,
        color=COLORS["paper"],
        zorder=4,
        label="Locked",
    )
    ax.axvspan(0.10, 0.30, color="#E8F1F7", alpha=0.8, zorder=0)
    ax.set_xticks(macro_alpha["alpha"])
    ax.set_xlabel("Omics fusion weight")
    ax.set_ylabel("Macro C-index")
    ax.set_ylim(0.694, 0.708)
    ax.set_title("Fixed-weight sensitivity")
    style_axis(ax, "y")
    ax.legend(frameon=False, loc="lower right")
    panel_label(ax, "B")

    ax = axes[2]
    macro_pca = pca[pca["dataset"] == "MACRO"].sort_values("requested_pca_components")
    positions = np.arange(len(macro_pca), dtype=float)
    ax.plot(
        positions,
        macro_pca["mean_fused_c_index"],
        color=COLORS["mcat"],
        marker="o",
        markersize=4.8,
    )
    locked_index = int(np.flatnonzero(macro_pca["requested_pca_components"].to_numpy() == 64)[0])
    ax.scatter(
        [locked_index],
        [macro_pca.iloc[locked_index]["mean_fused_c_index"]],
        marker="D",
        s=52,
        color=COLORS["paper"],
        zorder=4,
        label="Locked",
    )
    ax.set_xticks(positions)
    ax.set_xticklabels(macro_pca["requested_pca_components"].astype(int).astype(str))
    ax.set_xlabel("PCA components")
    ax.set_ylabel("Macro C-index")
    ax.set_ylim(0.696, 0.708)
    ax.set_title("Representation dimension")
    style_axis(ax, "y")
    ax.legend(frameon=False, loc="lower right")
    panel_label(ax, "C")

    fig.tight_layout(w_pad=2.0)
    return save_all(fig, "ficpf_ablation_sensitivity")


def make_robustness_figure():
    source = (
        ROOT
        / "result"
        / "ficpf-omics-robustness-20260821"
        / "dataset_summary.csv"
    )
    data = pd.read_csv(source)
    datasets = ["BLCA", "BRCA", "GBMLGG", "LUAD", "UCEC", "MACRO"]
    conditions = [
        ("feature_dropout", "Feature dropout", [0.0, 0.25, 0.50, 0.75], ["0", "25", "50", "75"]),
        ("gaussian_noise", "Gaussian noise", [0.0, 0.10, 0.25, 0.50], ["0", "0.10", "0.25", "0.50"]),
        ("patient_missing", "Patient-level missing", [0.0, 0.25, 0.50, 0.75, 1.0], ["0", "25", "50", "75", "100"]),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(10.7, 3.25), sharey=True)
    for panel_index, (ax, (condition, title, severity, tick_labels)) in enumerate(
        zip(axes, conditions)
    ):
        positions = np.arange(len(severity), dtype=float)
        for dataset in datasets:
            clean = data[
                (data["dataset"] == dataset) & (data["condition"] == "clean")
            ].iloc[0]
            values = [float(clean["delta_vs_pathology"])]
            condition_rows = data[
                (data["dataset"] == dataset) & (data["condition"] == condition)
            ].set_index("severity")
            values.extend(
                float(condition_rows.loc[value, "delta_vs_pathology"])
                for value in severity[1:]
            )
            is_macro = dataset == "MACRO"
            ax.plot(
                positions,
                values,
                color=DATASET_COLORS[dataset],
                marker="o",
                markersize=4.2 if is_macro else 3.2,
                linewidth=2.2 if is_macro else 1.0,
                alpha=1.0 if is_macro else 0.75,
                label="Macro" if is_macro else dataset,
                zorder=4 if is_macro else 2,
            )
        ax.axhline(0.0, color="#4A4D50", linewidth=0.8)
        ax.set_xticks(positions)
        ax.set_xticklabels(tick_labels)
        if condition == "gaussian_noise":
            ax.set_xlabel("Noise standard deviation")
        else:
            ax.set_xlabel("Missing rate (%)")
        ax.set_title(title)
        ax.set_ylim(-0.017, 0.021)
        style_axis(ax, "y")
        panel_label(ax, chr(ord("A") + panel_index))
    axes[0].set_ylabel("C-index gain vs pathology")
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        frameon=False,
        ncol=6,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.05),
        columnspacing=1.2,
        handletextpad=0.4,
    )
    fig.tight_layout(w_pad=1.8, rect=(0.0, 0.0, 1.0, 0.94))
    return save_all(fig, "ficpf_omics_robustness")


def make_curve_quality_figure():
    source = (
        ROOT
        / "result"
        / "frozen-pathomic-fusion-five-dataset-summary-20260820"
        / "curve_dataset_summary.csv"
    )
    data = pd.read_csv(source).set_index("dataset")
    datasets = ["BLCA", "BRCA", "GBMLGG", "LUAD", "UCEC"]
    data = data.loc[datasets]
    x = np.arange(len(datasets), dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 6.0))
    panels = [
        (
            axes[0, 0],
            "baseline_ibs",
            "ficpf_ibs",
            "Integrated Brier score",
            "IBS (lower is better)",
            (0.08, 0.29),
        ),
        (
            axes[0, 1],
            "baseline_km_calibration",
            "ficpf_km_calibration",
            "Kaplan-Meier calibration",
            "Calibration error (lower is better)",
            (0.0, 0.066),
        ),
        (
            axes[1, 0],
            "baseline_mae",
            "ficpf_mae",
            "Absolute-time accuracy",
            "Hinge MAE in months (lower is better)",
            (15.0, 57.0),
        ),
    ]

    for panel_index, (ax, baseline_col, ficpf_col, title, ylabel, ylim) in enumerate(
        panels
    ):
        baseline = data[baseline_col].to_numpy(dtype=float)
        candidate = data[ficpf_col].to_numpy(dtype=float)
        for xi, base_value, candidate_value in zip(x, baseline, candidate):
            ax.plot(
                [xi - 0.12, xi + 0.12],
                [base_value, candidate_value],
                color="#B8BABC",
                linewidth=1.0,
                zorder=1,
            )
        ax.scatter(
            x - 0.12,
            baseline,
            color=COLORS["matched"],
            marker="o",
            s=34,
            label="Matched VLSA",
            zorder=3,
        )
        ax.scatter(
            x + 0.12,
            candidate,
            color=COLORS["ficpf"],
            marker="D",
            s=38,
            label="FICPF",
            zorder=3,
        )
        ax.set_xticks(x)
        ax.set_xticklabels(datasets)
        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        ax.set_title(title)
        style_axis(ax, "y")
        panel_label(ax, chr(ord("A") + panel_index))

    ax = axes[1, 1]
    width = 0.34
    baseline_pass = data["baseline_d_calibration_pass_folds"].to_numpy(dtype=float)
    candidate_pass = data["ficpf_d_calibration_pass_folds"].to_numpy(dtype=float)
    ax.bar(
        x - width / 2,
        baseline_pass,
        width=width,
        color=COLORS["matched"],
        label="Matched VLSA",
    )
    ax.bar(
        x + width / 2,
        candidate_pass,
        width=width,
        color=COLORS["ficpf"],
        label="FICPF",
    )
    for xi, base_value, candidate_value in zip(x, baseline_pass, candidate_pass):
        ax.text(xi - width / 2, base_value + 0.09, f"{int(base_value)}", ha="center", fontsize=7.5)
        ax.text(xi + width / 2, candidate_value + 0.09, f"{int(candidate_value)}", ha="center", fontsize=7.5)
    ax.set_xticks(x)
    ax.set_xticklabels(datasets)
    ax.set_ylabel("Folds passing D-calibration (of 5)")
    ax.set_ylim(0.0, 5.7)
    ax.set_title("Distributional calibration")
    style_axis(ax, "y")
    panel_label(ax, "D")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        frameon=False,
        ncol=2,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        columnspacing=1.6,
    )
    fig.tight_layout(w_pad=2.0, h_pad=2.0, rect=(0.0, 0.0, 1.0, 0.96))
    return save_all(fig, "ficpf_survival_curve_quality")


def main():
    setup_style()
    outputs = {
        "main_comparison": make_main_figure(),
        "ablation_sensitivity": make_ablation_figure(),
        "survival_curve_quality": make_curve_quality_figure(),
        "omics_robustness": make_robustness_figure(),
    }
    method_outputs = [
        OUT / f"ficpf_method_overview.{extension}"
        for extension in ("png", "pdf", "svg")
    ]
    if all(path.exists() for path in method_outputs):
        outputs = {
            "method_overview": [str(path) for path in method_outputs],
            **outputs,
        }
    manifest = {
        "metric_contract": "SurvivalEVAL_mean_continuous_comparable_all_ties",
        "main_result_source": str(
            ROOT
            / "result"
            / "external_multimodal_comparison_20260821"
            / "unified"
        ),
        "ablation_source": str(
            ROOT / "result" / "ficpf-ablation-sensitivity-20260821"
        ),
        "pca_source": str(ROOT / "result" / "ficpf-pca-sensitivity-20260821"),
        "robustness_source": str(
            ROOT / "result" / "ficpf-omics-robustness-20260821"
        ),
        "survival_curve_source": str(
            ROOT
            / "result"
            / "frozen-pathomic-fusion-five-dataset-summary-20260820"
            / "curve_dataset_summary.csv"
        ),
        "method_implementation": str(
            ROOT / "tools" / "evaluate_frozen_pathomic_fusion_outer5.py"
        ),
        "outputs": outputs,
    }
    (OUT / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(outputs, indent=2))


if __name__ == "__main__":
    main()
