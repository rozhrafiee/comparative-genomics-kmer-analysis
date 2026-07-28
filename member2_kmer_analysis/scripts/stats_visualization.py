#!/usr/bin/env python3
"""Visualizations for stats_tests.py outputs — correlation scatter plots,
per-organism entropy/diversity distributions, ANOVA/Kruskal-Wallis summary,
and a Mann-Whitney pairwise significance heatmap.

Reads CSV outputs from results/ produced by stats_tests.py — no FASTA
parsing, so this script has no Biopython dependency and can be run any
time after stats_tests.py. Figures are saved to results/figures/,
alongside the rest of the pipeline's figures.

Usage:
    python scripts/stats_visualization.py
    python scripts/stats_visualization.py --results-dir results
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parent
DEFAULT_RESULTS_DIR = MODULE_DIR / "results"


DPI = 150


def format_label(name: str) -> str:
    return str(name).replace("_", " ")


def plot_correlation_matrix_heatmap(stats_dir: Path, figures_dir: Path) -> Path | None:
    """Heatmap of Pearson correlations for بخش پنجم variables."""
    path = stats_dir / "species_level_metrics.csv"
    if not path.exists():
        print(f"Skip correlation matrix: {path} not found.", file=sys.stderr)
        return None

    df = pd.read_csv(path)
    columns = ["genome_size", "gc_percent", "shannon_entropy", "kmer_diversity"]
    available = [col for col in columns if col in df.columns and df[col].notna().any()]
    if len(available) < 2:
        print("Skip correlation matrix: fewer than 2 numeric columns available.",
              file=sys.stderr)
        return None

    sub = df[available].astype(float)
    corr = sub.corr(method="pearson")
    labels = [format_label(col) for col in corr.columns]

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.heatmap(
        corr,
        annot=True,
        fmt=".2f",
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        square=True,
        linewidths=0.5,
        xticklabels=labels,
        yticklabels=labels,
        ax=ax,
        cbar_kws={"label": "Pearson r"},
    )
    ax.set_title("بخش پنجم — Correlation Matrix (species-level, n = organisms)")
    fig.tight_layout()
    out_path = figures_dir / "correlation_matrix.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_correlation_scatter(stats_dir: Path, figures_dir: Path) -> Path | None:
    path = stats_dir / "species_level_metrics.csv"
    corr_path = stats_dir / "correlation_results.csv"
    if not path.exists():
        print(f"Skip correlation scatter: {path} not found.", file=sys.stderr)
        return None
    df = pd.read_csv(path)
    corr_df = pd.read_csv(corr_path) if corr_path.exists() else pd.DataFrame()

    pairs = [("genome_size", "shannon_entropy"), ("genome_size", "kmer_diversity")]
    if "gc_percent" in df.columns and df["gc_percent"].notna().any():
        pairs += [("gc_percent", "shannon_entropy"), ("gc_percent", "kmer_diversity")]

    n = len(pairs)
    ncols = 2
    nrows = (n + 1) // 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
    axes = np.array(axes).reshape(-1)

    for ax, (x_col, y_col) in zip(axes, pairs):
        sub = df[[x_col, y_col, "organism"]].dropna()
        ax.scatter(sub[x_col], sub[y_col], s=90, color="steelblue", edgecolor="black", zorder=3)
        for _, row in sub.iterrows():
            ax.annotate(format_label(row["organism"]), (row[x_col], row[y_col]),
                        xytext=(5, 5), textcoords="offset points", fontsize=8)
        if len(sub) >= 3:
            slope, intercept = np.polyfit(sub[x_col], sub[y_col], 1)
            xs = np.linspace(sub[x_col].min(), sub[x_col].max(), 50)
            ax.plot(xs, slope * xs + intercept, color="darkorange", linestyle="--", zorder=2)

        r_txt = ""
        if not corr_df.empty:
            match = corr_df[(corr_df["variable_x"] == x_col) & (corr_df["variable_y"] == y_col)]
            if not match.empty and pd.notna(match.iloc[0]["pearson_r"]):
                r = match.iloc[0]["pearson_r"]
                p = match.iloc[0]["pearson_p"]
                r_txt = f"Pearson r={r:.2f}, p={p:.3f}"
        ax.set_title(f"{format_label(x_col)} vs {format_label(y_col)}\n{r_txt}", fontsize=10)
        ax.set_xlabel(format_label(x_col))
        ax.set_ylabel(format_label(y_col))
        ax.grid(alpha=0.3)

    for ax in axes[n:]:
        ax.axis("off")

    fig.suptitle("Species-Level Correlations (exploratory, n = number of organisms)", fontsize=12)
    fig.tight_layout()
    out_path = figures_dir / "correlation_scatter.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_metric_boxplot(stats_dir: Path, figures_dir: Path, metric: str, fname: str, title: str) -> Path | None:
    path = stats_dir / "per_record_metrics.csv"
    if not path.exists():
        print(f"Skip {fname}: {path} not found.", file=sys.stderr)
        return None
    df = pd.read_csv(path)
    df["organism_label"] = df["organism"].map(format_label)
    order = df.groupby("organism_label")[metric].median().sort_values(ascending=False).index

    fig, ax = plt.subplots(figsize=(9, 6))
    sns.boxplot(data=df, x="organism_label", y=metric, order=order,
                hue="organism_label", palette="Set2", legend=False, ax=ax)
    sns.stripplot(data=df, x="organism_label", y=metric, order=order,
                  color="black", size=3, alpha=0.5, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Organism")
    ax.set_ylabel(format_label(metric))
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    out_path = figures_dir / fname
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_anova_kruskal_summary(stats_dir: Path, figures_dir: Path) -> Path | None:
    path = stats_dir / "anova_kruskal_results.csv"
    if not path.exists():
        print(f"Skip ANOVA/Kruskal summary: {path} not found.", file=sys.stderr)
        return None
    df = pd.read_csv(path).dropna(subset=["anova_p", "kruskal_p"])
    if df.empty:
        print("Skip ANOVA/Kruskal summary: no valid rows.", file=sys.stderr)
        return None

    df["metric_label"] = df["metric"].map(format_label)
    melted = df.melt(
        id_vars="metric_label", value_vars=["anova_p", "kruskal_p"],
        var_name="test", value_name="p_value",
    )
    melted["neg_log10_p"] = -np.log10(melted["p_value"].clip(lower=1e-300))
    melted["test"] = melted["test"].map({"anova_p": "ANOVA", "kruskal_p": "Kruskal-Wallis"})

    fig, ax = plt.subplots(figsize=(7, 5))
    sns.barplot(data=melted, x="metric_label", y="neg_log10_p", hue="test", palette="muted", ax=ax)
    sig_line = -np.log10(0.05)
    ax.axhline(sig_line, color="red", linestyle="--", linewidth=1, label="p = 0.05")
    ax.set_title("ANOVA & Kruskal-Wallis: differences across organisms")
    ax.set_xlabel("Metric")
    ax.set_ylabel("-log10(p-value)")
    ax.legend()
    fig.tight_layout()
    out_path = figures_dir / "anova_kruskal_summary.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_mannwhitney_heatmap(stats_dir: Path, figures_dir: Path) -> Path | None:
    path = stats_dir / "mannwhitney_pairwise.csv"
    if not path.exists():
        print(f"Skip Mann-Whitney heatmap: {path} not found.", file=sys.stderr)
        return None
    df = pd.read_csv(path).dropna(subset=["p_value"])
    if df.empty:
        print("Skip Mann-Whitney heatmap: no valid rows.", file=sys.stderr)
        return None

    metrics = df["metric"].unique()
    fig, axes = plt.subplots(1, len(metrics), figsize=(7 * len(metrics), 6))
    axes = np.atleast_1d(axes)

    for ax, metric in zip(axes, metrics):
        sub = df[df["metric"] == metric]
        organisms = sorted(set(sub["organism_a"]).union(sub["organism_b"]))
        labels = [format_label(o) for o in organisms]
        mat = pd.DataFrame(np.nan, index=labels, columns=labels)
        for _, row in sub.iterrows():
            a, b = format_label(row["organism_a"]), format_label(row["organism_b"])
            mat.loc[a, b] = row["p_value"]
            mat.loc[b, a] = row["p_value"]
        sns.heatmap(mat, annot=True, fmt=".3f", cmap="coolwarm_r", vmin=0, vmax=0.1,
                    linewidths=0.5, ax=ax, cbar_kws={"label": "p-value"})
        ax.set_title(f"Pairwise Mann-Whitney U — {format_label(metric)}")
        ax.tick_params(axis="x", rotation=45)
        ax.tick_params(axis="y", rotation=0)

    fig.tight_layout()
    out_path = figures_dir / "mannwhitney_heatmap.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Visualize stats_tests.py outputs.")
    parser.add_argument(
        "--results-dir", type=Path, default=DEFAULT_RESULTS_DIR,
        help="Directory containing stats CSVs; figures/ is created under here.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        results_dir = args.results_dir
        figures_dir = results_dir / "figures"
        figures_dir.mkdir(parents=True, exist_ok=True)

        if not results_dir.exists():
            print(
                f"Error: {results_dir} not found. Run scripts/stats_tests.py first.",
                file=sys.stderr,
            )
            return 1

        saved = []
        saved.append(plot_correlation_matrix_heatmap(results_dir, figures_dir))
        saved.append(plot_correlation_scatter(results_dir, figures_dir))
        saved.append(plot_metric_boxplot(
            results_dir, figures_dir, "shannon_entropy",
            "entropy_by_organism_boxplot.png",
            "Shannon Entropy per Chromosome/Scaffold, Grouped by Organism",
        ))
        saved.append(plot_metric_boxplot(
            results_dir, figures_dir, "kmer_diversity",
            "diversity_by_organism_boxplot.png",
            "K-mer Diversity per Chromosome/Scaffold, Grouped by Organism",
        ))
        saved.append(plot_anova_kruskal_summary(results_dir, figures_dir))
        saved.append(plot_mannwhitney_heatmap(results_dir, figures_dir))

        for path in saved:
            if path:
                print(f"Figure saved to {path}")
        return 0
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())