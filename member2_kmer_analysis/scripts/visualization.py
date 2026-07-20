#!/usr/bin/env python3
"""Generate heatmaps, PCA plots, dendrograms, and comparison figures."""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from kmer_analysis import (
    DEFAULT_INPUT_DIR,
    DEFAULT_K,
    DEFAULT_OUTPUT_DIR,
    build_kmer_matrix,
    load_sequences,
    resolve_input_dir,
)
from similarity import cosine_similarity_matrix, euclidean_distance_matrix

FIGURES_DIR_NAME = "figures"
DPI = 150


def ensure_figures_dir(output_dir: Path) -> Path:
    figures_dir = output_dir / FIGURES_DIR_NAME
    figures_dir.mkdir(parents=True, exist_ok=True)
    return figures_dir


def plot_heatmap(dist_matrix: pd.DataFrame, figures_dir: Path) -> Path:
    """Euclidean distance heatmap."""
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(
        dist_matrix,
        annot=True,
        fmt=".3f",
        cmap="YlOrRd",
        square=True,
        linewidths=0.5,
        ax=ax,
        cbar_kws={"label": "Euclidean Distance"},
    )
    ax.set_title("K-mer Profile Distance Heatmap")
    plt.tight_layout()
    path = figures_dir / "heatmap.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_similarity_matrix(cos_matrix: pd.DataFrame, figures_dir: Path) -> Path:
    """Cosine similarity heatmap."""
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(
        cos_matrix,
        annot=True,
        fmt=".3f",
        cmap="viridis",
        square=True,
        linewidths=0.5,
        vmin=0,
        vmax=1,
        ax=ax,
        cbar_kws={"label": "Cosine Similarity"},
    )
    ax.set_title("K-mer Profile Cosine Similarity Matrix")
    plt.tight_layout()
    path = figures_dir / "similarity_matrix.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_pca(kmer_matrix: pd.DataFrame, figures_dir: Path) -> Path:
    """PCA scatter plot of k-mer profiles."""
    n_samples = kmer_matrix.shape[0]
    n_comp = min(2, n_samples, kmer_matrix.shape[1])
    scaler = StandardScaler()
    scaled = scaler.fit_transform(kmer_matrix.values)
    pca = PCA(n_components=n_comp)
    coords = pca.fit_transform(scaled)
    labels = kmer_matrix.index.tolist()

    fig, ax = plt.subplots(figsize=(8, 6))
    if n_comp >= 2:
        ax.scatter(coords[:, 0], coords[:, 1], s=100, c="steelblue", edgecolors="black")
        for i, label in enumerate(labels):
            ax.annotate(label, (coords[i, 0], coords[i, 1]), xytext=(5, 5), textcoords="offset points")
        var1 = pca.explained_variance_ratio_[0] * 100
        var2 = pca.explained_variance_ratio_[1] * 100
        ax.set_xlabel(f"PC1 ({var1:.1f}% variance)")
        ax.set_ylabel(f"PC2 ({var2:.1f}% variance)")
    else:
        ax.bar(range(len(labels)), coords[:, 0], color="steelblue")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_ylabel("PC1")

    ax.set_title("PCA of K-mer Frequency Profiles")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    path = figures_dir / "pca.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_dendrogram(dist_matrix: pd.DataFrame, figures_dir: Path) -> Path:
    """Hierarchical clustering dendrogram."""
    condensed = squareform(dist_matrix.values, checks=False)
    Z = linkage(condensed, method="ward")

    fig, ax = plt.subplots(figsize=(10, 6))
    dendrogram(
        Z,
        labels=dist_matrix.index.tolist(),
        leaf_rotation=45,
        leaf_font_size=10,
        ax=ax,
    )
    ax.set_title("Hierarchical Clustering Dendrogram (Ward Linkage)")
    ax.set_ylabel("Euclidean Distance")
    plt.tight_layout()
    path = figures_dir / "dendrogram.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_entropy_comparison(entropy_path: Path, figures_dir: Path) -> Path | None:
    """Bar chart comparing Shannon entropy across sequences."""
    if not entropy_path.exists():
        return None
    entropy_df = pd.read_csv(entropy_path)
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.barplot(
        data=entropy_df,
        x="sequence_id",
        y="shannon_entropy",
        hue="sequence_id",
        palette="muted",
        legend=False,
        ax=ax,
    )
    ax.set_title("Shannon Entropy Comparison (K-mer Distribution)")
    ax.set_xlabel("Sequence")
    ax.set_ylabel("Shannon Entropy (bits)")
    ax.tick_params(axis="x", rotation=45)
    plt.tight_layout()
    path = figures_dir / "entropy_comparison.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate analysis figures.")
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("-k", type=int, default=DEFAULT_K)
    parser.add_argument("--demo", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        input_dir = resolve_input_dir(args.input_dir, demo=args.demo)
        sequences = load_sequences(input_dir)
        kmer_matrix = build_kmer_matrix(sequences, args.k)
        dist_matrix = euclidean_distance_matrix(kmer_matrix)
        cos_matrix = cosine_similarity_matrix(kmer_matrix)

        args.output_dir.mkdir(parents=True, exist_ok=True)
        figures_dir = ensure_figures_dir(args.output_dir)

        saved = []
        saved.append(plot_heatmap(dist_matrix, figures_dir))
        saved.append(plot_similarity_matrix(cos_matrix, figures_dir))
        saved.append(plot_pca(kmer_matrix, figures_dir))
        saved.append(plot_dendrogram(dist_matrix, figures_dir))

        entropy_path = args.output_dir / "entropy.csv"
        entropy_fig = plot_entropy_comparison(entropy_path, figures_dir)
        if entropy_fig:
            saved.append(entropy_fig)
        else:
            print(
                "Note: entropy.csv not found. Run entropy.py first for entropy comparison figure.",
                file=sys.stderr,
            )

        for path in saved:
            print(f"Figure saved to {path}")
        return 0
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
