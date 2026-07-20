#!/usr/bin/env python3
"""PCA and hierarchical clustering on k-mer frequency profiles."""

import argparse
import sys
from pathlib import Path

import pandas as pd
from scipy.cluster.hierarchy import linkage
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
from similarity import euclidean_distance_matrix


def run_pca(kmer_matrix: pd.DataFrame, n_components: int = 2) -> pd.DataFrame:
    """Standardize k-mer matrix and run PCA."""
    n_samples = kmer_matrix.shape[0]
    n_comp = min(n_components, n_samples, kmer_matrix.shape[1])
    scaler = StandardScaler()
    scaled = scaler.fit_transform(kmer_matrix.values)
    pca = PCA(n_components=n_comp)
    coords = pca.fit_transform(scaled)
    columns = [f"PC{i + 1}" for i in range(n_comp)]
    result = pd.DataFrame(coords, index=kmer_matrix.index, columns=columns)
    variance_df = pd.DataFrame(
        {
            "component": columns,
            "explained_variance_ratio": pca.explained_variance_ratio_,
        }
    )
    return result, variance_df


def run_hierarchical_clustering(dist_matrix: pd.DataFrame, method: str = "ward") -> pd.DataFrame:
    """Perform hierarchical clustering from a distance matrix."""
    condensed = squareform(dist_matrix.values, checks=False)
    if method == "ward":
        linkage_matrix = linkage(condensed, method="ward")
    else:
        linkage_matrix = linkage(condensed, method=method)
    return pd.DataFrame(linkage_matrix, columns=["idx1", "idx2", "distance", "count"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="PCA and hierarchical clustering.")
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("-k", type=int, default=DEFAULT_K)
    parser.add_argument(
        "--linkage",
        default="ward",
        choices=["ward", "single", "complete", "average"],
    )
    parser.add_argument("--demo", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        input_dir = resolve_input_dir(args.input_dir, demo=args.demo)
        sequences = load_sequences(input_dir)
        kmer_matrix = build_kmer_matrix(sequences, args.k)
        dist_matrix = euclidean_distance_matrix(kmer_matrix)

        pca_coords, variance_df = run_pca(kmer_matrix)
        linkage_df = run_hierarchical_clustering(dist_matrix, method=args.linkage)

        args.output_dir.mkdir(parents=True, exist_ok=True)
        pca_path = args.output_dir / "pca_coordinates.csv"
        variance_path = args.output_dir / "pca_variance.csv"
        linkage_path = args.output_dir / "linkage_matrix.csv"

        pca_coords.to_csv(pca_path)
        variance_df.to_csv(variance_path, index=False)
        linkage_df.to_csv(linkage_path, index=False)

        print(f"PCA coordinates saved to {pca_path}")
        print(f"PCA variance saved to {variance_path}")
        print(f"Linkage matrix saved to {linkage_path}")
        return 0
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
