#!/usr/bin/env python3
"""Euclidean distance and cosine similarity between k-mer profiles."""

import argparse
import sys
from pathlib import Path

import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.metrics.pairwise import cosine_similarity

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


def euclidean_distance_matrix(kmer_matrix: pd.DataFrame) -> pd.DataFrame:
    """Pairwise Euclidean distance between k-mer frequency vectors."""
    labels = kmer_matrix.index.tolist()
    distances = cdist(kmer_matrix.values, kmer_matrix.values, metric="euclidean")
    return pd.DataFrame(distances, index=labels, columns=labels)


def cosine_similarity_matrix(kmer_matrix: pd.DataFrame) -> pd.DataFrame:
    """Pairwise cosine similarity between k-mer frequency vectors."""
    labels = kmer_matrix.index.tolist()
    sim = cosine_similarity(kmer_matrix.values)
    return pd.DataFrame(sim, index=labels, columns=labels)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="K-mer profile similarity and distance matrices."
    )
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
        dist_path = args.output_dir / "distance_matrix.csv"
        cos_path = args.output_dir / "cosine_similarity.csv"
        dist_matrix.to_csv(dist_path)
        cos_matrix.to_csv(cos_path)
        print(f"Distance matrix saved to {dist_path}")
        print(f"Cosine similarity matrix saved to {cos_path}")
        return 0
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
