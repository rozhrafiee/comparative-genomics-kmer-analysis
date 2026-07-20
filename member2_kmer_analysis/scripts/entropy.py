#!/usr/bin/env python3
"""Shannon entropy of k-mer distributions."""

import argparse
import math
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from kmer_analysis import (
    DEFAULT_INPUT_DIR,
    DEFAULT_K,
    DEFAULT_OUTPUT_DIR,
    count_kmers,
    create_demo_data,
    load_sequences,
    resolve_input_dir,
)


def shannon_entropy(counts: dict[str, int]) -> float:
    """Compute Shannon entropy (bits) from k-mer count dictionary."""
    total = sum(counts.values())
    if total == 0:
        return 0.0
    entropy = 0.0
    for count in counts.values():
        if count > 0:
            p = count / total
            entropy -= p * math.log2(p)
    return entropy


def compute_entropy_table(sequences: dict[str, str], k: int) -> pd.DataFrame:
    """Calculate Shannon entropy per sequence."""
    records = []
    for seq_id, sequence in sequences.items():
        counts = count_kmers(sequence, k)
        records.append(
            {
                "sequence_id": seq_id,
                "k": k,
                "total_kmers": sum(counts.values()),
                "unique_kmers": len(counts),
                "shannon_entropy": shannon_entropy(dict(counts)),
            }
        )
    return pd.DataFrame(records)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Shannon entropy of k-mer profiles.")
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
        entropy_df = compute_entropy_table(sequences, args.k)

        args.output_dir.mkdir(parents=True, exist_ok=True)
        output_path = args.output_dir / "entropy.csv"
        entropy_df.to_csv(output_path, index=False)
        print(f"Entropy results saved to {output_path}")
        return 0
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
