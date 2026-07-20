#!/usr/bin/env python3
"""TATA box motif detection in promoter regions."""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from kmer_analysis import (
    DEFAULT_INPUT_DIR,
    DEFAULT_OUTPUT_DIR,
    load_sequences,
    normalize_sequence,
    resolve_input_dir,
)

TATA_PATTERNS = [
    re.compile(r"TATAAA"),
    re.compile(r"TATAWA", re.IGNORECASE),  # W = A or T
    re.compile(r"TATA"),
]

DEFAULT_PROMOTER_LENGTH = 500


def find_tata_boxes(
    sequence: str,
    promoter_length: int = DEFAULT_PROMOTER_LENGTH,
) -> list[dict]:
    """Search for TATA box motifs in the 5' promoter region."""
    seq = normalize_sequence(sequence)
    promoter = seq[:promoter_length]
    hits = []
    for pattern in TATA_PATTERNS:
        for match in pattern.finditer(promoter):
            hits.append(
                {
                    "motif": match.group(),
                    "pattern": pattern.pattern,
                    "position": match.start(),
                    "position_1based": match.start() + 1,
                }
            )
    return hits


def analyze_tata_boxes(
    sequences: dict[str, str],
    promoter_length: int = DEFAULT_PROMOTER_LENGTH,
) -> pd.DataFrame:
    """Build TATA box summary table for all sequences."""
    records = []
    for seq_id, sequence in sequences.items():
        hits = find_tata_boxes(sequence, promoter_length)
        canonical = [h for h in hits if h["pattern"] == "TATAAA"]
        records.append(
            {
                "sequence_id": seq_id,
                "promoter_length": promoter_length,
                "tata_box_count": len(canonical),
                "tata_variant_count": len(hits),
                "has_canonical_tata": len(canonical) > 0,
                "first_tata_position": canonical[0]["position_1based"] if canonical else None,
                "first_tata_motif": canonical[0]["motif"] if canonical else None,
                "all_tata_positions": ";".join(
                    str(h["position_1based"]) for h in canonical
                ),
            }
        )
    return pd.DataFrame(records)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="TATA box motif analysis.")
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--promoter-length",
        type=int,
        default=DEFAULT_PROMOTER_LENGTH,
        help="Length of 5' region to scan (bp).",
    )
    parser.add_argument("--demo", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        input_dir = resolve_input_dir(args.input_dir, demo=args.demo)
        sequences = load_sequences(input_dir)
        tata_df = analyze_tata_boxes(sequences, args.promoter_length)

        args.output_dir.mkdir(parents=True, exist_ok=True)
        output_path = args.output_dir / "tata_box.csv"
        tata_df.to_csv(output_path, index=False)
        print(f"TATA box results saved to {output_path}")
        return 0
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
