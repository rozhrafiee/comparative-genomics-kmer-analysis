#!/usr/bin/env python3
"""K-mer counting and summary statistics for comparative genomics."""

import argparse
import sys
from collections import Counter
from itertools import product
from pathlib import Path

import pandas as pd
from Bio import SeqIO

SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parent
DEFAULT_INPUT_DIR = MODULE_DIR.parent / "data"
DEFAULT_DEMO_DIR = MODULE_DIR / "results" / "demo_input"
DEFAULT_OUTPUT_DIR = MODULE_DIR / "results"
DEFAULT_K = 4


def get_all_kmers(k: int) -> list[str]:
    """Return all possible DNA k-mers for a given k."""
    return ["".join(bases) for bases in product("ACGT", repeat=k)]


def normalize_sequence(sequence: str) -> str:
    """Uppercase and retain valid DNA bases only."""
    return "".join(base for base in sequence.upper() if base in "ACGT")


def count_kmers(sequence: str, k: int) -> Counter:
    """Count overlapping k-mers in a DNA sequence."""
    seq = normalize_sequence(sequence)
    counts: Counter = Counter()
    if len(seq) < k:
        return counts
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if "N" not in kmer:
            counts[kmer] += 1
    return counts


def load_sequences(input_dir: Path) -> dict[str, str]:
    """Load all FASTA sequences from a directory."""
    sequences: dict[str, str] = {}
    fasta_files = sorted(input_dir.glob("*.fa")) + sorted(input_dir.glob("*.fasta"))
    if not fasta_files:
        raise FileNotFoundError(
            f"No FASTA files found in {input_dir}. "
            "Place .fa or .fasta files there, or use --demo."
        )
    for fasta_path in fasta_files:
        try:
            for record in SeqIO.parse(fasta_path, "fasta"):
                seq_id = record.id
                if seq_id in sequences:
                    seq_id = f"{seq_id}_{fasta_path.stem}"
                sequences[seq_id] = str(record.seq)
        except Exception as exc:
            raise RuntimeError(f"Failed to parse {fasta_path}: {exc}") from exc
    if not sequences:
        raise ValueError(f"No sequences loaded from {input_dir}")
    return sequences


def kmer_frequency_vector(counts: Counter, k: int) -> dict[str, float]:
    """Convert k-mer counts to normalized frequency vector."""
    all_kmers = get_all_kmers(k)
    total = sum(counts.values())
    if total == 0:
        return {kmer: 0.0 for kmer in all_kmers}
    return {kmer: counts.get(kmer, 0) / total for kmer in all_kmers}


def build_kmer_matrix(sequences: dict[str, str], k: int) -> pd.DataFrame:
    """Build a k-mer frequency matrix (samples x k-mers)."""
    rows = []
    for seq_id, sequence in sequences.items():
        counts = count_kmers(sequence, k)
        rows.append({"sequence_id": seq_id, **kmer_frequency_vector(counts, k)})
    return pd.DataFrame(rows).set_index("sequence_id")


def summarize_kmers(sequences: dict[str, str], k: int) -> pd.DataFrame:
    """Generate per-sequence k-mer summary statistics."""
    records = []
    for seq_id, sequence in sequences.items():
        counts = count_kmers(sequence, k)
        total_kmers = sum(counts.values())
        unique_kmers = len(counts)
        most_common = counts.most_common(1)[0] if counts else ("", 0)
        records.append(
            {
                "sequence_id": seq_id,
                "sequence_length": len(normalize_sequence(sequence)),
                "k": k,
                "total_kmers": total_kmers,
                "unique_kmers": unique_kmers,
                "most_common_kmer": most_common[0],
                "most_common_count": most_common[1],
                "kmer_diversity": unique_kmers / total_kmers if total_kmers else 0.0,
            }
        )
    return pd.DataFrame(records)


def resolve_input_dir(input_dir: Path, demo: bool = False) -> Path:
    """Resolve FASTA input directory, optionally creating or falling back to demo data."""
    if demo:
        return create_demo_data()

    fasta_present = bool(
        list(input_dir.glob("*.fa")) + list(input_dir.glob("*.fasta"))
    )
    if fasta_present:
        return input_dir

    demo_present = bool(
        list(DEFAULT_DEMO_DIR.glob("*.fa")) + list(DEFAULT_DEMO_DIR.glob("*.fasta"))
    )
    if demo_present:
        return DEFAULT_DEMO_DIR

    return input_dir


def create_demo_data(output_dir: Path | None = None) -> Path:
    """Write demo FASTA sequences for pipeline testing."""
    output_dir = output_dir or DEFAULT_DEMO_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    demo_path = output_dir / "demo_sequences.fasta"
    sequences = {
        "species_A": (
            "ATCGATCGATCGTATAAAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        ),
        "species_B": (
            "GCTAGCTAGCTAGCTATAAAAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
            "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
        ),
        "species_C": (
            "TTTTAAAACCCCGGGGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
            "AAAATTTTCCCCGGGGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        ),
        "species_D": (
            "CGCGCGCGCGCGTATAAACGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG"
            "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"
        ),
    }
    with demo_path.open("w", encoding="utf-8") as handle:
        for seq_id, seq in sequences.items():
            handle.write(f">{seq_id}\n{seq}\n")
    return output_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="K-mer counting and summary.")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Directory containing FASTA files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for CSV outputs.",
    )
    parser.add_argument("-k", type=int, default=DEFAULT_K, help="K-mer size.")
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Generate demo FASTA data and run analysis.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        input_dir = resolve_input_dir(args.input_dir, demo=args.demo)
        if args.demo:
            print(f"Demo data written to {input_dir}")

        sequences = load_sequences(input_dir)
        summary = summarize_kmers(sequences, args.k)

        args.output_dir.mkdir(parents=True, exist_ok=True)
        output_path = args.output_dir / "kmer_summary.csv"
        summary.to_csv(output_path, index=False)
        print(f"K-mer summary saved to {output_path} ({len(summary)} sequences)")
        return 0
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
