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
MULTI_K_VALUES = (3, 4, 5)


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


def count_kmers_from_fasta(fasta_path: Path, k: int) -> Counter:
    """Count k-mers across all records in a (possibly multi-contig) FASTA file."""
    total_counts: Counter = Counter()
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            total_counts.update(count_kmers(str(record.seq), k))
    except Exception as exc:
        raise RuntimeError(f"Failed to count k-mers in {fasta_path}: {exc}") from exc
    return total_counts


def format_organism_label(sequence_id: str) -> str:
    """Convert stored sequence IDs into readable organism labels for plots."""
    cleaned = sequence_id.replace("_", " ").strip()
    # Drop common strain suffixes for cleaner figure labels.
    for suffix in (
        " OR74A",
        " AX4",
        " CBS 5749",
        " Bristol N2",
        " CLIB89(W29)",
    ):
        if cleaned.endswith(suffix):
            cleaned = cleaned[: -len(suffix)].strip()
    return cleaned


def read_organism_name(dataset_dir: Path) -> str:
    """Read organism scientific name from NCBI data_summary.tsv."""
    summary_files = list(dataset_dir.rglob("data_summary.tsv"))
    if not summary_files:
        return dataset_dir.name
    try:
        with summary_files[0].open(encoding="utf-8") as handle:
            lines = [line.strip() for line in handle if line.strip()]
        if len(lines) >= 2:
            organism = lines[1].split("\t")[0].strip()
            if organism:
                # Keep binomial name (Genus species); drop strain suffixes.
                parts = organism.split()
                if len(parts) >= 2:
                    organism = f"{parts[0]} {parts[1]}"
                return organism.replace(" ", "_")
    except OSError:
        pass
    return dataset_dir.name


def select_genomic_fasta(dataset_dir: Path) -> Path | None:
    """Pick RefSeq (GCF) or GenBank (GCA) genomic FASTA from an NCBI dataset folder."""
    fna_files = sorted(dataset_dir.rglob("*_genomic.fna"))
    gcf_files = [path for path in fna_files if path.name.startswith("GCF_")]
    gca_files = [path for path in fna_files if path.name.startswith("GCA_")]
    if gcf_files:
        return gcf_files[0]
    if gca_files:
        return gca_files[0]
    return None


def discover_ncbi_genomes(input_dir: Path) -> dict[str, Path]:
    """Map organism names to genomic FASTA paths under ncbi_dataset* folders."""
    genomes: dict[str, Path] = {}
    for dataset_dir in sorted(input_dir.glob("ncbi_dataset*")):
        if not dataset_dir.is_dir():
            continue
        fasta_path = select_genomic_fasta(dataset_dir)
        if fasta_path is None:
            continue
        organism = read_organism_name(dataset_dir)
        key = organism
        suffix = 2
        while key in genomes:
            key = f"{organism}_{suffix}"
            suffix += 1
        genomes[key] = fasta_path
    return genomes


def load_pooled_sequence(fasta_path: Path) -> str:
    """Concatenate all contigs/chromosomes from a genomic FASTA into one sequence."""
    parts: list[str] = []
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            parts.append(normalize_sequence(str(record.seq)))
    except Exception as exc:
        raise RuntimeError(f"Failed to parse {fasta_path}: {exc}") from exc
    if not parts:
        raise ValueError(f"No sequences found in {fasta_path}")
    return "".join(parts)


def load_sequences_by_record(input_dir: Path) -> dict[str, list[tuple[str, str]]]:
    """Load per-chromosome/scaffold sequences grouped by organism (unpooled).

    Unlike load_sequences()/load_pooled_sequence(), this keeps each FASTA
    record separate so per-organism *groups* of values (e.g. one Shannon
    entropy value per chromosome) can be built for statistical tests such
    as ANOVA / Kruskal-Wallis that need multiple samples per group.

    Returns: {organism: [(record_id, sequence), ...], ...}
    """
    ncbi_genomes = discover_ncbi_genomes(input_dir)
    grouped: dict[str, list[tuple[str, str]]] = {}

    if ncbi_genomes:
        for organism, fasta_path in ncbi_genomes.items():
            records: list[tuple[str, str]] = []
            try:
                for record in SeqIO.parse(fasta_path, "fasta"):
                    seq = normalize_sequence(str(record.seq))
                    if seq:
                        records.append((record.id, seq))
            except Exception as exc:
                raise RuntimeError(f"Failed to parse {fasta_path}: {exc}") from exc
            if records:
                grouped[organism] = records
        return grouped

    fasta_files = sorted(
        list(input_dir.glob("*.fa"))
        + list(input_dir.glob("*.fasta"))
        + list(input_dir.glob("*.fna"))
    )
    if not fasta_files:
        raise FileNotFoundError(
            f"No FASTA/FNA files or ncbi_dataset* folders found in {input_dir}. "
            "Place NCBI datasets under data/ncbi_dataset1..N, or use --demo."
        )
    for fasta_path in fasta_files:
        organism = fasta_path.stem
        records: list[tuple[str, str]] = []
        try:
            for record in SeqIO.parse(fasta_path, "fasta"):
                seq = normalize_sequence(str(record.seq))
                if seq:
                    records.append((record.id, seq))
        except Exception as exc:
            raise RuntimeError(f"Failed to parse {fasta_path}: {exc}") from exc
        if records:
            grouped[organism] = records
    if not grouped:
        raise ValueError(f"No sequences loaded from {input_dir}")
    return grouped


def has_sequence_data(input_dir: Path) -> bool:
    """Return True if input_dir contains usable FASTA or NCBI dataset folders."""
    if discover_ncbi_genomes(input_dir):
        return True
    extensions = ("*.fa", "*.fasta", "*.fna")
    return any(files for ext in extensions for files in input_dir.glob(ext))


def load_sequences(input_dir: Path) -> dict[str, str]:
    """Load sequences from NCBI dataset folders or flat FASTA files."""
    ncbi_genomes = discover_ncbi_genomes(input_dir)
    if ncbi_genomes:
        sequences: dict[str, str] = {}
        for organism, fasta_path in ncbi_genomes.items():
            print(f"Loading {organism} from {fasta_path.name}...")
            sequences[organism] = load_pooled_sequence(fasta_path)
        return sequences

    sequences: dict[str, str] = {}
    fasta_files = sorted(
        list(input_dir.glob("*.fa"))
        + list(input_dir.glob("*.fasta"))
        + list(input_dir.glob("*.fna"))
    )
    if not fasta_files:
        raise FileNotFoundError(
            f"No FASTA/FNA files or ncbi_dataset* folders found in {input_dir}. "
            "Place NCBI datasets under data/ncbi_dataset1..N, or use --demo."
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


def organism_kmer_set(sequence: str, k: int) -> set[str]:
    """Return the set of k-mers with non-zero count in a sequence."""
    return set(count_kmers(sequence, k).keys())


def build_kmer_top_by_k_table(
    sequences: dict[str, str], k_values: tuple[int, ...] = MULTI_K_VALUES
) -> pd.DataFrame:
    """Most common k-mer per organism for each k value."""
    records = []
    for k in k_values:
        for seq_id, sequence in sequences.items():
            counts = count_kmers(sequence, k)
            most_common = counts.most_common(1)[0] if counts else ("", 0)
            records.append(
                {
                    "k": k,
                    "organism": seq_id,
                    "most_common_kmer": most_common[0],
                    "most_common_count": most_common[1],
                }
            )
    return pd.DataFrame(records)


def find_species_unique_kmers(
    sequences: dict[str, str], k_values: tuple[int, ...] = MULTI_K_VALUES
) -> pd.DataFrame:
    """K-mers observed in exactly one organism (set difference across all sets)."""
    records = []
    for k in k_values:
        organism_sets = {
            org: organism_kmer_set(seq, k) for org, seq in sequences.items()
        }
        kmer_owners: dict[str, list[str]] = {}
        for org, kmers in organism_sets.items():
            for kmer in kmers:
                kmer_owners.setdefault(kmer, []).append(org)
        for kmer, owners in sorted(kmer_owners.items()):
            if len(owners) == 1:
                records.append(
                    {
                        "k": k,
                        "kmer": kmer,
                        "exclusive_organism": owners[0],
                    }
                )
    return pd.DataFrame(records)


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

    if has_sequence_data(input_dir):
        return input_dir

    demo_present = bool(
        list(DEFAULT_DEMO_DIR.glob("*.fa"))
        + list(DEFAULT_DEMO_DIR.glob("*.fasta"))
        + list(DEFAULT_DEMO_DIR.glob("*.fna"))
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
        "Caenorhabditis_elegans": (
            "ATCGATCGATCGTATAAAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        ),
        "Neurospora_crassa": (
            "GCTAGCTAGCTAGCTATAAAAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
            "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
        ),
        "Yarrowia_lipolytica": (
            "TTTTAAAACCCCGGGGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
            "AAAATTTTCCCCGGGGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        ),
        "Dictyostelium_discoideum": (
            "CGCGCGCGCGCGTATAAACGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG"
            "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"
        ),
        "Eremothecium_coryli": (
            "ACGTACGTACGTTATAAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
            "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"
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
        top_by_k = build_kmer_top_by_k_table(sequences)
        species_unique = find_species_unique_kmers(sequences)

        args.output_dir.mkdir(parents=True, exist_ok=True)
        output_path = args.output_dir / "kmer_summary.csv"
        summary.to_csv(output_path, index=False)
        print(f"K-mer summary saved to {output_path} ({len(summary)} sequences)")

        top_by_k_path = args.output_dir / "kmer_top_by_k.csv"
        top_by_k.to_csv(top_by_k_path, index=False)
        print(f"Multi-k top k-mers saved to {top_by_k_path} ({len(top_by_k)} rows)")

        species_unique_path = args.output_dir / "species_unique_kmers.csv"
        species_unique.to_csv(species_unique_path, index=False)
        print(
            f"Species-unique k-mers saved to {species_unique_path} "
            f"({len(species_unique)} rows)"
        )
        return 0
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())