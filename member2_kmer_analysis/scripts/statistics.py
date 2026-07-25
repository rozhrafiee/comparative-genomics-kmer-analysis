#!/usr/bin/env python3
"""Statistical analysis of k-mer diversity, Shannon entropy, and (optionally)
GC content across organisms — بخش پنجم (تحلیل آماری) of the project.

Two families of tests are produced:

1. Species-level correlation (n = number of organisms, one pooled value each):
   - Pearson & Spearman correlation between genome size and Shannon entropy
   - Pearson & Spearman correlation between genome size and k-mer diversity
   - The same two correlations against GC% (if member1's GC results are
     found under ../member1_gc_analysis/results/, read-only, never modified)

2. Group-level significance tests (samples = individual chromosomes/
   scaffolds within each organism, so each species has multiple data
   points to compare against the others):
   - One-way ANOVA (parametric) across organisms
   - Kruskal-Wallis (non-parametric) across organisms
   - Pairwise Mann-Whitney U between every pair of organisms, with a
     Bonferroni-corrected significance flag

NOTE ON STATISTICAL POWER: with only ~5 organisms, the species-level
correlations (test family 1) have very low power and should be read as
descriptive/exploratory, not as strong statistical evidence. This is
mentioned again in the printed output and should be repeated in your
written report.
"""

from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path

import pandas as pd
from scipy import stats

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from kmer_analysis import (
    DEFAULT_INPUT_DIR,
    DEFAULT_K,
    DEFAULT_OUTPUT_DIR,
    count_kmers,
    load_sequences_by_record,
    resolve_input_dir,
)
from entropy import shannon_entropy

# Default (read-only) location of member1's GC summary — never written to.
DEFAULT_GC_SUMMARY_PATH = (
    DEFAULT_INPUT_DIR.parent
    / "member1_gc_analysis"
    / "results"
    / "tables"
    / "student1_gc_genome_summary.csv"
)

ALPHA = 0.05


def _norm_name(name: str) -> str:
    """Normalize an organism label for cross-module matching
    (member2 uses underscores, member1 uses spaces)."""
    return name.replace("_", " ").strip().lower()


def build_per_record_metrics(input_dir: Path, k: int) -> pd.DataFrame:
    """Compute Shannon entropy and k-mer diversity for every individual
    chromosome/scaffold, grouped by organism. This is the per-record
    (unpooled) table used as the sample groups for ANOVA / Kruskal-Wallis /
    Mann-Whitney."""
    grouped = load_sequences_by_record(input_dir)
    rows = []
    for organism, records in grouped.items():
        for record_id, seq in records:
            counts = count_kmers(seq, k)
            total = sum(counts.values())
            unique = len(counts)
            rows.append(
                {
                    "organism": organism,
                    "record_id": record_id,
                    "record_length": len(seq),
                    "k": k,
                    "shannon_entropy": shannon_entropy(dict(counts)),
                    "kmer_diversity": (unique / total) if total else 0.0,
                }
            )
    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError("No per-record metrics could be computed from the input data.")
    return df


def load_species_level_table(output_dir: Path) -> pd.DataFrame:
    """Merge the already-computed species-level kmer_summary.csv and
    entropy.csv into one table (one row per organism)."""
    kmer_path = output_dir / "kmer_summary.csv"
    entropy_path = output_dir / "entropy.csv"
    if not kmer_path.exists():
        raise FileNotFoundError(
            f"{kmer_path} not found. Run kmer_analysis.py first."
        )
    if not entropy_path.exists():
        raise FileNotFoundError(
            f"{entropy_path} not found. Run entropy.py first."
        )

    kmer_df = pd.read_csv(kmer_path)[
        ["sequence_id", "sequence_length", "kmer_diversity"]
    ]
    entropy_df = pd.read_csv(entropy_path)[["sequence_id", "shannon_entropy"]]
    merged = kmer_df.merge(entropy_df, on="sequence_id", how="inner")
    merged = merged.rename(columns={"sequence_id": "organism", "sequence_length": "genome_size"})
    return merged


def attach_gc_content(species_df: pd.DataFrame, gc_summary_path: Path) -> pd.DataFrame:
    """Optionally attach member1's GC% per organism (read-only). If the file
    isn't found, returns species_df unchanged with a printed note."""
    if not gc_summary_path.exists():
        print(
            f"Note: GC summary not found at {gc_summary_path} — "
            "skipping GC-content correlations. (member1_gc_analysis results "
            "are read-only and never modified by this script.)",
            file=sys.stderr,
        )
        return species_df

    try:
        gc_df = pd.read_csv(gc_summary_path)
    except (OSError, pd.errors.ParserError) as exc:
        print(f"Note: could not read GC summary ({exc}) — skipping GC correlations.",
              file=sys.stderr)
        return species_df

    if "Organism" not in gc_df.columns or "Genome_GC_Percent" not in gc_df.columns:
        print(
            "Note: GC summary is missing expected columns "
            "('Organism', 'Genome_GC_Percent') — skipping GC correlations.",
            file=sys.stderr,
        )
        return species_df

    gc_df = gc_df[["Organism", "Genome_GC_Percent"]].copy()
    gc_df["_match_key"] = gc_df["Organism"].map(_norm_name)
    species_df = species_df.copy()
    species_df["_match_key"] = species_df["organism"].map(_norm_name)

    merged = species_df.merge(
        gc_df[["_match_key", "Genome_GC_Percent"]], on="_match_key", how="left"
    ).drop(columns="_match_key")
    merged = merged.rename(columns={"Genome_GC_Percent": "gc_percent"})

    n_matched = merged["gc_percent"].notna().sum()
    print(f"Matched GC% for {n_matched}/{len(merged)} organism(s) from member1_gc_analysis.")
    return merged


def compute_correlations(species_df: pd.DataFrame) -> pd.DataFrame:
    """Pearson and Spearman correlation for each variable pair of interest,
    computed across organisms (n = number of organisms)."""
    pairs = [("genome_size", "shannon_entropy"), ("genome_size", "kmer_diversity")]
    if "gc_percent" in species_df.columns:
        pairs += [("gc_percent", "shannon_entropy"), ("gc_percent", "kmer_diversity")]

    rows = []
    for x_col, y_col in pairs:
        sub = species_df[[x_col, y_col]].dropna()
        n = len(sub)
        if n < 3:
            rows.append(
                {
                    "variable_x": x_col, "variable_y": y_col, "n": n,
                    "pearson_r": None, "pearson_p": None,
                    "spearman_rho": None, "spearman_p": None,
                    "note": "n < 3, correlation not computed",
                }
            )
            continue
        pearson_r, pearson_p = stats.pearsonr(sub[x_col], sub[y_col])
        spearman_rho, spearman_p = stats.spearmanr(sub[x_col], sub[y_col])
        rows.append(
            {
                "variable_x": x_col, "variable_y": y_col, "n": n,
                "pearson_r": pearson_r, "pearson_p": pearson_p,
                "spearman_rho": spearman_rho, "spearman_p": spearman_p,
                "note": "low statistical power (small n)" if n < 10 else "",
            }
        )
    return pd.DataFrame(rows)


def compute_group_tests(per_record_df: pd.DataFrame, metric: str) -> dict:
    """One-way ANOVA and Kruskal-Wallis across organism groups for a given
    per-record metric column ('shannon_entropy' or 'kmer_diversity')."""
    groups = [
        sub[metric].values
        for _, sub in per_record_df.groupby("organism")
        if len(sub) > 0
    ]
    group_sizes = per_record_df.groupby("organism").size().to_dict()
    usable_groups = [g for g in groups if len(g) >= 1]

    result = {
        "metric": metric,
        "n_groups": len(usable_groups),
        "group_sizes": ";".join(f"{k}={v}" for k, v in group_sizes.items()),
    }

    groups_with_variance = [g for g in usable_groups if len(g) >= 2]
    if len(usable_groups) >= 2 and all(len(g) >= 2 for g in usable_groups):
        f_stat, anova_p = stats.f_oneway(*usable_groups)
        h_stat, kw_p = stats.kruskal(*usable_groups)
        result.update(
            {
                "anova_f": f_stat, "anova_p": anova_p,
                "kruskal_h": h_stat, "kruskal_p": kw_p,
                "note": "",
            }
        )
    else:
        result.update(
            {
                "anova_f": None, "anova_p": None,
                "kruskal_h": None, "kruskal_p": None,
                "note": (
                    "At least one organism has < 2 chromosomes/scaffolds "
                    "in this dataset, so a valid group test could not be run."
                ),
            }
        )
    return result


def compute_pairwise_mannwhitney(per_record_df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Pairwise Mann-Whitney U test between every pair of organisms for a
    given per-record metric, with a Bonferroni-corrected significance flag."""
    organisms = sorted(per_record_df["organism"].unique())
    pairs = list(combinations(organisms, 2))
    n_comparisons = len(pairs) if pairs else 1
    bonferroni_alpha = ALPHA / n_comparisons

    rows = []
    for org_a, org_b in pairs:
        a = per_record_df.loc[per_record_df["organism"] == org_a, metric].values
        b = per_record_df.loc[per_record_df["organism"] == org_b, metric].values
        if len(a) < 1 or len(b) < 1:
            continue
        try:
            u_stat, p_value = stats.mannwhitneyu(a, b, alternative="two-sided")
        except ValueError as exc:
            rows.append(
                {
                    "organism_a": org_a, "organism_b": org_b, "metric": metric,
                    "n_a": len(a), "n_b": len(b),
                    "u_stat": None, "p_value": None,
                    "significant_raw": None, "significant_bonferroni": None,
                    "note": str(exc),
                }
            )
            continue
        rows.append(
            {
                "organism_a": org_a, "organism_b": org_b, "metric": metric,
                "n_a": len(a), "n_b": len(b),
                "u_stat": u_stat, "p_value": p_value,
                "significant_raw": bool(p_value < ALPHA),
                "significant_bonferroni": bool(p_value < bonferroni_alpha),
                "note": f"bonferroni_alpha={bonferroni_alpha:.4f} across {n_comparisons} comparisons",
            }
        )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Correlation and significance testing across organisms."
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("-k", type=int, default=DEFAULT_K)
    parser.add_argument(
        "--gc-summary-path",
        type=Path,
        default=DEFAULT_GC_SUMMARY_PATH,
        help="Read-only path to member1's genome GC summary CSV (optional).",
    )
    parser.add_argument("--demo", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        input_dir = resolve_input_dir(args.input_dir, demo=args.demo)
        args.output_dir.mkdir(parents=True, exist_ok=True)

        # --- Species-level correlations (n = number of organisms) ---
        species_df = load_species_level_table(args.output_dir)
        species_df = attach_gc_content(species_df, args.gc_summary_path)
        species_path = args.output_dir / "species_level_metrics.csv"
        species_df.to_csv(species_path, index=False)
        print(f"Species-level metrics table saved to {species_path}")

        correlation_df = compute_correlations(species_df)
        correlation_path = args.output_dir / "correlation_results.csv"
        correlation_df.to_csv(correlation_path, index=False)
        print(f"Correlation results saved to {correlation_path}")
        if len(species_df) < 10:
            print(
                f"Note: only {len(species_df)} organism(s) — species-level "
                "correlations have low statistical power and are exploratory only.",
                file=sys.stderr,
            )

        # --- Group-level tests (samples = chromosomes/scaffolds per organism) ---
        per_record_df = build_per_record_metrics(input_dir, args.k)
        per_record_path = args.output_dir / "per_record_metrics.csv"
        per_record_df.to_csv(per_record_path, index=False)
        print(f"Per-record metrics table saved to {per_record_path}")

        group_test_rows = [
            compute_group_tests(per_record_df, "shannon_entropy"),
            compute_group_tests(per_record_df, "kmer_diversity"),
        ]
        group_test_df = pd.DataFrame(group_test_rows)
        group_test_path = args.output_dir / "anova_kruskal_results.csv"
        group_test_df.to_csv(group_test_path, index=False)
        print(f"ANOVA / Kruskal-Wallis results saved to {group_test_path}")

        mw_entropy = compute_pairwise_mannwhitney(per_record_df, "shannon_entropy")
        mw_diversity = compute_pairwise_mannwhitney(per_record_df, "kmer_diversity")
        mw_df = pd.concat([mw_entropy, mw_diversity], ignore_index=True)
        mw_path = args.output_dir / "mannwhitney_pairwise.csv"
        mw_df.to_csv(mw_path, index=False)
        print(f"Pairwise Mann-Whitney results saved to {mw_path}")

        return 0
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())