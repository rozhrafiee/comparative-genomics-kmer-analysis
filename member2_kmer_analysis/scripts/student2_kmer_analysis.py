#!/usr/bin/env python3
"""Student 2: k-mer, similarity, statistics, and TATA-box analysis.

Implements exactly what the project brief asks for (بخش سوم, چهارم, پنجم,
ششم) — no additional metrics, tests, or figures beyond that scope.

One script: load FASTAs -> compute -> write results/tables/ + results/figures/
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import Counter
from itertools import combinations, product
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import cdist, squareform
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import StandardScaler

PROJECT_ROOT = Path(__file__).resolve().parents[2]
MODULE_DIR = Path(__file__).resolve().parents[1]
INPUT_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = MODULE_DIR / "results"
GC_SUMMARY_PATH = (
    PROJECT_ROOT
    / "member1_gc_analysis"
    / "results"
    / "tables"
    / "student1_gc_genome_summary.csv"
)

# بخش سوم explicitly asks for 3-mer, 4-mer, and 5-mer answers.
K_VALUES = (3, 4, 5)
# k=4 is used for similarity / PCA / clustering / correlation-matrix figures
# (the brief's example table and every downstream question refer to a
# single k-mer feature space, not three separate ones).
DEFAULT_K = 4
PROMOTER_LENGTH = 500  # بخش ششم: TATA box is described as ~25-35bp upstream
ALPHA = 0.05
DPI = 150

TATA_CANONICAL = re.compile(r"TATAAA")


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def normalize_sequence(sequence: str) -> str:
    return "".join(base for base in sequence.upper() if base in "ACGT")


def label(name: str) -> str:
    return str(name).replace("_", " ")


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    print(f"  wrote {path.name}")


def read_organism_name(dataset_dir: Path) -> str:
    summary_files = list(dataset_dir.rglob("data_summary.tsv"))
    if not summary_files:
        return dataset_dir.name
    try:
        with summary_files[0].open(encoding="utf-8") as handle:
            lines = [line.strip() for line in handle if line.strip()]
        if len(lines) >= 2:
            organism = lines[1].split("\t")[0].strip()
            parts = organism.split()
            if len(parts) >= 2:
                return f"{parts[0]}_{parts[1]}"
            return organism.replace(" ", "_")
    except OSError:
        pass
    return dataset_dir.name


def select_genomic_fasta(dataset_dir: Path) -> Path | None:
    fna_files = sorted(dataset_dir.rglob("*_genomic.fna"))
    gcf = [p for p in fna_files if p.name.startswith("GCF_")]
    gca = [p for p in fna_files if p.name.startswith("GCA_")]
    if gcf:
        return gcf[0]
    if gca:
        return gca[0]
    return None


def discover_genomes(input_dir: Path) -> dict[str, Path]:
    genomes: dict[str, Path] = {}
    for dataset_dir in sorted(input_dir.glob("ncbi_dataset*")):
        if not dataset_dir.is_dir():
            continue
        fasta = select_genomic_fasta(dataset_dir)
        if fasta is None:
            continue
        org = read_organism_name(dataset_dir)
        key, n = org, 2
        while key in genomes:
            key = f"{org}_{n}"
            n += 1
        genomes[key] = fasta
    return genomes


def load_records(input_dir: Path) -> dict[str, list[tuple[str, str]]]:
    """{organism: [(record_id, ACGT sequence), ...]}"""
    genomes = discover_genomes(input_dir)
    if not genomes:
        raise FileNotFoundError(f"No ncbi_dataset* genomic FASTA under {input_dir}")
    grouped: dict[str, list[tuple[str, str]]] = {}
    for organism, fasta_path in genomes.items():
        print(f"Loading {label(organism)} from {fasta_path.name} ...")
        records: list[tuple[str, str]] = []
        for rec in SeqIO.parse(fasta_path, "fasta"):
            seq = normalize_sequence(str(rec.seq))
            if seq:
                records.append((rec.id, seq))
        if not records:
            raise ValueError(f"No sequences in {fasta_path}")
        grouped[organism] = records
    return grouped


def pool_sequences(grouped: dict[str, list[tuple[str, str]]]) -> dict[str, str]:
    return {org: "".join(seq for _, seq in recs) for org, recs in grouped.items()}


# ---------------------------------------------------------------------------
# K-mer core
# ---------------------------------------------------------------------------


def get_all_kmers(k: int) -> list[str]:
    return ["".join(bases) for bases in product("ACGT", repeat=k)]


def count_kmers(sequence: str, k: int) -> Counter:
    counts: Counter = Counter()
    if len(sequence) < k:
        return counts
    for i in range(len(sequence) - k + 1):
        counts[sequence[i : i + k]] += 1
    return counts


def shannon_entropy(counts: Counter) -> float:
    total = sum(counts.values())
    if total == 0:
        return 0.0
    h = 0.0
    for c in counts.values():
        if c > 0:
            p = c / total
            h -= p * math.log2(p)
    return h


def frequency_vector(counts: Counter, k: int) -> dict[str, float]:
    """Normalized k-mer frequencies (count / total) — the numeric vector
    representation each organism is reduced to for similarity / PCA."""
    total = sum(counts.values())
    kmers = get_all_kmers(k)
    if total == 0:
        return {km: 0.0 for km in kmers}
    return {km: counts.get(km, 0) / total for km in kmers}


def build_feature_matrix(organism_counts: dict[str, Counter], k: int) -> pd.DataFrame:
    rows = []
    for org, counts in organism_counts.items():
        rows.append({"organism": org, **frequency_vector(counts, k)})
    return pd.DataFrame(rows).set_index("organism")


# ---------------------------------------------------------------------------
# بخش سوم: K-mer questions
# ---------------------------------------------------------------------------


def build_most_common_table(organism_counts: dict[str, Counter], k: int) -> pd.DataFrame:
    """بخش سوم Q1-3: the single most frequent k-mer per organism, for
    k=3, k=4, k=5."""
    rows = []
    for org, counts in organism_counts.items():
        if counts:
            kmer, count = counts.most_common(1)[0]
        else:
            kmer, count = None, 0
        rows.append({"k": k, "organism": org, "kmer": kmer, "count": count})
    return pd.DataFrame(rows)


def build_diversity_entropy(organism_counts: dict[str, Counter], k: int) -> pd.DataFrame:
    """بخش سوم Q4 (unique k-mers), Q6 (diversity ranking) — plus Shannon
    entropy, needed by بخش پنجم's significance tests."""
    rows = []
    for org, counts in organism_counts.items():
        total = sum(counts.values())
        unique = len(counts)
        rows.append(
            {
                "k": k,
                "organism": org,
                "unique_kmers": unique,
                "total_kmers": total,
                "kmer_diversity": (unique / total) if total else 0.0,
                "shannon_entropy": shannon_entropy(counts),
            }
        )
    return pd.DataFrame(rows)


def build_unique_kmers(organism_counts: dict[str, Counter], k: int) -> pd.DataFrame:
    """بخش سوم Q5: k-mers observed in exactly one organism."""
    sets = {org: set(counts) for org, counts in organism_counts.items()}
    rows = []
    for org, kmers in sets.items():
        others = set().union(*(s for o, s in sets.items() if o != org))
        for kmer in sorted(kmers - others):
            rows.append(
                {
                    "k": k,
                    "kmer": kmer,
                    "exclusive_organism": org,
                    "count": organism_counts[org][kmer],
                }
            )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# بخش چهارم: similarity / PCA / clustering
# ---------------------------------------------------------------------------


def distance_matrix(feat: pd.DataFrame) -> pd.DataFrame:
    d = cdist(feat.values, feat.values, metric="euclidean")
    return pd.DataFrame(d, index=feat.index, columns=feat.index)


def cosine_matrix(feat: pd.DataFrame) -> pd.DataFrame:
    s = cosine_similarity(feat.values)
    return pd.DataFrame(s, index=feat.index, columns=feat.index)


def run_pca(feat: pd.DataFrame) -> pd.DataFrame:
    n_comp = min(2, feat.shape[0], feat.shape[1])
    scaled = StandardScaler().fit_transform(feat.values)
    pca = PCA(n_components=n_comp)
    coords = pca.fit_transform(scaled)
    cols = [f"PC{i + 1}" for i in range(n_comp)]
    out = pd.DataFrame(coords, index=feat.index, columns=cols)
    out.index.name = "organism"
    for i, col in enumerate(cols):
        out[f"{col}_variance_ratio"] = pca.explained_variance_ratio_[i]
    return out.reset_index()


# ---------------------------------------------------------------------------
# بخش ششم: TATA box
# ---------------------------------------------------------------------------


def analyze_tata(sequences: dict[str, str]) -> pd.DataFrame:
    """Scans the first PROMOTER_LENGTH bases of each pooled genome for the
    canonical TATAAA motif. NOTE: this is a raw-sequence scan, not a
    gene/promoter-annotation-based scan (no GFF/GTF used), so it cannot
    answer "% of genes with a TATA box" — only whether the motif is present
    near the start of the assembled sequence."""
    rows = []
    for org, seq in sequences.items():
        promoter = seq[:PROMOTER_LENGTH]
        canonical = list(TATA_CANONICAL.finditer(promoter))
        rows.append(
            {
                "organism": org,
                "promoter_length": PROMOTER_LENGTH,
                "tata_box_count": len(canonical),
                "has_canonical_tata": len(canonical) > 0,
                "first_tata_position": canonical[0].start() + 1 if canonical else None,
                "all_tata_positions": ";".join(str(m.start() + 1) for m in canonical),
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# بخش پنجم: statistics
# ---------------------------------------------------------------------------


def _norm_name(name: str) -> str:
    return name.replace("_", " ").strip().lower()


def attach_gc(species_df: pd.DataFrame, gc_path: Path) -> pd.DataFrame:
    if not gc_path.exists():
        print(f"Note: GC summary not found at {gc_path} — skipping GC columns.")
        return species_df
    gc = pd.read_csv(gc_path)
    if "Organism" not in gc.columns or "Genome_GC_Percent" not in gc.columns:
        print("Note: GC summary missing expected columns — skipping GC columns.")
        return species_df
    gc = gc[["Organism", "Genome_GC_Percent"]].copy()
    gc["_key"] = gc["Organism"].map(_norm_name)
    out = species_df.copy()
    out["_key"] = out["organism"].map(_norm_name)
    out = out.merge(gc[["_key", "Genome_GC_Percent"]], on="_key", how="left")
    out = out.drop(columns="_key").rename(columns={"Genome_GC_Percent": "gc_percent"})
    print(f"Matched GC% for {out['gc_percent'].notna().sum()}/{len(out)} organisms.")
    return out


def build_per_record_metrics(
    grouped: dict[str, list[tuple[str, str]]], k: int
) -> pd.DataFrame:
    """Per-chromosome/scaffold entropy and diversity, grouped by organism.
    Needed as the sample groups for ANOVA / Kruskal-Wallis / Mann-Whitney —
    a single pooled value per organism (n=1) cannot be used in a group
    significance test."""
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
                    "shannon_entropy": shannon_entropy(counts),
                    "kmer_diversity": (unique / total) if total else 0.0,
                }
            )
    return pd.DataFrame(rows)


def compute_statistical_results(
    species_df: pd.DataFrame, per_record_df: pd.DataFrame
) -> pd.DataFrame:
    """بخش پنجم: exactly the 3 named correlation pairs + TATA-vs-GC
    (بخش ششم Q5), plus ANOVA / Kruskal-Wallis / Mann-Whitney significance
    tests. Permutation Test is listed in the brief as optional
    ("می‌توان استفاده کرد") and is not included here."""
    rows: list[dict] = []

    # The 3 pairs named in بخش پنجم, plus TATA↔GC named in بخش ششم Q5.
    pairs = [
        ("genome_size", "kmer_diversity"),
    ]
    if "gc_percent" in species_df.columns:
        pairs += [
            ("gc_percent", "genome_size"),
            ("gc_percent", "shannon_entropy"),
        ]
    if "gc_percent" in species_df.columns and "tata_box_count" in species_df.columns:
        pairs.append(("tata_box_count", "gc_percent"))

    for x_col, y_col in pairs:
        sub = species_df[[x_col, y_col]].dropna()
        n = len(sub)
        if n < 3:
            rows.append(
                {
                    "analysis": "correlation",
                    "test": "Pearson/Spearman",
                    "variables": f"{x_col} vs {y_col}",
                    "n": n,
                    "statistic": None,
                    "p_value": None,
                    "extra": None,
                    "note": "n < 3",
                }
            )
            continue
        x_vals, y_vals = sub[x_col].to_numpy(), sub[y_col].to_numpy()
        pr, pp = stats.pearsonr(x_vals, y_vals)
        sr, sp = stats.spearmanr(x_vals, y_vals)
        rows.append(
            {
                "analysis": "correlation",
                "test": "Pearson",
                "variables": f"{x_col} vs {y_col}",
                "n": n,
                "statistic": pr,
                "p_value": pp,
                "extra": f"spearman_rho={sr:.6f};spearman_p={sp:.6g}",
                "note": "low power (small n)" if n < 10 else "",
            }
        )

    for metric in ("shannon_entropy", "kmer_diversity"):
        groups = [
            g[metric].to_numpy()
            for _, g in per_record_df.groupby("organism")
            if len(g) >= 2
        ]
        sizes = ";".join(
            f"{org}={n}" for org, n in per_record_df.groupby("organism").size().items()
        )
        if len(groups) < 2:
            rows.append(
                {
                    "analysis": "group_test",
                    "test": "ANOVA/Kruskal",
                    "variables": metric,
                    "n": len(per_record_df),
                    "statistic": None,
                    "p_value": None,
                    "extra": sizes,
                    "note": "insufficient groups",
                }
            )
            continue
        f_stat, anova_p = stats.f_oneway(*groups)
        h_stat, kw_p = stats.kruskal(*groups)
        rows.append(
            {
                "analysis": "group_test",
                "test": "ANOVA",
                "variables": metric,
                "n": len(per_record_df),
                "statistic": f_stat,
                "p_value": anova_p,
                "extra": f"kruskal_h={h_stat:.6f};kruskal_p={kw_p:.6g};{sizes}",
                "note": "",
            }
        )

        organisms = sorted(per_record_df["organism"].unique())
        for a, b in combinations(organisms, 2):
            va = per_record_df.loc[per_record_df["organism"] == a, metric].to_numpy()
            vb = per_record_df.loc[per_record_df["organism"] == b, metric].to_numpy()
            if len(va) < 1 or len(vb) < 1:
                continue
            u, p = stats.mannwhitneyu(va, vb, alternative="two-sided")
            rows.append(
                {
                    "analysis": "pairwise",
                    "test": "Mann-Whitney",
                    "variables": f"{metric}: {a} vs {b}",
                    "n": len(va) + len(vb),
                    "statistic": u,
                    "p_value": p,
                    "extra": None,
                    "note": "",
                }
            )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figures — exactly the 5 chart types named in the brief:
# Heatmap (k-mer frequency), PCA, Hierarchical Clustering, Dendrogram,
# Correlation Matrix. (Hierarchical Clustering itself is the linkage
# computation used to draw the Dendrogram — no separate figure needed
# for it beyond the dendrogram.)
# ---------------------------------------------------------------------------


def save_fig(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure {path.name}")


def generate_figures(
    feat: pd.DataFrame,
    dist: pd.DataFrame,
    pca_df: pd.DataFrame,
    species_df: pd.DataFrame,
    figure_dir: Path,
) -> None:
    # 1. Heatmap — k-mer frequency across species (بخش سوم / توضیحات اضافه)
    fig, ax = plt.subplots(figsize=(12, 4))
    sns.heatmap(feat.T, cmap="viridis", ax=ax, cbar_kws={"label": "frequency"})
    ax.set_title(f"K-mer frequency heatmap (k={DEFAULT_K})")
    ax.set_xlabel("Organism")
    ax.set_ylabel("K-mer")
    ax.set_xticklabels([label(t.get_text()) for t in ax.get_xticklabels()], rotation=25, ha="right")
    ax.set_yticks([])
    save_fig(fig, figure_dir / "kmer_frequency_heatmap.png")

    # 2. PCA (بخش چهارم / توضیحات اضافه)
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for _, row in pca_df.iterrows():
        ax.scatter(row["PC1"], row["PC2"], s=80)
        ax.annotate(label(row["organism"]), (row["PC1"], row["PC2"]), fontsize=9)
    v1 = pca_df["PC1_variance_ratio"].iloc[0] * 100
    v2 = pca_df["PC2_variance_ratio"].iloc[0] * 100 if "PC2" in pca_df.columns else 0
    ax.set_xlabel(f"PC1 ({v1:.1f}%)")
    ax.set_ylabel(f"PC2 ({v2:.1f}%)")
    ax.set_title(f"PCA of k-mer frequency profiles (k={DEFAULT_K})")
    save_fig(fig, figure_dir / "pca.png")

    # 3/4. Hierarchical Clustering + Dendrogram (بخش چهارم / توضیحات اضافه)
    condensed = squareform(dist.values, checks=False)
    z = linkage(condensed, method="ward")
    fig, ax = plt.subplots(figsize=(8, 5))
    dendrogram(z, labels=[label(i) for i in dist.index], leaf_rotation=25, ax=ax)
    ax.set_title(f"Hierarchical clustering (Ward, Euclidean k={DEFAULT_K})")
    ax.set_ylabel("Distance")
    save_fig(fig, figure_dir / "dendrogram.png")

    # 5. Correlation Matrix (بخش پنجم / توضیحات اضافه)
    cols = [c for c in ["genome_size", "gc_percent", "shannon_entropy", "kmer_diversity"] if c in species_df.columns]
    if len(cols) >= 2:
        corr = species_df[cols].corr(method="pearson")
        fig, ax = plt.subplots(figsize=(6, 5))
        sns.heatmap(corr, annot=True, fmt=".2f", cmap="coolwarm", vmin=-1, vmax=1, ax=ax)
        ax.set_title("Species-level Pearson correlation matrix")
        save_fig(fig, figure_dir / "correlation_matrix.png")


# ---------------------------------------------------------------------------
# Brief answers — printed AND written to a text file so they can be
# quoted directly in the report (a printout alone isn't citable).
# ---------------------------------------------------------------------------


def write_brief_answers(
    top_df: pd.DataFrame,
    unique_df: pd.DataFrame,
    div_ent: pd.DataFrame,
    dist: pd.DataFrame,
    tata_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    path: Path,
) -> None:
    lines = ["# Brief sections 3-6 — direct answers\n"]

    lines.append("## بخش سوم — K-mer\n")
    for k in K_VALUES:
        sub = top_df[top_df["k"] == k]
        lines.append(f"Most frequent {k}-mer per organism:")
        for _, r in sub.iterrows():
            lines.append(f"- {label(r['organism'])}: {r['kmer']} (n={int(r['count'])})")

    d4 = div_ent[div_ent["k"] == DEFAULT_K]
    lines.append(f"\nUnique k-mer counts (k={DEFAULT_K}):")
    for _, r in d4.iterrows():
        lines.append(f"- {label(r['organism'])}: {int(r['unique_kmers'])}")

    if unique_df.empty:
        lines.append("\nSpecies-exclusive k-mers: none found")
    else:
        lines.append(f"\nSpecies-exclusive k-mers: {len(unique_df)} total (see student2_unique_kmers.csv)")

    best = d4.loc[d4["kmer_diversity"].idxmax()]
    lines.append(
        f"\nHighest k-mer diversity (k={DEFAULT_K}): {label(best['organism'])} "
        f"({best['kmer_diversity']:.6g})"
    )

    lines.append("\n## بخش چهارم — Similarity\n")
    vals = dist.values.copy()
    np.fill_diagonal(vals, np.inf)
    i, j = np.unravel_index(np.argmin(vals), vals.shape)
    lines.append(f"Closest species (Euclidean): {label(dist.index[i])} - {label(dist.index[j])} ({vals[i, j]:.6f})")
    vals2 = dist.values.copy()
    np.fill_diagonal(vals2, -np.inf)
    i, j = np.unravel_index(np.argmax(vals2), vals2.shape)
    lines.append(f"Farthest species (Euclidean): {label(dist.index[i])} - {label(dist.index[j])} ({vals2[i, j]:.6f})")
    lines.append(
        "Animal clustering: cannot be tested — only one animal organism "
        "(C. elegans) is in this dataset."
    )
    lines.append(
        "Plant vs animal separation: cannot be tested — no plant organism "
        "is in this dataset."
    )

    lines.append("\n## بخش ششم — TATA Box\n")
    lines.append(f"Canonical TATAAA in first {PROMOTER_LENGTH} bp of pooled genome (not gene-annotation based):")
    for _, r in tata_df.iterrows():
        lines.append(f"- {label(r['organism'])}: {int(r['tata_box_count'])} match(es)")
    lines.append(
        "% of genes with a TATA box: cannot be answered — no gene/promoter "
        "annotation (GFF/GTF) was used, only a raw-sequence scan."
    )
    lines.append(
        "TATA+ genes having higher expression: cannot be answered — no "
        "expression data was collected in this project."
    )

    lines.append("\n## بخش پنجم — Statistics\n")
    corr = stats_df[stats_df["analysis"] == "correlation"]
    for _, r in corr.iterrows():
        lines.append(f"- Correlation {r['variables']}: r={r['statistic']}, p={r['p_value']} ({r['note']})")
    grp = stats_df[(stats_df["analysis"] == "group_test") & (stats_df["test"] == "ANOVA")]
    for _, r in grp.iterrows():
        lines.append(f"- ANOVA {r['variables']}: F={r['statistic']}, p={r['p_value']}")

    text = "\n".join(lines) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    print(text)
    print(f"  wrote {path.name}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument(
        "--gc-summary-path",
        type=Path,
        default=GC_SUMMARY_PATH,
        help="Read-only member1 GC summary CSV",
    )
    args = parser.parse_args()

    table_dir = args.output_dir / "tables"
    figure_dir = args.output_dir / "figures"
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    grouped = load_records(args.input_dir)
    sequences = pool_sequences(grouped)
    organisms = list(sequences.keys())
    print(f"Loaded {len(organisms)} organisms: {', '.join(label(o) for o in organisms)}")

    all_top: list[pd.DataFrame] = []
    all_unique: list[pd.DataFrame] = []
    all_div: list[pd.DataFrame] = []
    counts_by_k: dict[int, dict[str, Counter]] = {}

    for k in K_VALUES:
        print(f"\n=== Counting k={k} ===")
        organism_counts: dict[str, Counter] = {}
        for org, seq in sequences.items():
            print(f"  {label(org)} ...", flush=True)
            organism_counts[org] = count_kmers(seq, k)
        counts_by_k[k] = organism_counts
        all_top.append(build_most_common_table(organism_counts, k))
        all_unique.append(build_unique_kmers(organism_counts, k))
        all_div.append(build_diversity_entropy(organism_counts, k))

    top_df = pd.concat(all_top, ignore_index=True)
    unique_parts = [df for df in all_unique if not df.empty]
    unique_df = (
        pd.concat(unique_parts, ignore_index=True)
        if unique_parts
        else pd.DataFrame(columns=["k", "kmer", "exclusive_organism", "count"])
    )
    div_ent = pd.concat(all_div, ignore_index=True)

    print(f"\n=== Similarity / PCA (k={DEFAULT_K}) ===")
    feat = build_feature_matrix(counts_by_k[DEFAULT_K], DEFAULT_K)
    dist = distance_matrix(feat)
    cos = cosine_matrix(feat)
    pca_df = run_pca(feat)

    print("\n=== TATA box ===")
    tata_df = analyze_tata(sequences)

    print(f"\n=== Per-record metrics (for group tests, k={DEFAULT_K}) ===")
    per_record_df = build_per_record_metrics(grouped, DEFAULT_K)

    d4 = div_ent[div_ent["k"] == DEFAULT_K][
        ["organism", "total_kmers", "kmer_diversity", "shannon_entropy"]
    ].copy()
    d4["genome_size"] = d4["total_kmers"] + DEFAULT_K - 1
    species_df = d4[["organism", "genome_size", "kmer_diversity", "shannon_entropy"]]
    species_df = attach_gc(species_df, args.gc_summary_path)
    species_df = species_df.merge(
        tata_df[["organism", "tata_box_count"]], on="organism", how="left"
    )

    print("\n=== Statistical tests ===")
    stats_df = compute_statistical_results(species_df, per_record_df)

    print("\n=== Writing tables ===")
    write_csv(top_df, table_dir / "student2_kmer_most_common.csv")
    write_csv(div_ent, table_dir / "student2_kmer_diversity_entropy.csv")
    write_csv(unique_df, table_dir / "student2_unique_kmers.csv")
    dist.to_csv(table_dir / "student2_distance_matrix.csv")
    print("  wrote student2_distance_matrix.csv")
    cos.to_csv(table_dir / "student2_cosine_similarity.csv")
    print("  wrote student2_cosine_similarity.csv")
    write_csv(pca_df, table_dir / "student2_pca.csv")
    write_csv(tata_df, table_dir / "student2_tata_box.csv")
    write_csv(species_df, table_dir / "student2_species_metrics.csv")
    write_csv(stats_df, table_dir / "student2_statistical_results.csv")

    print("\n=== Writing figures ===")
    generate_figures(feat, dist, pca_df, species_df, figure_dir)

    write_brief_answers(
        top_df, unique_df, div_ent, dist, tata_df, stats_df,
        table_dir / "student2_brief_answers.md",
    )
    print(f"\nDone. Tables -> {table_dir}")
    print(f"Figures -> {figure_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)