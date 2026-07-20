#!/usr/bin/env python3
"""Reproducible Student 1 genome-size and GC-content analysis.

The program intentionally limits its scope to FASTA quality control, genome
characteristics, GC content, GC-related statistics, and visualisations.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import logging
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent

INPUT_DIR = PROJECT_ROOT / "data" / "raw"
OUTPUT_DIR = PROJECT_ROOT / "results"

try:
    from statsmodels.stats.multitest import multipletests
except Exception:  # pragma: no cover - optional at runtime
    multipletests = None


EXPECTED_SPECIES = {
    "Caenorhabditis elegans",
    "Arabidopsis thaliana",
    "Schizosaccharomyces pombe",
    "Drosophila melanogaster",
    "Anopheles gambiae",
}

IUPAC_AMBIGUOUS = set("RYSWKMBDHV")
VALID_SYMBOLS = set("ACGTN") | IUPAC_AMBIGUOUS
ORGANISM_ORDER = [
    "Caenorhabditis elegans",
    "Arabidopsis thaliana",
    "Schizosaccharomyces pombe",
    "Drosophila melanogaster",
    "Anopheles gambiae",
    "Neurospora crassa",
    "Yarrowia lipolytica",
    "Dictyostelium discoideum",
    "Eremothecium coryli",
    "Unknown organism",
]
COLORS = {
    name: color
    for name, color in zip(
        ORGANISM_ORDER,
        ["#4472C4", "#70AD47", "#ED7D31", "#A5A5A5", "#FFC000", "#5B9BD5", "#8064A2", "#C0504D", "#9BBB59", "#7F7F7F"],
    )
}


def configure_logging(log_path: Path) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("student1_gc")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(logging.StreamHandler(sys.stdout))
    return logger


def discover_fasta_files(input_dir: Path) -> List[Path]:
    """Return uploaded FASTA nucleotide files in deterministic order."""
    return sorted([*input_dir.glob("*.fna"), *input_dir.glob("*.fasta"), *input_dir.glob("*.fa")])


def extract_scientific_name(description: str) -> str:
    """Extract a binomial from a FASTA description without guessing."""
    words = re.findall(r"[A-Za-z]+", description)
    if len(words) >= 2:
        return f"{words[0]} {words[1]}"
    return "Unknown organism"


def parse_assembly_metadata(path: Path) -> Tuple[str, str, str]:
    """Extract accession and assembly label from a standard uploaded filename."""
    match = re.match(r"^(GCA_\d+\.\d+)_([^\s]+?)(?:_genomic)?\.f(?:n?a|asta)$", path.name, re.I)
    if not match:
        return "Not available in the provided file", "Not available in the provided file", "Not available in the provided file"
    accession, remainder = match.groups()
    return "Not available in the provided file", remainder, accession


def classify_sequence_record(sequence_id: str, description: str) -> str:
    """Classify a record using explicit header evidence only."""
    text = f"{sequence_id} {description}".lower()
    if re.search(r"mitochond|mitochondrial", text):
        return "mitochondrion"
    if re.search(r"chloroplast|plastid", text):
        return "chloroplast"
    if re.search(r"plasmid", text):
        return "plasmid"
    if re.search(r"unlocalized", text):
        return "unlocalized"
    if re.search(r"unplaced|chrun|supercont", text):
        return "unplaced"
    if re.search(r"scaffold", text):
        return "scaffold"
    if re.search(r"contig", text):
        return "contig"
    if re.search(r"chromosome|linkage group|chromosome", text):
        return "chromosome"
    return "unknown"


def iter_fasta_records(path: Path, logger: logging.Logger) -> Iterable[Dict[str, Any]]:
    """Stream one summary dictionary per FASTA record."""
    counts = Counter()
    sequence_id: Optional[str] = None
    description = ""
    saw_header = False
    line_number = 0
    duplicate_ids: set[str] = set()
    seen_ids: set[str] = set()

    def emit() -> Optional[Dict[str, Any]]:
        nonlocal counts, sequence_id, description
        if sequence_id is None:
            return None
        if sequence_id in seen_ids:
            duplicate_ids.add(sequence_id)
        seen_ids.add(sequence_id)
        total = int(sum(counts.values()))
        valid = int(sum(counts[s] for s in "ATGC"))
        gc = int(counts["G"] + counts["C"])
        organism = extract_scientific_name(description)
        return {
            "Sequence_ID": sequence_id,
            "Sequence_Description": description,
            "Scientific_Name": organism,
            "Sequence_Type": classify_sequence_record(sequence_id, description),
            "Sequence_Length": total,
            "A_Count": int(counts["A"]),
            "T_Count": int(counts["T"]),
            "G_Count": int(counts["G"]),
            "C_Count": int(counts["C"]),
            "N_Count": int(counts["N"]),
            "Other_Ambiguous_Count": int(sum(counts[s] for s in IUPAC_AMBIGUOUS)),
            "Unexpected_Symbol_Count": int(counts["OTHER"]),
            "Valid_ATGC_Length": valid,
            "GC_Count": gc,
            "GC_Percent_Valid_ATGC": (100.0 * gc / valid) if valid else np.nan,
            "GC_Percent_All_Bases": (100.0 * gc / total) if total else np.nan,
            "Duplicate_Sequence_ID": sequence_id in duplicate_ids,
        }

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line_number += 1
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                old = emit()
                if old is not None:
                    yield old
                header = line[1:].strip()
                if not header:
                    logger.warning("%s line %d has an empty FASTA header", path.name, line_number)
                    sequence_id, description = "", ""
                else:
                    parts = header.split(None, 1)
                    sequence_id = parts[0]
                    description = parts[1] if len(parts) == 2 else ""
                counts = Counter()
                saw_header = True
                continue
            if not saw_header:
                logger.warning("%s line %d contains sequence before first FASTA header", path.name, line_number)
                continue
            for char in line.upper():
                if char in "ACGTN" or char in IUPAC_AMBIGUOUS:
                    counts[char] += 1
                else:
                    counts["OTHER"] += 1
    final = emit()
    if final is not None:
        yield final
    if duplicate_ids:
        logger.warning("%s has duplicate sequence IDs: %s", path.name, ", ".join(sorted(duplicate_ids)[:10]))


def validate_fasta_file(path: Path, logger: logging.Logger) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Read and validate one FASTA file using streaming counting."""
    records = list(iter_fasta_records(path, logger))
    file_bytes = path.stat().st_size
    descriptions = [r["Sequence_Description"] for r in records]
    species = sorted(set(r["Scientific_Name"] for r in records if r["Scientific_Name"] != "Unknown organism"))
    source, assembly, accession = parse_assembly_metadata(path)
    summary = {
        "Input_File": path.name,
        "File_Size_Bytes": file_bytes,
        "File_Size_MB": file_bytes / 1_000_000,
        "Number_of_FASTA_Records": len(records),
        "FASTA_Record_Identifiers": " | ".join(r["Sequence_ID"] for r in records),
        "FASTA_Descriptions": " | ".join(descriptions),
        "Scientific_Names_In_Headers": " | ".join(species) if species else "Not available in the provided file",
        "Possible_Organism": species[0] if len(species) == 1 else "; ".join(species) if species else "Unknown organism",
        "Possible_Accession": accession,
        "Possible_Assembly_Version": assembly,
        "Possible_Source": source,
        "Header_Name_Consistency": "PASS" if len(species) <= 1 else "WARNING: multiple names",
        "Empty_File": "PASS" if file_bytes > 0 else "FAIL",
        "Zero_Length_Records": int(sum(r["Sequence_Length"] == 0 for r in records)),
        "Duplicate_Sequence_IDs": int(sum(r["Duplicate_Sequence_ID"] for r in records)),
        "Unexpected_Symbol_Count": int(sum(r["Unexpected_Symbol_Count"] for r in records)),
        "Lowercase_or_mixed_case_processed": "YES (case-insensitive)",
    }
    for record in records:
        record.update({"Input_File": path.name, "File_Size_Bytes": file_bytes, "File_Size_MB": file_bytes / 1_000_000})
    return records, summary


def ordered_names(values: Sequence[str]) -> List[str]:
    order = {name: i for i, name in enumerate(ORGANISM_ORDER)}
    return sorted(values, key=lambda x: (order.get(x, len(order)), x))


def infer_genome_rows(records: pd.DataFrame, file_audit: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for (organism, input_file), group in records.groupby(["Scientific_Name", "Input_File"], sort=False):
        audit = file_audit[file_audit["Input_File"] == input_file].iloc[0]
        source, assembly, accession = parse_assembly_metadata(Path(input_file))
        seq_types = group["Sequence_Type"]
        a = int(group["A_Count"].sum())
        t = int(group["T_Count"].sum())
        g = int(group["G_Count"].sum())
        c = int(group["C_Count"].sum())
        n = int(group["N_Count"].sum())
        other = int(group["Other_Ambiguous_Count"].sum())
        total = int(group["Sequence_Length"].sum())
        valid = int(group["Valid_ATGC_Length"].sum())
        gc = int(group["GC_Count"].sum())
        expected_match = organism in EXPECTED_SPECIES
        rows.append({
            "Organism": organism,
            "Scientific_Name": organism,
            "Input_File": input_file,
            "Expected_Species_Match": "YES" if expected_match else "NO - available file is a substitution",
            "Source": source,
            "Assembly_Version": assembly,
            "Accession": accession,
            "Number_of_FASTA_Records": len(group),
            "Number_of_Chromosomes": int((seq_types == "chromosome").sum()),
            "Number_of_Scaffolds": int((seq_types == "scaffold").sum()),
            "Number_of_Contigs": int((seq_types == "contig").sum()),
            "Number_of_Organellar_Records": int(seq_types.isin(["mitochondrion", "chloroplast"]).sum()),
            "Number_of_Unplaced_Records": int((seq_types == "unplaced").sum()),
            "Number_of_Unlocalized_Records": int((seq_types == "unlocalized").sum()),
            "Total_Sequence_Length": total,
            "Valid_ATGC_Length": valid,
            "Ambiguous_Base_Count": n + other,
            "Genome_GC_Count": gc,
            "Genome_GC_Percent": 100.0 * gc / valid if valid else np.nan,
            "Genome_GC_Percent_All_Bases": 100.0 * gc / total if total else np.nan,
            "File_Size_Bytes": int(audit["File_Size_Bytes"]),
            "File_Size_MB": float(audit["File_Size_MB"]),
            "A_Count": a,
            "T_Count": t,
            "G_Count": g,
            "C_Count": c,
            "N_Count": n,
            "Other_Ambiguous_Count": other,
            "Unexpected_Symbol_Count": int(group["Unexpected_Symbol_Count"].sum()),
            "Ambiguous_Percent": 100.0 * (n + other) / total if total else np.nan,
            "N_Percent": 100.0 * n / total if total else np.nan,
        })
    result = pd.DataFrame(rows)
    result["Organism_Order"] = result["Organism"].map({name: i for i, name in enumerate(ORGANISM_ORDER)}).fillna(999)
    return result.sort_values(["Organism_Order", "Organism"]).drop(columns="Organism_Order")


def calculate_gc_results(records: pd.DataFrame, genome: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for organism in ordered_names(genome["Organism"].tolist()):
        group = records[records["Scientific_Name"] == organism]
        all_gc = group["GC_Percent_Valid_ATGC"].dropna()
        main = group[group["Sequence_Type"] == "chromosome"]["GC_Percent_Valid_ATGC"].dropna()
        nuclear = group[~group["Sequence_Type"].isin(["mitochondrion", "chloroplast"])]
        for analysis_name, subset in [("All_FASTA_records", group), ("Main_chromosomes_only", group[group["Sequence_Type"] == "chromosome"]), ("Nuclear_genome_only", nuclear), ("Organellar_records", group[group["Sequence_Type"].isin(["mitochondrion", "chloroplast"])])]:
            vals = subset["GC_Percent_Valid_ATGC"].dropna()
            if vals.empty:
                continue
            lengths = subset.loc[vals.index, "Valid_ATGC_Length"]
            weighted = 100.0 * subset.loc[vals.index, "GC_Count"].sum() / lengths.sum() if lengths.sum() else np.nan
            q1, q3 = vals.quantile([0.25, 0.75]) if len(vals) > 1 else (np.nan, np.nan)
            mean = vals.mean()
            rows.append({
                "Organism": organism,
                "Analysis": analysis_name,
                "Record_Count": int(len(vals)),
                "Weighted_GC_Percent": weighted,
                "Unweighted_Mean_Record_GC_Percent": mean,
                "Median_Record_GC_Percent": vals.median(),
                "Standard_Deviation": vals.std(ddof=1) if len(vals) > 1 else 0.0,
                "Minimum": vals.min(),
                "Maximum": vals.max(),
                "Range": vals.max() - vals.min(),
                "Interquartile_Range": q3 - q1 if len(vals) > 1 else 0.0,
                "Coefficient_of_Variation_Percent": 100.0 * vals.std(ddof=1) / mean if len(vals) > 1 and mean else 0.0,
                "Ambiguous_Percent": float(subset["Other_Ambiguous_Count"].sum() + subset["N_Count"].sum()) / subset["Sequence_Length"].sum() * 100.0 if subset["Sequence_Length"].sum() else np.nan,
                "N_Percent": float(subset["N_Count"].sum()) / subset["Sequence_Length"].sum() * 100.0 if subset["Sequence_Length"].sum() else np.nan,
            })
    return pd.DataFrame(rows)


def bootstrap_correlation(x: np.ndarray, y: np.ndarray, method: str, seed: int = 42, n_boot: int = 5000) -> Tuple[float, float]:
    rng = np.random.default_rng(seed)
    values: List[float] = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(x), len(x))
        if len(np.unique(idx)) < 2:
            continue
        try:
            value = stats.pearsonr(x[idx], y[idx]).statistic if method == "Pearson" else stats.spearmanr(x[idx], y[idx]).statistic
            if np.isfinite(value):
                values.append(float(value))
        except Exception:
            continue
    if not values:
        return np.nan, np.nan
    return float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))


def run_statistical_tests(records: pd.DataFrame, genome: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    numeric = genome.dropna(subset=["Total_Sequence_Length", "Genome_GC_Percent"])
    x = numeric["Total_Sequence_Length"].to_numpy(float)
    y = numeric["Genome_GC_Percent"].to_numpy(float)
    if len(x) >= 3:
        pear = stats.pearsonr(x, y)
        spear = stats.spearmanr(x, y)
        for test_name, statistic, p_value, method in [("Pearson correlation", float(pear.statistic), float(pear.pvalue), "Pearson"), ("Spearman correlation", float(spear.statistic), float(spear.pvalue), "Spearman")]:
            low, high = bootstrap_correlation(x, y, method)
            rows.append({"Analysis": "Organism-level genome size vs GC", "Test": test_name, "Variables": "Total_Sequence_Length; Genome_GC_Percent", "Grouping_Variable": "Organism", "Sample_Size": len(x), "Statistic": statistic, "P_Value": p_value, "Adjusted_P_Value": np.nan, "Effect_Size": statistic, "Confidence_Interval": f"[{low:.6f}, {high:.6f}]" if np.isfinite(low) else "Not estimable", "Assumptions": "n=5; exploratory correlation; bootstrap percentile interval", "Interpretation": "Direction and magnitude are descriptive; p-values have low power with five genomes.", "Limitations": "Five organisms; correlation does not imply causation."})

    def add_row(analysis: str, test: str, variables: str, grouping: str, n: int, statistic: float, p: float, effect: Any, assumptions: str, interpretation: str, limitations: str, adjusted: Any = np.nan, ci: str = "Not estimable") -> None:
        rows.append({"Analysis": analysis, "Test": test, "Variables": variables, "Grouping_Variable": grouping, "Sample_Size": n, "Statistic": statistic, "P_Value": p, "Adjusted_P_Value": adjusted, "Effect_Size": effect, "Confidence_Interval": ci, "Assumptions": assumptions, "Interpretation": interpretation, "Limitations": limitations})

    for analysis_name, subset in [("All FASTA records", records), ("Main chromosomes only", records[records["Sequence_Type"] == "chromosome"]), ("Nuclear records only", records[~records["Sequence_Type"].isin(["mitochondrion", "chloroplast"])])]:
        subset = subset.dropna(subset=["GC_Percent_Valid_ATGC"])
        groups = [g["GC_Percent_Valid_ATGC"].to_numpy(float) for _, g in subset.groupby("Scientific_Name") if len(g) >= 1]
        names = [n for n, g in subset.groupby("Scientific_Name") if len(g) >= 1]
        if len(groups) < 2 or sum(len(g) for g in groups) < 3:
            continue
        shapiro_text = []
        for name, values in zip(names, groups):
            if 3 <= len(values) <= 5000:
                shapiro_text.append(f"{name}:Shapiro p={stats.shapiro(values).pvalue:.4g}")
        try:
            levene = stats.levene(*groups, center="median")
            levene_text = f"Levene p={levene.pvalue:.4g}"
        except Exception:
            levene = None
            levene_text = "Levene not estimable"
        normal_ok = all((stats.shapiro(g).pvalue > 0.05 if 3 <= len(g) <= 5000 else False) for g in groups)
        if normal_ok and levene is not None and levene.pvalue > 0.05:
            result = stats.f_oneway(*groups)
            grand = subset["GC_Percent_Valid_ATGC"].mean()
            ss_between = sum(len(g) * (np.mean(g) - grand) ** 2 for g in groups)
            ss_total = sum((v - grand) ** 2 for g in groups for v in g)
            eta = ss_between / ss_total if ss_total else np.nan
            add_row(analysis_name, "One-way ANOVA", "GC_Percent_Valid_ATGC", "Scientific_Name", len(subset), float(result.statistic), float(result.pvalue), eta, "; ".join(shapiro_text) + "; " + levene_text, "Parametric comparison used because normality and variance checks were acceptable.", "Records within a genome are not independent biological replicates; group sizes and record lengths differ.")
        else:
            result = stats.kruskal(*groups)
            n_total = len(subset)
            k = len(groups)
            epsilon = (float(result.statistic) - k + 1) / (n_total - k) if n_total > k else np.nan
            add_row(analysis_name, "Kruskal-Wallis", "GC_Percent_Valid_ATGC", "Scientific_Name", n_total, float(result.statistic), float(result.pvalue), epsilon, "; ".join(shapiro_text) + "; " + levene_text, "Non-parametric comparison selected because distribution/variance assumptions were not jointly acceptable.", "Records within a genome are not independent biological replicates; unequal group sizes and measurement precision remain.")
            if result.pvalue < 0.05 and multipletests is not None:
                pair_values: List[float] = []
                pair_meta: List[Tuple[str, str]] = []
                for i in range(len(names)):
                    for j in range(i + 1, len(names)):
                        u = stats.mannwhitneyu(groups[i], groups[j], alternative="two-sided")
                        pair_values.append(float(u.pvalue))
                        pair_meta.append((names[i], names[j]))
                adjusted = multipletests(pair_values, method="holm")[1] if pair_values else []
                for (a, b), raw_p, adj_p in zip(pair_meta, pair_values, adjusted):
                    add_row(analysis_name, "Pairwise Mann-Whitney U (Holm)", "GC_Percent_Valid_ATGC", f"{a} vs {b}", len(groups[names.index(a)]) + len(groups[names.index(b)]), np.nan, raw_p, np.nan, "Post-hoc only after significant Kruskal-Wallis.", "Adjusted p-value controls family-wise error across pairwise comparisons.", "Pseudo-replication and unequal record lengths.", adjusted=adj_p)
        if analysis_name == "All FASTA records":
            try:
                result_levene = stats.levene(*groups, center="median")
                add_row(analysis_name, "Levene variance test", "GC_Percent_Valid_ATGC", "Scientific_Name", len(subset), float(result_levene.statistic), float(result_levene.pvalue), np.nan, "Independent-group diagnostic only.", "Variance homogeneity diagnostic.", "Chromosome/record observations are nested within genomes.")
            except Exception:
                pass
    logger.info("Statistical testing completed on %d organism-level genomes and %d record-level observations", len(genome), len(records))
    return pd.DataFrame(rows)


def save_figure(fig: plt.Figure, figure_dir: Path, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(figure_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(figure_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def generate_figures(records: pd.DataFrame, genome: pd.DataFrame, stats_df: pd.DataFrame, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    g = genome.copy()
    g["Organism"] = pd.Categorical(g["Organism"], categories=ordered_names(g["Organism"].tolist()), ordered=True)
    g = g.sort_values("Organism")
    labels = g["Organism"].tolist()
    colors = [COLORS.get(x, "#7F7F7F") for x in labels]
    captions = []

    fig, ax = plt.subplots(figsize=(9, 5)); ax.bar(labels, g["Genome_GC_Percent"], color=colors); ax.set_ylabel("Whole-genome GC (%)"); ax.set_title("Whole-genome GC content by organism"); ax.tick_params(axis="x", rotation=35); save_figure(fig, figure_dir, "student1_gc_01_genome_gc_barplot"); captions.append("Figure 1. Whole-genome GC percentage calculated from summed G+C and valid A+T+G+C counts.")
    fig, ax = plt.subplots(figsize=(9, 5)); ax.bar(labels, g["Total_Sequence_Length"] / 1e6, color=colors); ax.set_ylabel("Genome size (Mb; total sequence length)"); ax.set_title("Genome size by organism"); ax.tick_params(axis="x", rotation=35); save_figure(fig, figure_dir, "student1_gc_02_genome_size_barplot"); captions.append("Figure 2. Total uploaded FASTA sequence length in megabases.")
    fig, ax = plt.subplots(figsize=(7, 5)); ax.scatter(g["Total_Sequence_Length"] / 1e6, g["Genome_GC_Percent"], c=colors, s=80); [ax.annotate(label, (x, y), xytext=(4, 4), textcoords="offset points", fontsize=8) for label, x, y in zip(labels, g["Total_Sequence_Length"] / 1e6, g["Genome_GC_Percent"])]; ax.set_xlabel("Genome size (Mb)"); ax.set_ylabel("Whole-genome GC (%)"); ax.set_title("Genome size versus whole-genome GC"); save_figure(fig, figure_dir, "student1_gc_03_size_gc_scatter"); captions.append("Figure 3. Each point is one supplied genome; annotations identify organisms.")
    fig, ax = plt.subplots(figsize=(7, 5)); ax.scatter(g["Total_Sequence_Length"] / 1e6, g["Genome_GC_Percent"], c=colors, s=80);
    if len(g) >= 2:
        slope, intercept = np.polyfit(g["Total_Sequence_Length"] / 1e6, g["Genome_GC_Percent"], 1); xx = np.linspace((g["Total_Sequence_Length"] / 1e6).min(), (g["Total_Sequence_Length"] / 1e6).max(), 100); ax.plot(xx, slope * xx + intercept, color="#222222", linestyle="--", label="OLS trend"); ax.legend()
    ax.set_xlabel("Genome size (Mb)"); ax.set_ylabel("Whole-genome GC (%)"); ax.set_title("Genome size versus GC with fitted trend"); save_figure(fig, figure_dir, "student1_gc_04_size_gc_trend"); captions.append("Figure 4. Ordinary least-squares trend is descriptive and is not a causal model.")
    fig, ax = plt.subplots(figsize=(10, 5)); data = [records.loc[records["Scientific_Name"] == label, "GC_Percent_Valid_ATGC"].dropna() for label in labels]; ax.boxplot(data, tick_labels=labels, patch_artist=True, boxprops={"facecolor": "#D9EAF7"}); ax.set_ylabel("Record GC (%)"); ax.set_title("Record-level GC content by organism"); ax.tick_params(axis="x", rotation=35); save_figure(fig, figure_dir, "student1_gc_05_record_gc_boxplot"); captions.append("Figure 5. Record-level GC distributions; records nested within an organism are not independent biological replicates.")
    fig, ax = plt.subplots(figsize=(10, 5)); violin = ax.violinplot(data, showmeans=False, showmedians=True); [body.set_facecolor("#9DC3E6") for body in violin["bodies"]]; ax.set_xticks(range(1, len(labels) + 1), labels, rotation=35); ax.set_ylabel("Record GC (%)"); ax.set_title("Record-level GC distributions with medians"); save_figure(fig, figure_dir, "student1_gc_06_record_gc_violin_boxplot"); captions.append("Figure 6. Violin distributions with median markers; sparse groups should be interpreted cautiously.")
    major = records[records["Sequence_Type"] == "chromosome"].copy(); major["Short_Label"] = major["Scientific_Name"] + " | " + major["Sequence_ID"]; major = major.sort_values("Sequence_Length", ascending=False).head(40); fig, ax = plt.subplots(figsize=(11, 6)); ax.barh(major["Short_Label"], major["Sequence_Length"] / 1e6, color=[COLORS.get(x, "#7F7F7F") for x in major["Scientific_Name"]]); ax.set_xlabel("Record length (Mb)"); ax.set_title("Main chromosome or linkage-group lengths"); ax.invert_yaxis(); save_figure(fig, figure_dir, "student1_gc_07_major_record_lengths"); captions.append("Figure 7. Lengths of records classified from headers as chromosomes or linkage groups; the 40 longest are shown for readability.")
    corr_cols = ["Total_Sequence_Length", "Genome_GC_Percent", "Valid_ATGC_Length", "Ambiguous_Percent", "N_Percent"]
    corr = g[corr_cols].corr(); corr.to_csv(figure_dir.parent / "tables" / "student1_gc_organism_correlation_matrix.csv", index=True)
    fig, ax = plt.subplots(figsize=(8, 5)); feature_names = ["Genome size", "Valid length", "Ambiguous %", "N %"]; vals = [corr.loc["Genome_GC_Percent", col] for col in ["Total_Sequence_Length", "Valid_ATGC_Length", "Ambiguous_Percent", "N_Percent"]]; ax.bar(feature_names, vals, color="#5B9BD5"); ax.axhline(0, color="black", linewidth=0.8); ax.set_ylim(-1, 1); ax.set_ylabel("Pearson r with genome GC (%)"); ax.set_title("Organism-level feature correlations with GC"); ax.tick_params(axis="x", rotation=25); save_figure(fig, figure_dir, "student1_gc_08_organism_feature_correlations"); captions.append("Figure 8. Pearson correlations among organism-level GC and selected genome characteristics.")
    fig, ax = plt.subplots(figsize=(9, 5)); ax.bar(labels, g["N_Percent"], color=colors, label="N %"); ax.bar(labels, g["Ambiguous_Percent"] - g["N_Percent"], bottom=g["N_Percent"], color="#A5A5A5", label="Other ambiguous %"); ax.set_ylabel("Percentage of total sequence"); ax.set_title("Ambiguous-base composition by genome"); ax.tick_params(axis="x", rotation=35); ax.legend(); save_figure(fig, figure_dir, "student1_gc_10_ambiguous_percent"); captions.append("Figure 10. N and other IUPAC ambiguous symbols as percentages of total sequence length.")
    fig, axes = plt.subplots(1, max(1, len(labels)), figsize=(4 * max(1, len(labels)), 4), squeeze=False); axes = axes.ravel();
    for ax, label in zip(axes, labels):
        vals = records.loc[(records["Scientific_Name"] == label) & records["GC_Percent_Valid_ATGC"].notna(), "GC_Percent_Valid_ATGC"]
        if len(vals) >= 3: stats.probplot(vals, dist="norm", plot=ax)
        else: ax.text(0.5, 0.5, "n < 3", ha="center", va="center")
        ax.set_title(label); ax.set_xlabel("Theoretical quantiles"); ax.set_ylabel("Ordered GC (%)")
    save_figure(fig, figure_dir, "student1_gc_11_qqplots"); captions.append("Figure 11. Q-Q diagnostics for record-level GC distributions; groups with fewer than three records are marked n < 3.")
    fig, ax = plt.subplots(figsize=(8, 5)); ax.scatter(records["Sequence_Length"] / 1e6, records["GC_Percent_Valid_ATGC"], c=[COLORS.get(x, "#7F7F7F") for x in records["Scientific_Name"]], alpha=0.65, s=22); ax.set_xlabel("Record length (Mb)"); ax.set_ylabel("Record GC (%)"); ax.set_title("GC percentage versus record length"); save_figure(fig, figure_dir, "student1_gc_12_gc_vs_sequence_length"); captions.append("Figure 12. Record-level GC percentage versus sequence length.")
    main_g = records[records["Sequence_Type"] == "chromosome"]; fig, ax = plt.subplots(figsize=(10, 5)); main_data = [main_g.loc[main_g["Scientific_Name"] == label, "GC_Percent_Valid_ATGC"].dropna() for label in labels if (main_g["Scientific_Name"] == label).any()]; main_labels = [label for label in labels if (main_g["Scientific_Name"] == label).any()]; ax.boxplot(main_data, tick_labels=main_labels, patch_artist=True, boxprops={"facecolor": "#C6E0B4"}); ax.set_ylabel("Main-chromosome GC (%)"); ax.set_title("GC comparison for confidently identified main chromosomes"); ax.tick_params(axis="x", rotation=35); save_figure(fig, figure_dir, "student1_gc_13_main_chromosome_gc"); captions.append("Figure 13. GC values for records explicitly labelled chromosome or linkage group.")
    organ = records[records["Sequence_Type"].isin(["mitochondrion", "chloroplast"])]; nuclear = records[~records["Sequence_Type"].isin(["mitochondrion", "chloroplast"])]; organ_rows = []
    for label in labels:
        for category, subset in [("Nuclear records", nuclear[nuclear["Scientific_Name"] == label]), ("Organellar records", organ[organ["Scientific_Name"] == label])]:
            if not subset.empty and subset["Valid_ATGC_Length"].sum() > 0: organ_rows.append((label, category, 100 * subset["GC_Count"].sum() / subset["Valid_ATGC_Length"].sum()))
    fig, ax = plt.subplots(figsize=(10, 5)); odf = pd.DataFrame(organ_rows, columns=["Organism", "Category", "GC"]);
    if not odf.empty:
        pivot = odf.pivot(index="Organism", columns="Category", values="GC").reindex(labels); pivot.plot(kind="bar", ax=ax, color=["#4472C4", "#ED7D31"])
    ax.set_ylabel("Weighted GC (%)"); ax.set_title("Nuclear versus organellar GC where available"); ax.tick_params(axis="x", rotation=35); save_figure(fig, figure_dir, "student1_gc_14_nuclear_organellar_gc"); captions.append("Figure 14. Weighted nuclear and organellar GC; an organellar category is shown only where supplied headers identify it.")
    (figure_dir / "student1_gc_figure_captions.txt").write_text("\n".join(captions) + "\n", encoding="utf-8")


def write_quality_report(records: pd.DataFrame, genome: pd.DataFrame, audit: pd.DataFrame, stats_df: pd.DataFrame, results_dir: Path) -> None:
    checks: List[str] = []
    def check(condition: bool, passed: str, failed: str) -> None:
        checks.append(f"[PASS] {passed}" if condition else f"[FAIL] {failed}")
    check((records["Sequence_Length"] == records[["A_Count", "T_Count", "G_Count", "C_Count", "N_Count", "Other_Ambiguous_Count", "Unexpected_Symbol_Count"]].sum(axis=1)).all(), "All symbol counts sum to total sequence length.", "A symbol count does not equal total sequence length.")
    check((records["Valid_ATGC_Length"] == records[["A_Count", "T_Count", "G_Count", "C_Count"]].sum(axis=1)).all(), "Valid ATGC length equals A+T+G+C.", "Valid ATGC length mismatch.")
    check((records["GC_Count"] == records[["G_Count", "C_Count"]].sum(axis=1)).all(), "GC count equals G+C.", "GC count mismatch.")
    check(records["GC_Percent_Valid_ATGC"].dropna().between(0, 100).all(), "All primary GC percentages are between 0 and 100.", "A primary GC percentage is outside 0-100.")
    check((genome["Total_Sequence_Length"] == genome["Valid_ATGC_Length"] + genome["Ambiguous_Base_Count"] + genome["Unexpected_Symbol_Count"]).all(), "Genome totals equal valid, ambiguous, and unexpected symbol totals.", "A genome total does not equal its component totals.")
    check((audit["Empty_File"] == "PASS").all(), "No empty FASTA files.", "An empty FASTA file was detected.")
    check((audit["Zero_Length_Records"] == 0).all(), "No zero-length FASTA records.", "A zero-length FASTA record was detected.")
    check((audit["Duplicate_Sequence_IDs"] == 0).all(), "No duplicate sequence identifiers.", "Duplicate sequence identifiers were detected.")
    check((audit["Unexpected_Symbol_Count"] == 0).all(), "No unexpected non-IUPAC sequence symbols.", "Unexpected sequence symbols were detected.")
    warnings = ["[WARNING] Four supplied files are not among the five species named in the project brief; their header-derived identities are retained and reported.", "[WARNING] Source database is not explicitly available in the provided files.", "[WARNING] Chromosome/record observations are nested within genomes and are not independent biological replicates.", "[WARNING] Organism-level correlation uses n=5 available genomes and has low statistical power."]
    text = ["Student 1 GC analysis quality-control report", "", "Checks:"] + checks + ["", "Warnings:"] + warnings + ["", f"Records analysed: {len(records)}", f"Genomes analysed: {len(genome)}", f"Statistical result rows: {len(stats_df)}", "", "Corrective actions: FASTA letters were counted case-insensitively; ambiguous symbols were excluded from the primary GC denominator; no record was silently removed."]
    (results_dir / "logs" / "student1_gc_quality_control_report.txt").write_text("\n".join(text) + "\n", encoding="utf-8")


def write_csv(df: pd.DataFrame, path: Path, columns: Optional[List[str]] = None) -> None:
    if columns:
        for col in columns:
            if col not in df.columns:
                df[col] = np.nan
        df = df[columns]
    df.to_csv(path, index=False, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=INPUT_DIR,
        help="Input FASTA directory",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Output directory",
    )
    args = parser.parse_args()
    out = args.output_dir
    table_dir, figure_dir, log_dir = out / "tables", out / "figures", out / "logs"
    for directory in [table_dir, figure_dir, log_dir]: directory.mkdir(parents=True, exist_ok=True)
    logger = configure_logging(log_dir / "student1_gc_analysis.log")
    files = discover_fasta_files(args.input_dir)
    if not files: raise SystemExit("No FASTA files found")
    all_records: List[Dict[str, Any]] = []; audits: List[Dict[str, Any]] = []
    for path in files:
        logger.info("Validating %s", path.name)
        recs, audit = validate_fasta_file(path, logger); all_records.extend(recs); audits.append(audit)
    records = pd.DataFrame(all_records); audit_df = pd.DataFrame(audits)
    genome = infer_genome_rows(records, audit_df)
    gc_results = calculate_gc_results(records, genome)
    stat_results = run_statistical_tests(records, genome, logger)
    sequence_columns = ["Scientific_Name", "Input_File", "Sequence_ID", "Sequence_Description", "Sequence_Type", "Sequence_Length", "A_Count", "T_Count", "G_Count", "C_Count", "N_Count", "Other_Ambiguous_Count", "Unexpected_Symbol_Count", "Valid_ATGC_Length", "GC_Count", "GC_Percent_Valid_ATGC", "GC_Percent_All_Bases", "File_Size_Bytes", "File_Size_MB", "Duplicate_Sequence_ID"]
    genome_columns = ["Organism", "Scientific_Name", "Input_File", "Expected_Species_Match", "Source", "Assembly_Version", "Accession", "Number_of_FASTA_Records", "Number_of_Chromosomes", "Number_of_Scaffolds", "Number_of_Contigs", "Number_of_Organellar_Records", "Number_of_Unplaced_Records", "Number_of_Unlocalized_Records", "Total_Sequence_Length", "Valid_ATGC_Length", "Ambiguous_Base_Count", "Unexpected_Symbol_Count", "Genome_GC_Count", "Genome_GC_Percent", "Genome_GC_Percent_All_Bases", "File_Size_Bytes", "File_Size_MB", "A_Count", "T_Count", "G_Count", "C_Count", "N_Count", "Other_Ambiguous_Count", "Ambiguous_Percent", "N_Percent"]
    gc_columns = ["Organism", "Analysis", "Record_Count", "Weighted_GC_Percent", "Unweighted_Mean_Record_GC_Percent", "Median_Record_GC_Percent", "Standard_Deviation", "Minimum", "Maximum", "Range", "Interquartile_Range", "Coefficient_of_Variation_Percent", "Ambiguous_Percent", "N_Percent"]
    stat_columns = ["Analysis", "Test", "Variables", "Grouping_Variable", "Sample_Size", "Statistic", "P_Value", "Adjusted_P_Value", "Effect_Size", "Confidence_Interval", "Assumptions", "Interpretation", "Limitations"]
    write_csv(records, table_dir / "student1_gc_sequence_summary.csv", sequence_columns)
    write_csv(genome, table_dir / "student1_gc_genome_summary.csv", genome_columns)
    write_csv(gc_results, table_dir / "student1_gc_results.csv", gc_columns)
    write_csv(stat_results, table_dir / "student1_gc_statistical_results.csv", stat_columns)
    write_csv(audit_df, table_dir / "student1_gc_input_file_audit.csv")
    generate_figures(records, genome, stat_results, figure_dir)
    write_quality_report(records, genome, audit_df, stat_results, out)
    logger.info("Completed analysis: %d records, %d organisms, %d statistical rows", len(records), len(genome), len(stat_results))


if __name__ == "__main__":
    main()
