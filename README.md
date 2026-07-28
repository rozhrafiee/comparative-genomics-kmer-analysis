# Comparative Genomics: GC Content & K-mer Analysis

Coursework on five eukaryotic NCBI assemblies. Each student has **one analysis script** (same pattern).

| Module | Focus | Script | Owner |
|--------|--------|--------|--------|
| [`member1_gc_analysis/`](member1_gc_analysis/) | Data summary, GC content, GC stats & figures | `scripts/student1_gc_analysis.py` | Sadra Kasai |
| [`member2_kmer_analysis/`](member2_kmer_analysis/) | K-mers, similarity, PCA, clustering, stats, TATA | `scripts/student2_kmer_analysis.py` | Rozhina Rafiee |

## Organisms

| Dataset | Organism |
|---------|----------|
| `data/ncbi_dataset1/` | *Caenorhabditis elegans* |
| `data/ncbi_dataset2/` | *Neurospora crassa* |
| `data/ncbi_dataset3/` | *Yarrowia lipolytica* |
| `data/ncbi_dataset4/` | *Dictyostelium discoideum* |
| `data/ncbi_dataset5/` | *Eremothecium coryli* |

Shared genomes live under [`data/`](data/). Member 1 may also use a local `data/raw/` copy of the FASTAs.

## Setup

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
pip install -r member1_gc_analysis/requirements.txt
pip install -r member2_kmer_analysis/requirements.txt
```

## Run

**Member 1 — GC (brief §§1–2):**
```bash
cd member1_gc_analysis
python scripts/student1_gc_analysis.py --input-dir data/raw --output-dir results
```

**Member 2 — K-mer (brief §§3–6):**
```bash
cd member2_kmer_analysis
python scripts/student2_kmer_analysis.py
```

Member 2 reads genomes from repo `data/` and (read-only) Member 1’s GC summary for correlations.

## Layout

```
data/                          # NCBI assemblies (shared)
member1_gc_analysis/
  scripts/student1_gc_analysis.py
  results/tables|figures|logs/
  report/student1_gc_report_fa.md
member2_kmer_analysis/
  scripts/student2_kmer_analysis.py
  results/tables/              # all CSVs
  results/figures/             # all plots
  report/student2_kmer_report.md
```

## Outputs

| Path | Contents |
|------|----------|
| `member1_gc_analysis/results/tables/` | Genome / sequence / GC / stats CSVs |
| `member1_gc_analysis/results/figures/` | GC figures |
| `member1_gc_analysis/report/student1_gc_report_fa.md` | Member 1 report (Persian) |
| `member2_kmer_analysis/results/tables/` | K-mer, distance, PCA, TATA, stats CSVs |
| `member2_kmer_analysis/results/figures/` | Heatmaps, PCA, dendrogram, stats plots |
| `member2_kmer_analysis/report/student2_kmer_report.md` | Member 2 report (English) |
