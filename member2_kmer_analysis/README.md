# Member 2 — K-mer Analysis Module

Comparative genomics pipeline for k-mer profiling, Shannon entropy, similarity metrics, clustering, TATA box detection, and visualization.

## Overview

This module analyzes DNA sequences from FASTA files and produces:

| Output | Description |
|--------|-------------|
| `kmer_summary.csv` | Per-sequence k-mer counts and diversity |
| `entropy.csv` | Shannon entropy of k-mer distributions |
| `distance_matrix.csv` | Pairwise Euclidean distance matrix |
| `cosine_similarity.csv` | Pairwise cosine similarity matrix |
| `tata_box.csv` | TATA box motif detection in promoter regions |

| Figure | Description |
|--------|-------------|
| `figures/heatmap.png` | Euclidean distance heatmap |
| `figures/similarity_matrix.png` | Cosine similarity heatmap |
| `figures/pca.png` | PCA scatter plot |
| `figures/dendrogram.png` | Hierarchical clustering dendrogram |
| `figures/entropy_comparison.png` | Shannon entropy bar chart |

## Setup

```bash
cd member2_kmer_analysis
pip install -r requirements.txt
```

## Input Data

Place NCBI datasets under the project-level `data/` directory using the standard NCBI folder layout:

```
comparative-genomics-kmer-analysis/
├── data/
│   ├── ncbi_dataset1/    # e.g. Caenorhabditis elegans
│   ├── ncbi_dataset2/    # e.g. Neurospora crassa
│   ├── ncbi_dataset3/    # e.g. Yarrowia lipolytica
│   ├── ncbi_dataset4/    # e.g. Dictyostelium discoideum
│   └── ncbi_dataset5/    # e.g. Eremothecium coryli
└── member2_kmer_analysis/
```

The pipeline auto-discovers `ncbi_dataset*` folders, reads organism names from `data_summary.tsv`, and uses the RefSeq (`GCF_*`) or GenBank (`GCA_*`) `*_genomic.fna` assembly in each folder. All contigs/chromosomes in an assembly are pooled into one k-mer profile per species.

Alternatively, pass a custom path with `--input-dir` or place flat `.fa`/`.fasta`/`.fna` files directly in `data/`.

### Demo Mode

Run with `--demo` to generate sample sequences in `results/demo_input/demo_sequences.fasta`:

```bash
python scripts/kmer_analysis.py --demo
```

## Usage

Run scripts from the `member2_kmer_analysis/` directory:

```bash
# 1. K-mer counting and summary
python scripts/kmer_analysis.py --demo

# 2. Shannon entropy
python scripts/entropy.py

# 3. Similarity matrices (Euclidean + cosine)
python scripts/similarity.py

# 4. PCA and hierarchical clustering
python scripts/clustering.py

# 5. TATA box motif analysis
python scripts/tata_box.py

# 6. Generate all figures (run entropy.py first for entropy figure)
python scripts/visualization.py
```

### Full pipeline (demo)

```bash
python scripts/kmer_analysis.py --demo
python scripts/entropy.py
python scripts/similarity.py
python scripts/clustering.py
python scripts/tata_box.py
python scripts/visualization.py
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--input-dir` | `../../data` | Directory with FASTA files |
| `--output-dir` | `results/` | Output directory for CSVs and figures |
| `-k` | `4` | K-mer size |
| `--demo` | off | Create demo FASTA data |
| `--promoter-length` | `500` | TATA scan window (tata_box.py only) |
| `--linkage` | `ward` | Clustering method (clustering.py only) |

## Scripts

| Script | Purpose |
|--------|---------|
| `kmer_analysis.py` | K-mer counting, frequency vectors, summary statistics |
| `entropy.py` | Shannon entropy (bits) of k-mer distributions |
| `similarity.py` | Euclidean distance and cosine similarity matrices |
| `clustering.py` | PCA dimensionality reduction and hierarchical clustering |
| `tata_box.py` | Canonical TATA box (TATAAA) detection in 5' regions |
| `visualization.py` | Heatmaps, PCA plot, dendrogram, entropy comparison |

## Methods

### K-mer Analysis
Overlapping k-mers (default k=4) are counted per sequence. Frequencies are normalized to produce comparable profiles across sequences of different lengths.

### Shannon Entropy
H = −Σ pᵢ log₂(pᵢ) computed over the k-mer frequency distribution. Higher entropy indicates more uniform k-mer usage.

### Similarity Metrics
- **Euclidean distance**: L2 norm between k-mer frequency vectors
- **Cosine similarity**: Angular similarity (0–1) between frequency vectors

### Clustering
- **PCA**: Standardized k-mer matrix projected to 2 principal components
- **Hierarchical clustering**: Ward linkage on Euclidean distance matrix

### TATA Box
Scans the first 500 bp (configurable) for the canonical eukaryotic promoter motif `TATAAA` and variants.

## Dependencies

- Biopython — FASTA parsing
- NumPy / Pandas — numerical and tabular data
- SciPy — distance metrics, hierarchical clustering
- Scikit-learn — PCA, cosine similarity, scaling
- Matplotlib / Seaborn — visualization

## Output Structure

```
member2_kmer_analysis/
├── results/
│   ├── kmer_summary.csv
│   ├── entropy.csv
│   ├── distance_matrix.csv
│   ├── cosine_similarity.csv
│   ├── tata_box.csv
│   ├── pca_coordinates.csv
│   ├── pca_variance.csv
│   ├── linkage_matrix.csv
│   └── figures/
│       ├── heatmap.png
│       ├── similarity_matrix.png
│       ├── pca.png
│       ├── dendrogram.png
│       └── entropy_comparison.png
├── scripts/
└── requirements.txt
```
