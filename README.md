# Comparative Genomics: GC Content & K-mer Analysis

Coursework repository for comparative genomics across five eukaryotic assemblies. The project is split into two independent modules that share the same NCBI genome packages under `data/`.

| Module | Focus |
|--------|--------|
| [`member1_gc_analysis/`](member1_gc_analysis/) | Genome characteristics, GC content, statistics, and visualization |
| [`member2_kmer_analysis/`](member2_kmer_analysis/) | K-mer profiles, entropy, similarity, clustering, TATA-box scan, and plots |

## Organisms

Assemblies are read from NCBI Datasets packages. Header inspection identifies:

| Dataset folder | Organism |
|----------------|----------|
| `data/ncbi_dataset1/` | *Caenorhabditis elegans* |
| `data/ncbi_dataset2/` | *Neurospora crassa* |
| `data/ncbi_dataset3/` | *Yarrowia lipolytica* |
| `data/ncbi_dataset4/` | *Dictyostelium discoideum* |
| `data/ncbi_dataset5/` | *Eremothecium coryli* |

Each package follows the NCBI layout (`ncbi_dataset/data/.../*_genomic.fna`) with accompanying `data_summary.tsv` and assembly reports.

## Repository layout

```text
comparative-genomics-kmer-analysis/
├── data/                      # Shared NCBI genome packages
│   ├── ncbi_dataset1/
│   ├── ncbi_dataset2/
│   ├── ncbi_dataset3/
│   ├── ncbi_dataset4/
│   └── ncbi_dataset5/
├── member1_gc_analysis/       # Student 1 — GC & genome-size analysis
├── member2_kmer_analysis/     # Student 2 — k-mer analysis pipeline
├── .gitignore
└── README.md
```

## Requirements

- Python **3.10+**
- Separate virtual environments (recommended) per module, or one env with both requirement sets

```bash
python -m venv .venv

# Windows
.venv\Scripts\activate

# macOS / Linux
source .venv/bin/activate

python -m pip install --upgrade pip
```

Install dependencies for the module you need:

```bash
pip install -r member1_gc_analysis/requirements.txt
# and/or
pip install -r member2_kmer_analysis/requirements.txt
```

**Member 1:** NumPy, pandas, SciPy, Matplotlib, statsmodels, Biopython  
**Member 2:** Biopython, NumPy, pandas, SciPy, scikit-learn, Matplotlib, seaborn

## Quick start

### Member 1 — GC analysis

Place flat `.fna` / `.fasta` / `.fa` files in `member1_gc_analysis/data/raw/` (copy or symlink assemblies from `data/ncbi_dataset*/` if needed), then:

```bash
cd member1_gc_analysis
python scripts/student1_gc_analysis.py --input-dir data/raw --output-dir results
```

See [`member1_gc_analysis/README.md`](member1_gc_analysis/README.md) for GC definitions, outputs, and troubleshooting.

### Member 2 — K-mer analysis

Scripts expect to be run from `member2_kmer_analysis/` and discover genomes under the shared `../data/ncbi_dataset*` folders (or use `--demo`):

```bash
cd member2_kmer_analysis

python scripts/kmer_analysis.py          # or: python scripts/kmer_analysis.py --demo
python scripts/entropy.py
python scripts/similarity.py
python scripts/clustering.py
python scripts/tata_box.py
python scripts/visualization.py
```

See [`member2_kmer_analysis/README.md`](member2_kmer_analysis/README.md) for flags (`-k`, `--input-dir`, `--output-dir`, etc.) and method notes.

## Outputs

| Location | Contents |
|----------|----------|
| `member1_gc_analysis/results/tables/` | Genome/GC summaries, statistics, audit CSVs |
| `member1_gc_analysis/results/figures/` | GC and size plots (PNG + PDF) |
| `member1_gc_analysis/report/` | Persian report and oral presentation guide |
| `member2_kmer_analysis/results/` | K-mer, entropy, distance/similarity, PCA, TATA CSVs |
| `member2_kmer_analysis/results/figures/` | Heatmaps, PCA, dendrogram, entropy chart |
| `member2_kmer_analysis/report/` | Persian report and oral presentation guide |

## Data notes

- Genomic FASTA files under `data/ncbi_dataset*/` are part of this repository’s working dataset.
- Large archives (`.zip`, `.tar.gz`) and common sequencing/annotation dump formats are listed in `.gitignore` so accidental bulk downloads are not committed.
- Regenerated analysis tables and figures can be overwritten by re-running the module scripts.

## License / course use

Prepared for a comparative genomics coursework assignment. Cite NCBI Datasets for the source assemblies when presenting or submitting results.
