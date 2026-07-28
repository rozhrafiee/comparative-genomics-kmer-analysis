# Comparative Genomics: GC Content & K-mer Analysis

Coursework repository for comparative genomics across five eukaryotic assemblies.

| Module | Focus | Owner |
|--------|--------|--------|
| `member1_gc_analysis/` | GC content, genome metadata, GC figures/stats | Sadra Kasai |
| `member2_kmer_analysis/` | K-mers, entropy, similarity, clustering, TATA, stats | Rozhina Rafiee |
| `presentation/` | Seminar deck + Persian speaking guide | Both |

## Organisms (selected NCBI assemblies)

| Dataset | Organism |
|---------|----------|
| `data/ncbi_dataset1/` | *Caenorhabditis elegans* |
| `data/ncbi_dataset2/` | *Neurospora crassa* |
| `data/ncbi_dataset3/` | *Yarrowia lipolytica* |
| `data/ncbi_dataset4/` | *Dictyostelium discoideum* |
| `data/ncbi_dataset5/` | *Eremothecium coryli* |

## Setup

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
pip install -r member1_gc_analysis/requirements.txt
pip install -r member2_kmer_analysis/requirements.txt
```

## Run

**Member 1 (GC):**
```bash
cd member1_gc_analysis
python scripts/student1_gc_analysis.py --input-dir data/raw --output-dir results
```

**Member 2 (k-mer):**
```bash
cd member2_kmer_analysis
python scripts/student2_kmer_analysis.py
```

**Presentation:** open `presentation/presentation.html`  
**Speaking guide:** `presentation/راهنمای_ارائه.md`

## Outputs

| Path | Contents |
|------|----------|
| `member1_gc_analysis/results/` | GC tables, figures, logs |
| `member1_gc_analysis/report/student1_gc_report_fa.md` | Member 1 report |
| `member2_kmer_analysis/results/tables/` | All CSV outputs (k-mer, similarity, PCA, TATA, stats) |
| `member2_kmer_analysis/results/figures/` | All figures |
| `member2_kmer_analysis/report/student2_kmer_report.md` | Member 2 report |
