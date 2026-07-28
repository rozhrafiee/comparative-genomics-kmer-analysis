# Comparative Genomics — Student 2 Report (K-mer Analysis)

**Student:** Rozhina Rafiee  
**Scope:** Brief §§3–6 (k-mer, similarity, statistics, TATA). GC / data collection = Student 1 (read-only).

**Run:**
```bash
cd member2_kmer_analysis
python scripts/student2_kmer_analysis.py
```

Outputs: `results/tables/`, `results/figures/`.

---

## 1. Data

Five NCBI assemblies under `data/ncbi_dataset1…5` (student-chosen; brief species list was examples only):

| Dataset | Assembly | Organism |
|---|---|---|
| ncbi_dataset1 | GCF_000002985.6_WBcel235 | *C. elegans* |
| ncbi_dataset2 | GCF_000182925.2_NC12 | *N. crassa* |
| ncbi_dataset3 | GCF_001761485.1_ASM176148v1 | *Y. lipolytica* |
| ncbi_dataset4 | GCF_000004695.1_dicty_2.7 | *D. discoideum* |
| ncbi_dataset5 | GCA_000710315.1 | *E. coryli* |

Genome metadata / GC% come from Student 1 tables (not recomputed here).

---

## 2. Methods (one script)

`scripts/student2_kmer_analysis.py`:

1. Load genomic FASTA (GCF preferred, else GCA); keep A/C/G/T only.
2. Count overlapping k-mers for **k = 3, 4, 5** (pooled per species).
3. Main comparison at **k = 4**: frequency vectors → Euclidean distance, cosine similarity, PCA, Ward dendrogram.
4. Shannon entropy and k-mer diversity (`unique / total`).
5. Species-exclusive k-mers = set difference across the five species.
6. TATA: scan first **500 bp** of pooled sequence for `TATAAA`.
7. Stats: Pearson/Spearman (+ permutation) at species level; ANOVA / Kruskal–Wallis / Mann–Whitney on per-chromosome metrics.

---

## 3. K-mer results (§3)

**Most frequent k-mer (rank 1):**

| Organism | k=3 | k=4 | k=5 |
|---|---|---|---|
| *C. elegans* | TTT | AAAA | AAAAA |
| *N. crassa* | TTT | TTTT | TTTTT |
| *Y. lipolytica* | CAA | TTTT | AAAAA |
| *D. discoideum* | AAA | AAAA | AAAAA |
| *E. coryli* | AAA | AAAA | AAAAA |

| Question | Answer |
|---|---|
| Unique 4-mers | All 256 possible in every species |
| Species-exclusive k-mers (k=3–5) | None |
| Highest diversity (k=4) | *E. coryli* (2.84×10⁻⁵) |
| Highest entropy (k=4) | *N. crassa* (7.965 bits) |
| Lowest entropy (k=4) | *D. discoideum* (6.917 bits) |

---

## 4. Similarity (§4)

From Euclidean distance on k=4 frequency vectors:

| Question | Answer |
|---|---|
| Closest pair | *N. crassa* – *Y. lipolytica* (0.0154) |
| Farthest pair | *Y. lipolytica* – *D. discoideum* (0.1023) |
| Animal cluster? | Not testable — only one animal (*C. elegans*) |
| Plant vs animals? | Not testable — no plant in this dataset |
| Match biology? | Closest pair are both fungi; *D. discoideum* (A/T-rich) is distant — consistent with composition |

PCA (k=4): **PC1 ≈ 70.4%**, **PC2 ≈ 12.4%** of variance.

---

## 5. Statistics (§5)

Species-level correlations (n=5 → low power; exploratory):

| Pair | Pearson r | p |
|---|---:|---:|
| genome size ↔ k-mer diversity | −0.743 | 0.150 |
| GC% ↔ genome size | −0.249 | 0.687 |
| GC% ↔ Shannon entropy | **0.966** | **0.007** |
| TATA count ↔ GC% | 0.112 | 0.858 |

Group tests on per-record (chromosome/scaffold) values — significant between species, but records within a genome are nested (not independent replicates):

| Metric | ANOVA F | p |
|---|---:|---:|
| Shannon entropy | 19.81 | 1.1×10⁻¹¹ |
| k-mer diversity | 13.38 | 1.4×10⁻⁸ |

GC difference between species is Student 1’s analysis; here GC is only used for correlations above.

---

## 6. TATA box (§6)

Raw-sequence scan of first 500 bp (not annotated promoters):

| Organism | Canonical TATAAA count |
|---|---:|
| *C. elegans* | 0 |
| *N. crassa* | 5 |
| *Y. lipolytica* | 0 |
| *D. discoideum* | 2 |
| *E. coryli* | 0 |

| Question | Answer |
|---|---|
| % of genes with TATA? | Cannot answer — no GFF/GTF gene annotation |
| Frequency differs by species? | Yes in this window (*N. crassa* 5, *D. discoideum* 2, others 0) |
| TATA+ genes higher expression? | Cannot answer — no expression data |
| Similar patterns in close species? | Closest fungi pair both lack TATAAA in this window; limited evidence |
| Linked to GC%? | Weak / non-significant (r≈0.11, p≈0.86) |

---

## 7. Outputs

**Script:** `scripts/student2_kmer_analysis.py`

**Tables (`results/tables/`):**  
`student2_kmer_top_by_k.csv`, `student2_kmer_diversity_entropy.csv`, `student2_unique_kmers.csv`, `student2_kmer_feature_matrix.csv`, `student2_distance_matrix.csv`, `student2_cosine_similarity.csv`, `student2_pca.csv`, `student2_tata_box.csv`, `student2_species_metrics.csv`, `student2_per_record_metrics.csv`, `student2_statistical_results.csv`

**Figures (`results/figures/`):**  
k-mer heatmap, distance heatmap, cosine matrix, PCA, dendrogram, entropy / diversity / top-k bars, correlation matrix & scatter, entropy & diversity boxplots, ANOVA summary, Mann–Whitney heatmap

---

## 8. Conclusion

K-mer profiles separate these five genomes mainly by relative oligomer frequencies (not exclusive k-mers at k≤5). The closest pair is fungal (*N. crassa* / *Y. lipolytica*); *D. discoideum* is the most distant. GC% tracks Shannon entropy strongly at the species level (small n). TATA hits are descriptive only for a short 5′ window. Plant/animal brief questions do not apply to this organism set.
