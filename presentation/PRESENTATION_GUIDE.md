# Oral Presentation Guide (for the team)

Use this guide while presenting. It is **separate** from the Gamma deck (`Gamma_Presentation_EN.docx`).

**How to use:** for each part, show the listed plots, say the short talking points, and keep the optional deeper quotes only if the teacher asks.

Suggested total time: **10–15 minutes** (≈1–2 min per major part).

---

## Before you start (30 seconds)

**Say:**
> We compared five NCBI genome assemblies using GC content and alignment-free k-mer profiles in Python. The brief’s species list were examples; our real set is what’s in the repository.

**Show:** title slide only (no plot yet).

**If asked why organisms differ from the brief:**
> Four FASTA files did not match the example list. We kept the real files, documented the mismatch, and analyzed those five genomes.

---

## Part 0 — Organisms (1 minute)

| Species | One-line intro to say |
|---------|------------------------|
| *C. elegans* | Model nematode; only species that matched the brief; largest genome (~100 Mb). |
| *N. crassa* | Model filamentous fungus; highest k-mer entropy; close to *Y. lipolytica*. |
| *Y. lipolytica* | Industrial yeast; highest GC (~49%). |
| *D. discoideum* | Social amoeba; AT-rich; lowest GC/entropy; separates in PCA. |
| *E. coryli* | Compact yeast/fungus; smallest assembly (~9 Mb); scaffolds only. |

**Plots:** none required here (or briefly flash genome-size bar if you want a visual hook).

---

## Part 1 — Data collection (1–2 minutes)

### Plots to show
1. `student1_gc_02_genome_size_barplot.png` — **main**
2. `student1_gc_07_major_record_lengths.png` — optional if asked about chromosomes

### What to say (quote-style talking points)
> All assemblies come from NCBI Datasets as genomic FASTA. We report accession, assembly label, record count, total size, and file size in one table.

> *C. elegans* is about 100 Mb; *E. coryli* is about 9 Mb. Record lengths differ a lot within genomes.

### Numbers to have ready
| Species | Accession | Assembly | Size | File |
|---------|-----------|----------|------|------|
| *C. elegans* | GCA_000002985.3 | WBcel235 | 100.27 Mb | 101.5 MB |
| *N. crassa* | GCA_000182925.2 | NC12 | 41.10 Mb | 41.6 MB |
| *Y. lipolytica* | GCA_001761485.1 | ASM176148v1 | 20.55 Mb | 20.8 MB |
| *D. discoideum* | GCA_000004695.1 | dicty_2.7 | 34.13 Mb | 34.6 MB |
| *E. coryli* | GCA_000710315.1 | E. coryli | 9.09 Mb | 9.2 MB |

### Code line (if teacher asks “how?”)
> Streaming FASTA parse in Python; metadata from headers and NCBI package files — not manual guessing.

---

## Part 2 — GC content (2–3 minutes)

### Plots to show (in this order)
1. `student1_gc_01_genome_gc_barplot.png` — **must**
2. `student1_gc_03_size_gc_scatter.png` — **must**
3. `student1_gc_05_record_gc_boxplot.png` — **must** (answers “same GC per chromosome?”)
4. Optional: `student1_gc_06_record_gc_violin_boxplot.png`, `student1_gc_04_size_gc_trend.png`, `student1_gc_13_main_chromosome_gc.png`

### What to say
> Primary metric: GC% = 100 × (G+C) / (A+T+G+C). Ambiguous bases are excluded from the denominator. Whole-genome GC is length-weighted from summed counts.

> Highest GC: *Y. lipolytica* ≈ 49%. Lowest: *D. discoideum* ≈ 22%.

> Size versus GC correlation is weak and non-significant with five genomes (Pearson ≈ −0.25, p ≈ 0.69). So we do **not** claim a size–GC law from this sample.

> Boxplots show chromosomes/records inside one species do **not** share one GC value.

### Skip unless asked
- Ambiguous % plot (`student1_gc_10_…`) — *E. coryli* has more Ns
- Nuclear/organellar (`student1_gc_14_…`) — only where headers label organelles

---

## Part 3 — K-mers & entropy (2 minutes)

### Plots to show
1. `entropy_comparison.png` — **must**

### What to say
> A k-mer is a substring of length k. We used overlapping 4-mers, counted only A/C/G/T, and normalized frequencies so genome length does not dominate.

> All five species use all 256 possible 4-mers. So there is no 4-mer unique to only one species at k=4.

> Most common 4-mer is AAAA in AT-rich genomes and TTTT in the two fungi.

> Shannon entropy is highest in *N. crassa* (~7.97 bits) and lowest in *D. discoideum* (~6.92). Lower entropy means a less even k-mer distribution.

### If asked about 3-mer / 5-mer
> The same script supports `-k 3` and `-k 5`. Our reported figures and matrices are for the primary run at k=4.

### If asked “does similar k-mer mean related evolution?”
> It means similar oligonucleotide composition. It is evidence of compositional closeness, not a formal phylogenetic proof.

---

## Part 4 — Similarity, PCA, clustering (2–3 minutes)

### Plots to show (in this order)
1. `heatmap.png` — Euclidean distance — **must**
2. `similarity_matrix.png` — cosine — **must**
3. `pca.png` — **must**
4. `dendrogram.png` — **must**

### What to say
> Each genome becomes a 256-dimensional frequency vector. We compare vectors with Euclidean distance and cosine similarity.

> Closest pair: *N. crassa* and *Y. lipolytica* (distance ≈ 0.015, cosine ≈ 0.97).

> *D. discoideum* is farthest from those fungi and sits apart on PCA. PC1 alone explains about 70% of variance.

> The Ward dendrogram recovers the same fungal cluster. This is hierarchical clustering of composition — not a calibrated species tree.

### Brief questions that need careful answers
| Teacher question | Short answer |
|------------------|--------------|
| Do animals cluster? | Only *C. elegans* is an animal model here; brief animals were not in our files. |
| Does the plant separate? | *A. thaliana* was an example only; not in our dataset. |
| Match biology? | Broadly yes for fungi-close / amoeba-separate; still not formal phylogeny. |

---

## Part 5 — Statistics (1–2 minutes)

### Plots to show
1. `student1_gc_08_organism_feature_correlations.png` — **must**
2. `student1_gc_11_qqplots.png` — optional (shows we checked assumptions)

### What to say
> We used Pearson and Spearman for size vs GC, and Kruskal–Wallis / ANOVA for GC differences across species after checking distribution assumptions.

> GC differs among species at the record level, but records are nested inside genomes, so p-values are not five independent biological replicates.

> Unique k-mer count at k=4 is saturated (256/256). Entropy is the better complexity readout. Low GC in *Dicty* lines up with low entropy.

---

## Part 6 — TATA box (1 minute)

### Plots to show
None (results are tabular). Stay on a text/table card in Gamma.

### What to say
> We scanned the first 500 bp of each concatenated assembly for the canonical TATAAA motif.

> Counts: *N. crassa* 5, *D. discoideum* 2, others 0 canonical hits in that window.

> We cannot honestly answer “what percent of genes have a TATA box?” or “do TATA genes express more?” without gene annotation and expression data.

> Missing in a 500 bp window does **not** mean missing genome-wide.

---

## Closing (30 seconds)

**Say:**
> Across GC and k-mers, *Y. lipolytica* is GC-rich, *D. discoideum* is AT-rich and compositionally distinct, and the two fungi are the closest pair. Methods are reproducible in Python; limits of sample size and TATA scanning are explicit.

**Show:** takeaways slide / thank-you.

---

## Quick “which plot for which brief question?”

| Brief question | Best plot / evidence |
|----------------|----------------------|
| Genome sizes / assemblies | Size barplot + data table |
| GC per genome; max/min | GC barplot |
| Size–GC relationship | Scatter (+ optional trend) |
| Same GC across chromosomes? | Boxplot / violin |
| K-mer diversity / entropy | Entropy bar chart |
| Closest / farthest species | Heatmap + cosine matrix |
| Clusters / outliers | PCA + dendrogram |
| Feature relationships | Correlation matrix |
| TATA differences | TATA table (no plot) |

---

## Plot file locations

**GC (Member 1):** `member1_gc_analysis/results/figures/`  
**K-mer (Member 2):** `member2_kmer_analysis/results/figures/`

**Must-show set (11 plots):**
1. `student1_gc_01_genome_gc_barplot.png`
2. `student1_gc_02_genome_size_barplot.png`
3. `student1_gc_03_size_gc_scatter.png`
4. `student1_gc_05_record_gc_boxplot.png`
5. `student1_gc_08_organism_feature_correlations.png`
6. `entropy_comparison.png`
7. `heatmap.png`
8. `similarity_matrix.png`
9. `pca.png`
10. `dendrogram.png`
11. (optional depth) `student1_gc_06_…violin…` or `student1_gc_13_main_chromosome_gc.png`
