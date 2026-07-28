# Brief sections 3-6 — direct answers

## بخش سوم — K-mer

Most frequent 3-mer per organism:
- Caenorhabditis elegans: TTT (n=6287034)
- Neurospora crassa: TTT (n=875051)
- Yarrowia lipolytica: CAA (n=416599)
- Dictyostelium discoideum: AAA (n=3131027)
- Eremothecium coryli: AAA (n=254667)
Most frequent 4-mer per organism:
- Caenorhabditis elegans: AAAA (n=2977271)
- Neurospora crassa: TTTT (n=280073)
- Yarrowia lipolytica: TTTT (n=139253)
- Dictyostelium discoideum: AAAA (n=1845792)
- Eremothecium coryli: AAAA (n=80029)
Most frequent 5-mer per organism:
- Caenorhabditis elegans: AAAAA (n=1354147)
- Neurospora crassa: TTTTT (n=106524)
- Yarrowia lipolytica: AAAAA (n=53846)
- Dictyostelium discoideum: AAAAA (n=1247746)
- Eremothecium coryli: AAAAA (n=25068)

Unique k-mer counts (k=4):
- Caenorhabditis elegans: 256
- Neurospora crassa: 256
- Yarrowia lipolytica: 256
- Dictyostelium discoideum: 256
- Eremothecium coryli: 256

Species-exclusive k-mers: none found

Highest k-mer diversity (k=4): Eremothecium coryli (2.83516e-05)

## بخش چهارم — Similarity

Closest species (Euclidean): Neurospora crassa - Yarrowia lipolytica (0.015417)
Farthest species (Euclidean): Yarrowia lipolytica - Dictyostelium discoideum (0.102305)
Animal clustering: cannot be tested — only one animal organism (C. elegans) is in this dataset.
Plant vs animal separation: cannot be tested — no plant organism is in this dataset.

## بخش ششم — TATA Box

Canonical TATAAA in first 500 bp of pooled genome (not gene-annotation based):
- Caenorhabditis elegans: 0 match(es)
- Neurospora crassa: 5 match(es)
- Yarrowia lipolytica: 0 match(es)
- Dictyostelium discoideum: 2 match(es)
- Eremothecium coryli: 0 match(es)
% of genes with a TATA box: cannot be answered — no gene/promoter annotation (GFF/GTF) was used, only a raw-sequence scan.
TATA+ genes having higher expression: cannot be answered — no expression data was collected in this project.

## بخش پنجم — Statistics

- Correlation genome_size vs kmer_diversity: r=-0.7428574247967478, p=0.15034542863238057 (low power (small n))
- Correlation gc_percent vs genome_size: r=-0.24863063315958012, p=0.6867261177155706 (low power (small n))
- Correlation gc_percent vs shannon_entropy: r=0.966094518363889, p=0.007456200687136317 (low power (small n))
- Correlation tata_box_count vs gc_percent: r=0.11158921529669427, p=0.858215618111856 (low power (small n))
- ANOVA shannon_entropy: F=19.80908387399981, p=1.0865033261035296e-11
- ANOVA kmer_diversity: F=13.375484121159129, p=1.415747767518784e-08
