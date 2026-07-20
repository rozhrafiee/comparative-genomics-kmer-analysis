# Student 1 GC and Genome-Size Analysis

This project contains only the Student 1 work for the comparative-genomics assignment: FASTA preparation and validation, genome characteristics, GC-content analysis, GC-related statistics, genome-size/GC visualisation, the Persian report, and the Persian oral guide. The input files are used as provided; no replacement genomes or metadata were downloaded.

## Scope and data identity

The project brief names *Caenorhabditis elegans*, *Arabidopsis thaliana*, *Schizosaccharomyces pombe*, *Drosophila melanogaster*, and *Anopheles gambiae*. Header inspection identified the five supplied files as:

| Input file | Header-derived organism | Match to brief |
|---|---|---|
| `GCA_000002985.3_WBcel235_genomic.fna` | *Caenorhabditis elegans* | Yes |
| `GCA_000182925.2_NC12_genomic.fna` | *Neurospora crassa* | No; substitution |
| `GCA_001761485.1_ASM176148v1_genomic.fna` | *Yarrowia lipolytica* | No; substitution |
| `GCA_000004695.1_dicty_2.7_genomic.fna` | *Dictyostelium discoideum* | No; substitution |
| `GCA_000710315.1_Eremothecium_coryli_genomic.fna` | *Eremothecium coryli* | No; substitution |

The mapping is derived from FASTA descriptions, not filenames alone. The exact identifiers, descriptions, file sizes, assembly labels visible in filenames, and unavailable metadata are recorded in `results/tables/student1_gc_input_file_audit.csv`.

## Project structure

```text
student1_gc_analysis/
├── data/raw/                         # read-only working copies of supplied .fna files
├── scripts/student1_gc_analysis.py   # reproducible analysis program
├── results/
│   ├── tables/                       # CSV outputs
│   ├── figures/                      # 300-DPI PNG and vector PDF figures
│   └── logs/                         # run log and QC report
├── report/                           # Persian report and oral guide
├── requirements.txt
└── README.md
```

## Requirements and installation

Python 3.10 or newer is recommended. Create an isolated environment and install the listed packages:

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The program streams FASTA records and does not load complete genomes into memory. Biopython is listed for compatibility with common FASTA workflows; the supplied script uses a streaming standard-library parser together with pandas, NumPy, SciPy, statsmodels, and Matplotlib.

## Reproduce the analysis

From the project root:

```bash
python scripts/student1_gc_analysis.py \
  --input-dir data/raw \
  --output-dir results
```

The command overwrites the generated Student 1 results. It writes a log to `results/logs/student1_gc_analysis.log` and a quality-control report to `results/logs/student1_gc_quality_control_report.txt`.

## GC definitions

The primary metric is `GC_Percent_Valid_ATGC = 100 * (G + C) / (A + T + G + C)`. `N` and other IUPAC ambiguous symbols are excluded from this denominator. `GC_Percent_All_Bases` uses total sequence length as the denominator and is reported only as a secondary diagnostic. Whole-genome GC is calculated from summed nucleotide counts, not from the arithmetic mean of record percentages. File-size MB values use decimal megabytes (bytes/1,000,000).

Record classifications are made only from header evidence: `chromosome`, `mitochondrion`, `chloroplast`, `scaffold`, `contig`, `plasmid`, `unplaced`, `unlocalized`, or `unknown`. No record is silently removed. The main-chromosome sensitivity analysis includes only records explicitly classified as `chromosome`; nuclear analysis excludes records explicitly classified as mitochondrial or chloroplast.

## Output files

### Tables

- `student1_gc_sequence_summary.csv`: one row per FASTA record with header, classification, length, A/T/G/C/N/ambiguous counts, valid length, GC counts, two GC metrics, and file size.
- `student1_gc_genome_summary.csv`: one row per header-derived organism with record counts, classification counts, total and valid lengths, whole-genome GC, ambiguity percentages, file size, and available filename-derived accession/assembly labels.
- `student1_gc_results.csv`: all-record, main-chromosome, nuclear, and organellar GC summaries, including weighted and unweighted values, median, spread, and ambiguity percentages.
- `student1_gc_statistical_results.csv`: organism-level Pearson/Spearman tests, record-level assumption diagnostics, ANOVA or Kruskal-Wallis results, and post-hoc results when applicable.
- `student1_gc_input_file_audit.csv`: file sizes, record identifiers, relevant header descriptions, header-derived organism names, possible accession/assembly strings, and QC fields.
- `student1_gc_organism_correlation_matrix.csv`: organism-level correlation matrix for genome size, GC, valid length, and ambiguity percentages.

### Figures

Each figure is saved as a 300-DPI PNG and a PDF. `student1_gc_figure_captions.txt` contains scientific captions. The figures cover whole-genome GC, genome size, size–GC scatter and trend, record-level GC distributions, main-record lengths, a compact GC-feature correlation view, ambiguity percentages, Q-Q diagnostics, record-length/GC, main-chromosome GC, and nuclear/organellar GC where available.

A heatmap or dendrogram is intentionally not generated because those analyses belong to the separate teammate scope. No code, output, or placeholder for that scope is included.

## Adding or replacing an input genome

1. Copy a new `.fna`, `.fasta`, or `.fa` file into `data/raw/` without changing its contents.
2. Rerun the command above.
3. Inspect `student1_gc_input_file_audit.csv` and the QC report before interpreting the results.

The organism name is taken from the first two alphabetic words of each FASTA description and is cross-checked for consistency across records. If the header does not provide a confident name, the output uses `Unknown organism`; it does not guess. Source database metadata is reported as `Not available in the provided file` unless explicitly present. Filename accessions and assembly labels are retained as possible metadata, not treated as proof of source database.

## Common errors and troubleshooting

- **No FASTA files found:** check `--input-dir` and file extensions.
- **Large memory use:** the parser is streaming; avoid opening multiple runs simultaneously.
- **Missing source/accession metadata:** this is a data limitation, not an analysis error; leave the explicit unavailable value in the tables.
- **Unexpected symbols:** inspect the QC report. Unexpected symbols are counted separately and excluded from valid ATGC length.
- **No main-chromosome records:** the sensitivity analysis is omitted for that organism rather than treating scaffolds as chromosomes.
- **Small-sample inference:** organism-level correlation has `n=5` available genomes. Record-level tests are descriptive and subject to pseudo-replication because records are nested within organisms.

## Scientific interpretation limits

The supplied genomes are not the exact five organisms listed in the brief, so results describe the five available header-identified genomes. A higher GC percentage can be discussed in relation to the three-versus-two hydrogen-bond distinction of G–C versus A–T pairs, but it does not prove universal whole-genome stability, gene-expression changes, or evolutionary causation. Correlation is not causation, and chromosome-level p-values should not be interpreted as five independent biological replicates.
