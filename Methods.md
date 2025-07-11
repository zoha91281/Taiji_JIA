## Methods and Workflow

### Data Acquisition and Preprocessing

Multi-omic datasets were obtained from the NCBI Gene Expression Omnibus (GEO) under accession number **[GSE164213](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164213)**. These included peripheral blood CD4+ T cells from children with active juvenile idiopathic arthritis (JIA), covering:

- **ATAC-seq**: to measure chromatin accessibility
- **RNA-seq**: for gene expression profiling
- **HiChIP**: for 3D chromatin architecture and regulatory looping

#### Step-by-Step:

1. **ATAC-seq Processing**
   - Download `bedGraph` files from GEO.
   - Convert them to `narrowPeak` format using [MACS3 v3.0.0a6](https://github.com/macs3-project/MACS):
     - Default MACS3 parameters were used to model local Poisson noise and detect statistically significant open chromatin regions.

2. **RNA-seq Processing**
   - Retrieve raw count matrices (`.txt` or `.csv`) from GEO.
   - Use `biomaRt` (R package) to map Ensembl transcript IDs to gene symbols:
     ```r
     library(biomaRt)
     ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
     genes <- getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol"), mart = ensembl)
     ```
   - Additional formatting and conversion were handled with:
     - `readr::read_csv()` for fast file import
     - `plyr::mapvalues()` to align Ensembl IDs to gene symbols
   - Final output: standardized **gene-level expression matrices** for each sample.

3. **HiChIP / 3D Contacts**
   - Download HiChIP output or processed contacts from GEO.
   - Run [EpiTensor](https://github.com/gersteinlab/EpiTensor) to integrate epigenomic marks and infer enhancer–promoter interactions.
   - Filter to retain the **top 10% most confident long-range contacts** to prioritize meaningful enhancer–gene loops.

---

### Regulatory Network Construction with Taiji

We used the [**Taiji framework**](https://github.com/OSU-BMBL/taiji) to integrate ATAC-seq, RNA-seq, and HiChIP data into a unified transcriptional regulatory network. The process:

- **Inputs:**
  - ATAC-seq peak files (`narrowPeak`)
  - Gene expression matrices (RNA-seq)
  - Chromatin contact predictions (HiChIP/EpiTensor)
  - Custom Taiji YAML config (see [`taiji_config/config.yaml`](config.yaml))

- **Execution:**
  - Taiji maps TF motifs to accessible regions and links them to gene targets using 3D contact data.
  - It calculates **influence scores** using the PageRank algorithm.
  - Run time: ~5 days on an HPC cluster.

---

### Pathway Enrichment and Network Visualization

1. **Pathway Analysis**
   - Top gene targets from TF–gene networks were analyzed using:
     - `ReactomePA::enrichPathway()`
     - `ClusterProfiler::enrichGO()`
   - Enriched pathways included cytokine signaling, cell cycle regulation, and immune cell migration.

2. **Heatmap Visualization**
   - Influence scores (from Taiji) were z-score scaled and visualized using `pheatmap`.
   - Heatmaps were stratified by:
     - All samples
     - TNF-treated only
     - Control + medium-treated

3. **Network Graphs**
   - Regulatory modules were visualized with `igraph` to show hierarchical TF–target relationships.

---

### Statistical Analysis

To explore **tissue-specific transcriptional patterns**, we applied:

- **Pairwise Wilcoxon rank-sum tests** on TF influence scores.
- Focused on healthy control samples to compare:
  - Hand vs. Knee
  - Hand vs. Hip
  - Hip vs. Knee

Script: [`scripts/stats/wilcoxon_tests.R`](scripts/stats/wilcoxon_tests.R)

---

### Computational Environment and Reproducibility

All analyses were performed in a reproducible, Conda-managed R environment:

#### Conda Setup

```bash
conda create -n r_bio_env -c bioconda -c conda-forge \
  r-base=4.2.2 \
  bioconductor-reactomepa \
  bioconductor-clusterprofiler \
  -y
conda activate r_bio_env
install.packages("readr")
install.packages("plyr")
install.packages("pheatmap")
install.packages("igraph")
