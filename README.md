# Taiji_JIA

**Integrative Transcription Factor Network Analysis in Juvenile Idiopathic Arthritis (JIA)**

This repository contains code, configuration, and processed data for analyzing transcription factor (TF) influence in CD4+ T cells from JIA patients using the Taiji systems biology framework. By integrating ATAC-seq, RNA-seq, and Hi-C (via EpiTensor), we construct gene regulatory networks and apply PageRank to prioritize key TFs contributing to disease.

---

## Project Summary

Juvenile Idiopathic Arthritis (JIA) is a heterogeneous pediatric autoimmune disease affecting children under 16. We hypothesize that aberrant transcription factor (TF) regulation in CD4+ T cells contributes to disease persistence and immune priming. Using a multi-omic approach and the Taiji framework, we model gene regulatory networks under basal conditions to uncover TFs that differentiate disease states and tissue-specific signatures.

---

## Data Sources

- **GSE164213** (GEO): ATAC-seq, RNA-seq, and HiChIP profiles from peripheral blood CD4+ T cells from children with active JIA.
- Processed data includes:
  - Normalized gene expression (RNA-seq)
  - Chromatin accessibility peaks (ATAC-seq narrowPeaks)
  - High-confidence 3D enhancer–promoter contacts (HiChIP via EpiTensor)
- Metadata includes sample IDs, anatomical site, and treatment (control, medium, TNF).

---

## 🧠 Methodology Overview

1. **Data Preprocessing**
   - Raw count matrices and peak files downloaded and converted to compatible formats.
   - Metadata curated and matched across assays.

2. **Taiji Network Construction**
   - TF–gene networks constructed by combining:
     - TF binding to accessible chromatin (from ATAC-seq)
     - Gene expression level (from RNA-seq)
     - Spatial enhancer–promoter contact (from HiChIP/EpiTensor)
   - Influence scores for each TF calculated using the PageRank algorithm.

3. **TF Prioritization and Visualization**
   - Top 250 most variable TFs selected based on variance across samples.
   - Z-score scaling applied and heatmaps generated for each condition group (all, TNF only, control/medium).

4. **Statistical Comparison**
   - Pairwise Wilcoxon rank-sum tests performed on TF influence scores across tissue sites (hand, hip, knee) in control samples.
