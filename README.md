# taiji_sri

**Integrative Transcription Factor Network Analysis in Juvenile Idiopathic Arthritis (JIA)**

This repository contains code, configuration, and processed data for analyzing transcription factor (TF) influence in CD4+ T cells from JIA patients using the Taiji systems biology framework. By integrating ATAC-seq, RNA-seq, and Hi-C (via EpiTensor), we construct gene regulatory networks and apply PageRank to prioritize key TFs contributing to disease.

## 🔬 Project Summary

Juvenile Idiopathic Arthritis (JIA) is a heterogeneous pediatric autoimmune disease. We hypothesize that CD4+ T cells from JIA patients are primed for inflammation via altered TF-driven regulatory programs. This project uses Taiji to reconstruct tissue- and condition-specific TF–gene networks and statistically compare TF influence scores across disease states and anatomical sites.

## 🧬 Data Sources

- **GSE164213** (GEO): ATAC-seq, RNA-seq, and HiChIP from peripheral blood CD4+ T cells.
- All data is processed for basal (unstimulated) conditions unless otherwise stated.
