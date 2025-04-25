# Adrenal Medulla Single-Cell Analysis

This repository contains scripts for a single-cell transcriptomic analysis of human adrenal medulla (AM) cells.  
The goal was to identify gene expression changes along differentiation trajectories using Slingshot tool.

---

## Contents

| File                   | Description |
|------------------------|-------------|
| `full_preprocessing.R` | Full preprocessing pipeline for the complete dataset and extraction of adrenal medulla cells|
| `scrublet.ipynb`       | Python notebook using **Scrublet** for doublet detection on raw data. |
| `am_processing.R`      | Processing of adrenal medulla single-cell data, filtering artefact cells and integration. |
| `trajectories.R`       | Differentiation trajectories analysis using **Slingshot**, or **tradeSeq**, and identification of pseudotime-dependent gene expression. |
| `helpers.R`            | Custom functions used in analysis. |


---

## Abstract

<div align="justify">

Neuroblastoma is one of the most common pediatric extracranial tumors that is found
adrenal medulla. Evidence suggests that neuroblastoma may arise during developmental stages
from one of the cell types originating in the developing adrenal medulla. This particular
analysis focused on the single-cell RNA-sequencing data of the development of human adrenal
glands. During this study, different cell types that populate developing adrenal glands were
identified. Specifically, analysis was focused on 3 cell types originating in the adrenal medulla:
Schwann cell precursors (SCPs), sympathoblasts and chromaffin cells. Utilizing trajectory
analysis with Slingshot, specific genesets associated with transitions between cell states. The
transition from SCPs to chromaffin cells was marked by strong associations with *DLK1*, *PENK*,
*CDKN1C*, and *PLP1*. *CD24*, *ELAVL4*, *SPARC*, and *METRN* were specifically linked to the
trajectory of SCPs and sympathoblasts. Finally, the trajectory from sympathoblasts to
chromaffin cells involved genes such as *CHGA*, *HTATSF1*, *JUNB*, and *PLXNA4*.