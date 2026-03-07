# GPCR Receptor Expression Across Evolution

> A comparative single-cell RNA-seq study of G protein-coupled receptor (GPCR) expression across three invertebrate species — *Aedes aegypti*, *Camponotus floridanus* (Ant), and *Drosophila melanogaster* — with orthogroup-based cross-species comparisons and phylogenetic analysis.

**MSc Bioinformatics Dissertation Project** | University of Leicester | 2023  
**Author:** Abhijith Suresh Babu | Pharm.D, MSc Bioinformatics

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Species Studied](#species-studied)
- [Workflow](#workflow)
- [Key Analyses](#key-analyses)
- [Tools & Dependencies](#tools--dependencies)
- [Results Summary](#results-summary)
- [Data](#data)
- [How to Run](#how-to-run)

---

## Overview

GPCRs are the largest family of membrane receptors and key drug targets. This project investigates how GPCR expression is conserved or diverged across evolutionary lineages by applying single-cell RNA sequencing (scRNA-seq) to three arthropod species. The pipeline integrates:

- scRNA-seq clustering and cell type annotation
- GPCR gene list filtering per species
- Orthogroup (OG) mapping for cross-species homology
- Dot plots and feature plots for receptor subtype expression (serotonin, dopamine)
- Phylogenetic tree visualisation using FigTree
- Cross-species GPCR expression comparison

---

## Repository Structure

```
GPCR-Receptor-Expression-Across-Evolution/
│
├── Aedes/                          # Aedes aegypti analysis
│   ├── Aedes.R                     # Main R analysis script
│   ├── data/                       # Filtered/unfiltered GPCR lists, OG-mapped TSVs, XLSX
│   │   ├── Aedes_aegypti_GPCR_list.tsv
│   │   ├── Aedes_GPcr_Filtered.tsv
│   │   ├── Aedes_GPcr_Filtered_with OG.tsv
│   │   ├── Aedes_MarkerGenes.tsv
│   │   ├── Aedes_dopamine_.tsv
│   │   ├── Aedes_serotonin_.tsv
│   │   └── ...
│   ├── markers/                    # Per-cluster marker gene files
│   └── results/                    # Cluster plots, dot plots, UMAP outputs
│
├── Ant/                            # Camponotus floridanus analysis
│   ├── data/                       # GPCR filtered lists, OG-mapped TSVs, marker genes
│   │   ├── ant_25KGPcr_filtered.tsv
│   │   ├── ant_25K_GPCRfiltered_OG.tsv
│   │   ├── ant_25K_MarkerGenes.tsv
│   │   ├── ant_25K_dopamine_.tsv
│   │   ├── ant_25K_serotonin_.tsv
│   │   └── ...
│   ├── markers/                    # Per-cluster marker gene files
│   └── results/                    # Cluster plots, dot plots, scatter plots
│
├── Feature-Plots/                  # Cross-species GPCR feature plots
│   ├── Aedes/
│   ├── Ant/
│   └── Drosophila/
│
├── Phylogeny/                      # FigTree phylogenetic tree screenshots
│   └── Screenshot at 2023-07-13 *.png
│
├── data/                           # Cross-species integrated data
│   ├── cross-species/              # XLSX files for Aedes, Ant, Drosophila OG mapping
│   │   ├── Aedes GPCR with OG.xlsx
│   │   ├── ANT GPCR OG.xlsx
│   │   ├── Droso GPCR OG.xlsx
│   │   ├── GPCR EXPRESSION.xlsx
│   │   ├── Markers clusters.xlsx
│   │   └── ...
│   └── GRaph.py                   # Python script for cross-species graph generation
│
├── GPCR_Receptor_Expression_Evolution.pdf   # Full dissertation report
└── README.md
```

---

## Species Studied

| Species | Common Name | Relevance |
|---|---|---|
| *Aedes aegypti* | Yellow fever mosquito | Vector biology; neuromodulation in host-seeking |
| *Camponotus floridanus* | Florida carpenter ant | Social insect; complex behaviours mediated by GPCRs |
| *Drosophila melanogaster* | Fruit fly | Model organism; well-annotated GPCR repertoire |

---

## Workflow

```
Raw scRNA-seq Data
        │
        ▼
 Seurat / R Pipeline (Aedes.R)
        │
        ├── QC & Filtering
        ├── Normalisation & Scaling
        ├── PCA → UMAP → Clustering
        ├── Marker Gene Identification
        │
        ▼
 GPCR Gene List Filtering
        │
        ├── Species-specific GPCR lists (TSV)
        ├── Orthogroup (OG) Mapping via OrthoFinder
        │
        ▼
 Expression Analysis
        │
        ├── Dot Plots (serotonin, dopamine receptors per cluster)
        ├── Feature Plots (spatial expression on UMAP)
        ├── Scatter Plots
        │
        ▼
 Cross-Species Comparison
        │
        ├── OG-based homology mapping
        ├── GPCR EXPRESSION.xlsx (integrated)
        ├── GRaph.py (visualisation)
        │
        ▼
 Phylogenetic Analysis
        │
        └── FigTree visualisation of GPCR evolutionary relationships
```

---

## Key Analyses

### 1. scRNA-seq Clustering
Each species' dataset was processed through a standard Seurat pipeline: quality filtering, normalisation, dimensionality reduction (PCA + UMAP), and unsupervised clustering. Clusters were annotated based on marker genes.

### 2. GPCR Filtering
A curated GPCR gene list for each species was cross-referenced against the scRNA-seq expression matrix. Genes were filtered for expression thresholds to retain biologically relevant receptors.

### 3. Orthogroup Mapping
GPCRs from each species were mapped to orthogroups using OrthoFinder output, enabling direct comparison of homologous receptors across species regardless of naming differences.

### 4. Receptor Subtype Expression
Serotonin and dopamine receptor subtypes were specifically investigated across clusters using dot plots (showing both expression level and percentage of expressing cells) and feature plots (showing spatial distribution on UMAP).

### 5. Phylogenetic Analysis
FigTree was used to visualise phylogenetic relationships among GPCR sequences, with annotated trees showing evolutionary conservation and divergence across the three species.

---

## Tools & Dependencies

### R Packages
```r
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(patchwork)
```

### Python
```python
# GRaph.py dependencies
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
```

### External Tools
- **OrthoFinder** — orthogroup inference for cross-species mapping
- **FigTree** — phylogenetic tree visualisation
- **Excel / LibreOffice** — cluster annotation (ODT files)

---

## Results Summary

- Identified distinct GPCR expression profiles across cell clusters in all three species
- Serotonin and dopamine receptors showed **cluster-specific expression patterns**, suggesting cell-type-specific neuromodulatory roles
- Orthogroup mapping revealed **conserved GPCR families** across Aedes, Ant, and Drosophila, with species-specific expansions in certain receptor classes
- Phylogenetic analysis confirmed evolutionary conservation of core GPCR subfamilies in arthropods

---

## Data

Raw scRNA-seq datasets were sourced from publicly available repositories. Processed outputs (filtered TSVs, marker gene tables, OG-mapped files) are available in the species-specific `data/` folders. The full cross-species integrated dataset is in `data/cross-species/`.

> **Note:** Large raw data files (e.g. `.rds` Seurat objects) are not included due to file size constraints. Processed TSV and XLSX outputs are provided for reproducibility.

---

## How to Run

### R Analysis (Aedes example)
```bash
# Clone the repository
git clone https://github.com/abhijithsureshbabu97/GPCR-Receptor-Expression-Across-Evolution.git
cd GPCR-Receptor-Expression-Across-Evolution

# Open R and run the Aedes script
Rscript Aedes/Aedes.R
```

### Python Cross-Species Graph
```bash
cd data/
python GRaph.py
```

---

## Citation

If you use this code or data, please cite:

> Suresh Babu, A. (2023). *GPCR Receptor Expression Across Evolution: A comparative scRNA-seq analysis in Aedes aegypti, Camponotus floridanus, and Drosophila melanogaster.* MSc Bioinformatics Dissertation, University of Leicester.

---

## Contact

**Abhijith Suresh Babu**  
MSc Bioinformatics, University of Leicester  
Pharm.D, RVS College of Pharmacy  
GitHub: [@abhijithsureshbabu97](https://github.com/abhijithsureshbabu97)
