# A Unified Computational Framework for Bulk RNA-seq Analysis in Colorectal Cancer

[![R](https://img.shields.io/badge/Language-R_4.0+-198CE7.svg)](https://www.r-project.org/)
[![DESeq2](https://img.shields.io/badge/Bioc-DESeq2-F05032.svg)](https://bioconductor.org/packages/DESeq2/)
[![WGCNA](https://img.shields.io/badge/CRAN-WGCNA-00599C.svg)](https://cran.r-project.org/package=WGCNA)
[![sva](https://img.shields.io/badge/Bioc-sva_(ComBat--seq)-8A2BE2.svg)](https://bioconductor.org/packages/sva/)
[![org.Hs.eg.db](https://img.shields.io/badge/Bioc-org.Hs.eg.db-3CB371.svg)](https://bioconductor.org/packages/org.Hs.eg.db/)
[![clusterProfiler](https://img.shields.io/badge/Bioc-clusterProfiler-FF1493.svg)](https://bioconductor.org/packages/clusterProfiler/)
[![enrichplot](https://img.shields.io/badge/Bioc-enrichplot-FF8C00.svg)](https://bioconductor.org/packages/enrichplot/)
[![ComplexHeatmap](https://img.shields.io/badge/Bioc-ComplexHeatmap-DC143C.svg)](https://bioconductor.org/packages/ComplexHeatmap/)
[![biomaRt](https://img.shields.io/badge/Bioc-biomaRt-00CED1.svg)](https://bioconductor.org/packages/biomaRt/)
[![BiocParallel](https://img.shields.io/badge/Bioc-BiocParallel-4682B4.svg)](https://bioconductor.org/packages/BiocParallel/)
[![tidyverse](https://img.shields.io/badge/CRAN-tidyverse-5F9EA0.svg)](https://cran.r-project.org/package=tidyverse)
[![data.table](https://img.shields.io/badge/CRAN-data.table-D2691E.svg)](https://cran.r-project.org/package=data.table)
[![circlize](https://img.shields.io/badge/CRAN-circlize-9932CC.svg)](https://cran.r-project.org/package=circlize)
[![RColorBrewer](https://img.shields.io/badge/CRAN-RColorBrewer-2E8B57.svg)](https://cran.r-project.org/package=RColorBrewer)
[![ggrepel](https://img.shields.io/badge/CRAN-ggrepel-BDB76B.svg)](https://cran.r-project.org/package=ggrepel)
[![scales](https://img.shields.io/badge/CRAN-scales-CD5C5C.svg)](https://cran.r-project.org/package=scales)
[![patchwork](https://img.shields.io/badge/CRAN-patchwork-4CAF50.svg)](https://cran.r-project.org/package=patchwork)
[![ggExtra](https://img.shields.io/badge/CRAN-ggExtra-8B4513.svg)](https://cran.r-project.org/package=ggExtra)
[![ggplotify](https://img.shields.io/badge/CRAN-ggplotify-20B2AA.svg)](https://cran.r-project.org/package=ggplotify)
[![svglite](https://img.shields.io/badge/CRAN-svglite-FF6347.svg)](https://cran.r-project.org/package=svglite)
[![License](https://img.shields.io/badge/License-MIT-4CAF50.svg)](https://opensource.org/licenses/MIT)

## Project Overview
This repository contains a comprehensive and highly rigorous computational pipeline designed for the transcriptomic evaluation of Colorectal Cancer (CRC). Colorectal cancer represents a highly prevalent gastrointestinal malignancy characterized by complex underlying molecular mechanisms and profound tumor heterogeneity. 

A major challenge in transcriptomic meta-analyses is the limited statistical power of individual studies and the technical batch effects that arise when combining datasets across different sequencing platforms or institutions. By systematically integrating independent bulk RNA-seq datasets sourced from the NCBI Gene Expression Omnibus (GEO), this project aims to identify robust, highly confident transcriptomic signatures and co-expression modules associated with CRC tumorigenesis. 

The core of this project relies on a unified R script written to process raw count matrices, apply negative binomial regression to eliminate cross-study technical batch effects, perform differential expression analysis, and execute systems-level weighted gene co-expression network analysis (WGCNA).

## Datasets Analyzed
To ensure maximum biological consistency and statistical power, this study utilizes a strictly curated meta-cohort comprising **54 individual samples** (27 normal adjacent/control tissues and 27 primary CRC tumors). 

* **GSE144259:** 6 samples (3 Normal, 3 Tumor). 
  * *Exclusion Criteria:* 3 metastatic samples (CRC1M, CRC2M, CRC3M) were explicitly excluded to focus solely on primary tumor biology.
* **GSE50760:** 36 samples (18 Normal, 18 Tumor). 
  * *Exclusion Criteria:* 18 metastatic samples designated by the suffix '-3' were systematically excluded prior to count matrix integration.
* **GSE87096:** 12 samples (6 Normal, 6 Tumor). 
  * *Exclusion Criteria:* Filtered to retain only standard RNA-seq samples; 24 MeDIP/hMeDIP DNA methylation profiling samples were excluded to maintain uniform data modalities.

## Analytical Methodology & Core Pipeline

The entire end-to-end computational workflow is executed via a single, highly modular R script located at **`r script/colaon_bulk_rna.R`**. 

### 1. Preprocessing & Technical Batch Correction (ComBat-seq)
Raw counts from the three independent GEO cohorts were filtered for common genes and merged into a single combined count matrix. Cross-study technical variance was mitigated directly at the raw count level using `ComBat_seq`. Unlike standard ComBat, which operates on normalized data, ComBat-seq utilizes a negative binomial regression model to adjust batch effects while preserving the integer nature of the RNA-seq counts, maintaining the strict statistical assumptions required for downstream modeling.

### 2. Differential Expression Modeling (DESeq2)
Differential Gene Expression (DGE) modeling was performed on the combined count matrix using `DESeq2`, contrasting the primary Tumor condition against the Normal control. Biotype filtering via custom regex matching was applied to remove non-coding elements (LOC, MIR, snoRNAs, LINC), isolating functionally relevant protein-coding targets.

### 3. Weighted Gene Co-Expression Network Analysis (WGCNA)
To move beyond single-gene DGE, systems-level transcriptomic architecture was modeled using `WGCNA`. Variance-stabilized (VST) counts, filtered for the top 5000 highly variable genes, were used to construct a scale-free topological network (soft threshold power = 14). Genes were clustered into co-expression modules, which were subsequently correlated with the clinical trait (Diseased vs. Control) to identify primary driver modules in CRC pathogenesis.

### 4. Functional Enrichment & ID Translation
Significant protein-coding transcripts and hub genes were subjected to over-representation analysis (ORA) across Gene Ontology (BP, CC, MF) and KEGG Pathways using `clusterProfiler`. A dedicated annotation module utilizing `biomaRt` was implemented to accurately map complex STRING database ENSP IDs back to universal Gene Symbols and Entrez IDs for network continuity.

---

## Visualizations

### Quality Control: Principal Component Analysis (PCA)
The PCA visualizations validate the efficacy of the ComBat-seq batch correction. The 54 samples successfully cluster according to biological condition (PC1) rather than their respective GEO source study.

<p align="center">
  <img src="results_CRC/03_plots/QC_PCA/PCA_by_Batch.png" width="48%" alt="PCA by Batch">
  <img src="results_CRC/03_plots/QC_PCA/PCA_by_Condition.png" width="48%" alt="PCA by Condition">
</p>

### Differential Expression Landscape
Global shifts in CRC gene expression are mapped via a Volcano plot, accompanied by a hierarchical clustering Heatmap highlighting the consistent, normalized expression patterns of the top 50 highly significant differentially expressed genes across the meta-cohort.

<p align="center">
  <img src="results_CRC/03_plots/DEG/Volcano_plot.png" width="48%" alt="Volcano Plot">
  <img src="results_CRC/03_plots/Heatmap/Heatmap_top25up_top25dn.png" width="48%" alt="Top 50 DEGs Heatmap">
</p>

### Functional Enrichment (GO & KEGG)
Elucidation of functional roles across Gene Ontology categories and KEGG pathways driving the transcriptomic signature, compiled into a single high-density collage.

<table align="center" style="width:100%; border:none;">
  <tr>
    <td width="50%"><img src="results_CRC/03_plots/enrichment/png/GO_BP_CRC.png" alt="GO Biological Process" width="100%"></td>
    <td width="50%"><img src="results_CRC/03_plots/enrichment/png/GO_CC_CRC.png" alt="GO Cellular Component" width="100%"></td>
  </tr>
  <tr>
    <td width="50%"><img src="results_CRC/03_plots/enrichment/png/GO_MF_CRC.png" alt="GO Molecular Function" width="100%"></td>
    <td width="50%"><img src="results_CRC/03_plots/enrichment/png/KEGG_CRC.png" alt="KEGG Pathways" width="100%"></td>
  </tr>
</table>

### Systems Biology: WGCNA
Scale-free topology determination, hierarchical clustering, and module-trait relationship heatmaps identifying highly correlated gene clusters within the CRC meta-cohort.

<p align="center">
  <img src="results_CRC/03_plots/WGCNA/png/Fig1_Network_Topology.png" width="100%" alt="WGCNA Network Topology">
</p>
<p align="center">
  <img src="results_CRC/03_plots/WGCNA/png/Fig2_Module_Trait_Scatters.png" width="100%" alt="WGCNA Module Trait Relationships">
</p>

---

## Repository Structure
* **`r script/colaon_bulk_rna.R`** : Core R script executing the full analytical workflow from raw data ingestion to visualization.
* `count_matrix_metadata/` : Uncompressed `.tsv.gz` raw data files and merged clinical metadata.
* `results_CRC/01_objects/` : Serialized R objects (`.rds`) for rapid local reloading of DESeq2 and WGCNA environments.
* `results_CRC/02_tables/` : Compiled DEG outputs, module assignments, raw/adjusted count matrices, and ORA results.
* `results_CRC/03_plots/` : High-fidelity graphical outputs. *Note: Heavy, uncompressed `.tiff` files generated for publication are excluded from version control via `.gitignore`.*

---
