# BaP-Lung Cancer Network Toxicology Pipeline

**Identification of Candidate Molecular Mediators in Benzo[a]pyrene-Associated Lung Cancer through Network Toxicology, Ensemble Machine Learning, and SHAP Analysis**

## Overview

This repository contains the complete computational pipeline for our study identifying ADH1B as a candidate mediator of BaP-associated lung carcinogenesis. The pipeline integrates:

- Network toxicology (ChEMBL, SEA, SwissTargetPrediction)
- WGCNA co-expression analysis
- 113 machine learning models (12 base thms, pairwislty， combinations)
- SHAP feature importance analysis
- Molecular docking (AutoDock Vina)
- Molecular dynamics simulation (GROMACS 2022, 200 ns)
- CIBERSORT immune infiltration analysis

## Pipeline Steps

| Script | Step | Description |
|--------|------|-------------|
| `01_ML_data_prep.R` | Data Preparation | Merge 6 GEO datasets, ComBat batch correction, train/test split |
| `02_refer_ML_functions.R` | ML Framework | Core ML functions (Lasso, Ridge, Enet, Stepglm, SVM, glmBoost, RF, GBM, LDA, XGBoost, Naive Bayes) |
| `03_ML_113_models.R` | Model Training | Build and evaluate 113 pairwise-combination models, select optimal |
| `04_SHAP_analysis.R` | SHAP Analysis | Compute SHAP values, generate summary/beeswarm/dependence plots |
| `05_volcano_ROC.R` | Visualization | Volcano plots for DEGs, ROC curves for hub genes |
| `06_CIBERSORT_pipeline.R` | Immune Infiltration | CIBERSORT deconvolution, barplot, differential boxplot, correlation heatmap |
| `07_PPI_network.R` | Network Analysis | STRING PPI network construction and Cytoscape visualization |

## Data Sources

Six GEO datasets: GSE10072, GSE115002, GSE19804, GSE32863, GSE43458, GSE68465

BaP target databases: ChEMBL, SEA, SwissTargetPrediction

## Quick Start

```bash
git clone https://github.com/YOUR_USERNAME/BaP-LungCancer-NetworkTox.git
cd BaP-LungCancer-NetworkTox
Rscript scripts/00_master.R
```

Data files are downloaded automatically from NCBI GEO (see `data/README.md`).

## Requirements

- **R >= 4.2.1**
- **R packages**: limma, sva, WGCNA, clusterProfiler, org.Hs.eg.db, enrichplot, pheatmap, ggplot2, VennDiagram, glmnet, gbm, xgboost, randomForest, e1071, mboost, plsRglm, caret, IOBR, reshape2, ggpubr, corrplot
- **GROMACS 2022** (molecular dynamics simulation; optional — pre-computed results included)
- **AutoDock Vina 1.5.7** (molecular docking; optional — pre-computed results included)
- **Cytoscape 3.10.3** (network visualization)

## Key Results

- **522** DEGs (163 up, 359 down) from 6 merged GEO datasets (227 samples)
- **30** candidate mediators (BaP targets ∩ DEGs ∩ WGCNA turquoise module)
- **113** machine learning models evaluated; Lasso + Stepglm [both] selected
- **8** hub candidates: ADH1B, ADRB1, CDKN1A, GMNN, HTR2B, MMP1, NET1, NPY1R
- **ADH1B**: top SHAP value (0.190), strongest BaP docking score (-6.6 kcal/mol), stable 200-ns MD complex

## Citation

Manuscript submitted for publication. If you use this code, please cite:本代码原作者为勾雨轩,编写此代码时就读于新疆第二医学院麻醉学本科 大学一年级  (Xinjiang Second Medicial College)
> [AUTHOR :勾雨轩 ] Identification of Candidate Molecular Mediators in Benzo[a]pyrene-Associated Lung Cancer through Network Toxicology, Ensemble Machine Learning, and SHAP Analysis. (2026)如需使用本论文代码或需进行复现,请您确保在代码编写者同意之后再进行复现或改写且不得用于商业用途(如要用于商业用途请联系作者)

## License

MIT
