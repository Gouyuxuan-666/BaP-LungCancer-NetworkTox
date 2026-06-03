# ============================================================
# BaP-Lung Cancer Network Toxicology Pipeline
# Master runner: sequentially executes all pipeline steps
# Usage: Rscript 00_master.R
# ============================================================

cat("========================================\n")
cat("BaP-Lung Cancer Network Toxicology\n")
cat("Pipeline started:", format(Sys.time()), "\n")
cat("========================================\n\n")

# --- Configuration ---
DATA_DIR  <- "data"
OUT_DIR   <- "output"
FIG_DIR   <- "figures"

dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(FIG_DIR, showWarnings = FALSE)

# --- Step 1: Data Preparation ---
cat("[Step 1/7] Data preparation: merging GEO datasets, ComBat correction...\n")
source("scripts/01_ML_data_prep.R")
cat("  -> data.train.txt, data.test.txt generated\n\n")

# --- Step 2: ML Model Construction (113 models) ---
cat("[Step 2/7] Building 113 machine learning models...\n")
source("scripts/02_refer_ML_functions.R")
source("scripts/03_ML_113_models.R")
cat("  -> model.AUCmatrix.txt, model.genes.txt generated\n\n")

# --- Step 3: SHAP Analysis ---
cat("[Step 3/7] SHAP feature importance analysis...\n")
source("scripts/04_SHAP_analysis.R")
cat("  -> SHAP plots generated\n\n")

# --- Step 4: Volcano Plot & ROC Curves ---
cat("[Step 4/7] Volcano plot and ROC curves for hub genes...\n")
source("scripts/05_volcano_ROC.R")
cat("  -> vol.pdf, ROC.genes.pdf generated\n\n")

# --- Step 5: PPI Network ---
cat("[Step 5/7] STRING PPI network analysis...\n")
source("scripts/07_PPI_network.R")
cat("  -> PPI network visualized\n\n")

# --- Step 6: CIBERSORT Immune Infiltration (optional, needs raw intensity data) ---
cat("[Step 6/7] CIBERSORT immune infiltration (skipped if raw data unavailable)...\n")
tryCatch({
  source("scripts/06_CIBERSORT_pipeline.R")
  cat("  -> CIBERSORT results generated\n\n")
}, error = function(e) {
  cat("  -> Skipped: raw microarray intensity data required\n\n")
})

# --- Done ---
cat("========================================\n")
cat("Pipeline completed:", format(Sys.time()), "\n")
cat("Results in:", OUT_DIR, "\n")
cat("Figures in:", FIG_DIR, "\n")
cat("========================================\n")
