################################################################################
# config.R — Pseudobulk DESeq2 Pipeline Configuration
#
# Edit this file for each new project. Do not edit pipeline.R or modules/.
#
# Usage:
#   Rscript pipeline.R config.R
################################################################################

# ── Project ───────────────────────────────────────────────────────────────────

PROJECT_NAME <- "my_project"          # prefix for all output files
OUTPUT_DIR   <- "."                   # root output directory

# ── Input ─────────────────────────────────────────────────────────────────────

# Path to annotated Seurat object.
# Preferred format: .rds (single object, no variable name needed)
#   saveRDS(your_seurat_obj, "object.rds")
# Also accepts: .RData / .rda (must contain exactly one Seurat object)
SEURAT_PATH <- "/path/to/seurat_object.rds"

# ── Cell type and sample metadata columns ─────────────────────────────────────

CELLTYPE_COL <- "celltype.level2"  # colData column with cell type labels
SAMPLE_COL   <- "Library"          # colData column with sample/library IDs
GROUP_COL    <- "Age"              # colData column used as group in prepSCE

# ── Pseudobulk QC ─────────────────────────────────────────────────────────────

# Minimum cells per sample for a cell type to be included
MIN_CELLS_PER_SAMPLE  <- 15

# Minimum number of samples meeting the cell threshold
MIN_SAMPLES_PASSING   <- 4

# Minimum number of samples a gene must be detected in (count > 0)
MIN_SAMPLES_EXPRESSED <- 3

# ── DESeq2 design ─────────────────────────────────────────────────────────────

# Biological variables of interest — DESeq2 results are extracted for each one.
# All variables are included in the design formula together so they mutually
# adjust for each other.
#
# BIO_VARS is a named list. Each entry has:
#   type      — "numeric" (continuous) or "categorical" (factor)
#   reference — for categorical only: the baseline level (e.g. "CTL")
#               ignored for numeric variables
#
# Examples:
#   Single numeric variable (age):
#     BIO_VARS <- list(age = list(type = "numeric"))
#
#   Single categorical variable (treatment):
#     BIO_VARS <- list(treatment = list(type = "categorical", reference = "CTL"))
#
#   Multiple biological variables (age + sex):
#     BIO_VARS <- list(
#       age = list(type = "numeric"),
#       sex = list(type = "categorical", reference = "F")
#     )
BIO_VARS <- list(
  age = list(type = "numeric")
)

# Technical covariates included in the DESeq2 design but NOT tested.
# Must be column names in SAMPLE_METADATA_FILE. Set to character(0) if none.
COVARIATES <- character(0)           # e.g. c("batch", "rin")

# ── Sample metadata ───────────────────────────────────────────────────────────

# Path to a tab-delimited file mapping sample IDs to metadata variables.
# Required columns: sample_id, all names in BIO_VARS, and all in COVARIATES.
SAMPLE_METADATA_FILE <- "sample_metadata.tsv"

# ── SVA ───────────────────────────────────────────────────────────────────────

# Whether to run SVA to detect and account for hidden technical variation.
# When TRUE, surrogate variables are added to the DESeq2 design formula.
# SVA is never used to modify counts — only as additional covariates.
USE_SVA <- FALSE

# SVA estimation method: "be" (formal test, slower) or "leek" (faster)
SVA_METHOD <- "be"

# ── FDR threshold ─────────────────────────────────────────────────────────────

FDR_THRESHOLD <- 0.05

# ── Cell type ordering for strip plot ─────────────────────────────────────────
# Set to NULL to use alphabetical order
CELL_TYPE_ORDER <- NULL
# e.g. c("Granulosa", "Theca", "Stroma", "Myeloid", "B", "T")

# Strip plot y-axis limits.
# NULL = auto-detect from data (recommended — no points will be cut off).
# Set to a fixed range if you want consistent axes across multiple runs,
# e.g. for comparison between datasets.
STRIP_YLIM <- NULL               # e.g. c(-2, 2)

# ── Colors ────────────────────────────────────────────────────────────────────

# Named vector mapping sample group values to colors for MDS/PCA plots.
# Names should match the values of the first entry in BIO_VARS.
# Set to NULL for automatic color assignment.
SAMPLE_COLORS <- NULL
# e.g. c("4" = "deeppink1", "20" = "deeppink4")
