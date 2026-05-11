# Pseudobulk DESeq2 Pipeline

An R pipeline for pseudobulk differential expression analysis from annotated Seurat objects, using **muscat** for aggregation and **DESeq2** for differential expression.

> **Note:** This pipeline was developed with assistance from [Claude](https://claude.ai) (Anthropic).

---

## Overview

```
Seurat object → muscat aggregation → DESeq2 (+ optional SVA) → plots + tables
```

Supports:
- **Numeric or categorical** primary variables (e.g. age in months, or treatment vs control)
- **Optional SVA** for detecting and accounting for hidden technical variation
- **Batch and covariate correction** via DESeq2 design formula
- **MDS and PCA plots** per cell type
- **Heatmaps** of significant genes per cell type
- **Strip plot** summarizing log2FC across all cell types

---

## Repository Structure

```
├── pipeline.R                  # Orchestrator — runs all steps
├── config.R                    # Template config (edit per project)
├── sample_metadata.tsv         # Template sample sheet
├── config_AC.R                 # Example: Benayoun lab AC dataset
├── config_VCD.R                # Example: Benayoun lab VCD dataset
├── sample_metadata_AC.tsv      # AC sample metadata
├── sample_metadata_VCD.tsv     # VCD sample metadata
└── modules/
    ├── check_packages.R        # Package installation and loading
    ├── utils.R                 # Shared helpers (stamp, colors, design)
    ├── pseudobulk.R            # Step 1: muscat aggregation
    ├── deseq2_analysis.R       # Step 2: DESeq2 per cell type
    └── summary_plots.R         # Step 3: strip plot + result tables
```

---

## Prerequisites

All packages are auto-installed by `check_packages.R` on first run.

| Package | Source |
|---------|--------|
| `Seurat` | CRAN |
| `ggplot2`, `dplyr`, `pheatmap`, `RColorBrewer` | CRAN |
| `muscat`, `DESeq2`, `sva`, `limma` | Bioconductor |
| `SingleCellExperiment`, `ComplexHeatmap`, `DEGreport` | Bioconductor |

---

## Quick Start

### 1. Edit config.R and sample_metadata.tsv

```r
# Minimum required settings in config.R:
SEURAT_PATH      <- "/path/to/your/seurat_object.rds"  # preferred format
PRIMARY_VAR      <- "age"                  # column in sample_metadata.tsv
PRIMARY_VAR_TYPE <- "numeric"              # "numeric" or "categorical"
COVARIATES       <- c("batch")             # technical covariates (or character(0))
USE_SVA          <- FALSE                  # TRUE to enable SVA
```

### 2. Run

```bash
Rscript pipeline.R config.R
```

---

## Configuration Reference

### Primary variable type

| Setting | Use case | log2FC interpretation |
|---------|----------|----------------------|
| `"numeric"` | Continuous variable (age in months, score) | Per unit increase |
| `"categorical"` | Discrete groups (treatment vs control) | Group vs reference |

For categorical, set `REFERENCE_LEVEL` to the baseline group (e.g. `"CTL"`).

### SVA usage

When `USE_SVA = TRUE`:
- SVA runs on log2 CPM-normalized counts to detect surrogate variables
- Detected SVs are added as covariates in the **DESeq2 design formula**
- SVA-corrected values (via `limma::removeBatchEffect`) are computed **only for visualization** (MDS, PCA, heatmaps)

When `USE_SVA = FALSE`:
- Known technical variables (`COVARIATES`) are included directly in the DESeq2 design
- This is the recommended starting point — add SVA only if unexplained variation is visible in PCA

### Output directory layout

```
OUTPUT_DIR/
├── DESeq2_results/
│   ├── YYYY-MM-DD_<PROJECT>_<CellType>_DESeq2_results.txt
│   ├── YYYY-MM-DD_<PROJECT>_PB_VST_counts.RData
│   └── YYYY-MM-DD_<PROJECT>_PB_DESeq2_results.RData
├── Heatmaps/
│   └── YYYY-MM-DD_<PROJECT>_<CellType>_heatmap_FDR5.pdf
├── MDS_Plots/
│   └── YYYY-MM-DD_<PROJECT>_<CellType>_MDS_plot.pdf
├── PCA_Plots/
│   └── YYYY-MM-DD_<PROJECT>_<CellType>_PCA_plot.pdf
├── Summary/
│   ├── YYYY-MM-DD_<PROJECT>_stripplot_significant_genes.pdf
│   └── YYYY-MM-DD_<PROJECT>_significant_gene_summary.txt
├── YYYY-MM-DD_<PROJECT>_muscat_pseudobulk_MDS.pdf
├── YYYY-MM-DD_<PROJECT>_cell_counts_per_sample.txt
└── YYYY-MM-DD_<PROJECT>_session_info.txt
```

---

## Notes

- **Interaction terms** are not currently supported. All biological variables in `BIO_VARS` are included as main effects only (e.g. `~ batch + age + sex`). If you need to test whether the effect of one variable depends on another (e.g. does the age effect differ by sex), this would require adding interaction terms to the design formula and is planned for a future version.
- **Cell type filtering**: a cell type is included only if it has `>= MIN_CELLS_PER_SAMPLE` cells in `>= MIN_SAMPLES_PASSING` samples. Adjust these thresholds in `config.R` based on your dataset size.
- **Covariate filtering**: covariates with fewer than 2 unique values within a cell type's sample set are automatically dropped to avoid `model.matrix` errors (e.g. batch is always constant within a single sample, so it is never included in per-sample SCT steps).

---

## Acknowledgements

Pipeline code generated with assistance from [Claude](https://claude.ai) (claude-sonnet-4-6, Anthropic).

### References

- **muscat**: Crowell et al. (2020). *muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data.* Nature Communications. https://doi.org/10.1038/s41467-020-19894-4
- **DESeq2**: Love et al. (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology. https://doi.org/10.1186/s13059-014-0550-8
- **SVA**: Leek et al. (2012). *The sva package for removing batch effects and other unwanted variation in high-throughput experiments.* Bioinformatics. https://doi.org/10.1093/bioinformatics/bts034
