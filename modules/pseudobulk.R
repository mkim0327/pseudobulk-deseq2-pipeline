################################################################################
# modules/pseudobulk.R — Pseudobulk aggregation with muscat
#
# Converts Seurat object to SCE, aggregates counts per cell type × sample,
# generates a global MDS plot, and filters cell types by minimum cell counts.
#
# Expects in environment:
#   SEURAT_PATH, CELLTYPE_COL, SAMPLE_COL, GROUP_COL,
#   MIN_CELLS_PER_SAMPLE, MIN_SAMPLES_PASSING,
#   PROJECT_NAME, OUTPUT_DIR, stamp(), stamp_pdf()
#
# Returns: list(pb, counts_pb, sce)
################################################################################

run_pseudobulk <- function() {

  message("\n── Step 1: Load and aggregate pseudobulk counts ─────────────────────────────")

  # ── Load Seurat object ────────────────────────────────────────────────────────

  if (!file.exists(SEURAT_PATH))
    stop("Seurat object not found: ", SEURAT_PATH)

  ext <- tolower(tools::file_ext(SEURAT_PATH))
  if (!ext %in% c("rds", "rdata", "rda"))
    stop("SEURAT_PATH must be a .rds, .RData, or .rda file. Got: .", ext)

  message("  Loading: ", SEURAT_PATH)

  if (ext == "rds") {
    seurat_obj <- readRDS(SEURAT_PATH)
  } else {
    # RData fallback for backwards compatibility
    env <- new.env()
    load(SEURAT_PATH, envir = env)
    obj_names <- ls(env)
    seurat_objs <- obj_names[sapply(obj_names, function(n)
      inherits(get(n, envir = env), "Seurat"))]
    if (length(seurat_objs) == 0)
      stop("No Seurat object found in: ", SEURAT_PATH)
    if (length(seurat_objs) > 1)
      stop("Multiple Seurat objects found in: ", SEURAT_PATH,
           "\nPlease convert to .rds format: saveRDS(your_obj, 'object.rds')")
    seurat_obj <- get(seurat_objs[1], envir = env)
    rm(env)
  }
  gc()

  if (!inherits(seurat_obj, "Seurat"))
    stop("Loaded object is not a Seurat object. Check SEURAT_PATH.")

  DefaultAssay(seurat_obj) <- "RNA"

  # ── Validate metadata columns before conversion ───────────────────────────────

  meta_cols    <- colnames(seurat_obj@meta.data)
  required_sce <- c(CELLTYPE_COL, SAMPLE_COL, GROUP_COL)
  missing_sce  <- setdiff(required_sce, meta_cols)

  if (length(missing_sce) > 0)
    stop(
      "Required column(s) not found in Seurat metadata:\n",
      "  Missing  : ", paste(missing_sce, collapse = ", "), "\n",
      "  Available: ", paste(sort(meta_cols), collapse = ", "), "\n",
      "Fix: update CELLTYPE_COL, SAMPLE_COL, or GROUP_COL in config.R."
    )

  # ── Convert to SCE and aggregate ─────────────────────────────────────────────

  message("  Converting to SingleCellExperiment...")
  sce <- as.SingleCellExperiment(seurat_obj)
  rm(seurat_obj); gc()

  sce <- prepSCE(sce,
                 kid  = CELLTYPE_COL,
                 gid  = GROUP_COL,
                 sid  = SAMPLE_COL,
                 drop = TRUE)

  nk <- length(kids <- levels(sce$cluster_id))
  ns <- length(sids <- levels(sce$sample_id))
  names(kids) <- kids
  names(sids) <- sids
  message("  Cell types: ", nk, "  |  Samples: ", ns)
  message("  Sample IDs from Seurat '", SAMPLE_COL, "' column: ",
          paste(sort(sids), collapse = ", "))

  pb <- aggregateData(sce, assay = "counts", fun = "sum",
                      by = c("cluster_id", "sample_id"))

  # ── Global MDS plot ───────────────────────────────────────────────────────────

  pb_mds <- pbMDS(pb)
  pdf(stamp_pdf("", "muscat_pseudobulk_MDS"))
  print(pb_mds)
  dev.off()
  message("  Global MDS plot saved.")

  # ── Filter cell types by minimum cell count ───────────────────────────────────

  cell_per_samp    <- t(table(sce$cluster_id, sce$sample_id))
  celltype_pass    <- colnames(cell_per_samp)[
    apply(cell_per_samp >= MIN_CELLS_PER_SAMPLE, 2, sum) >= MIN_SAMPLES_PASSING
  ]

  message("  Cell types passing QC (>= ", MIN_CELLS_PER_SAMPLE, " cells in >= ",
          MIN_SAMPLES_PASSING, " samples): ", length(celltype_pass))
  message("  ", paste(celltype_pass, collapse = ", "))

  if (length(celltype_pass) == 0)
    stop("No cell types passed QC thresholds. ",
         "Consider lowering MIN_CELLS_PER_SAMPLE or MIN_SAMPLES_PASSING.")

  counts_pb <- pb@assays@data[celltype_pass]

  # Save cell count table
  write.table(cell_per_samp,
              file = stamp("cell_counts_per_sample", ext = "txt"),
              sep = "\t", quote = FALSE)

  invisible(list(pb = pb, counts_pb = counts_pb, sce = sce))
}
