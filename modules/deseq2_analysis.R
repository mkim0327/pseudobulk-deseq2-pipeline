################################################################################
# modules/deseq2_analysis.R — Per-cell-type pseudobulk DESeq2
#
# For each cell type:
#   1. Filter lowly expressed genes
#   2. Cast biological variables to correct types (numeric/categorical)
#      Numeric variables are centered and scaled to improve DESeq2 convergence
#      and avoid collinearity with the intercept.
#   3. Optionally run SVA — SVs added as design covariates only, never modifying
#      counts. removeBatchEffect used for visualization only.
#   4. Run DESeq2 with all biological variables + technical covariates + SVs
#   5. Extract results for each biological variable separately
#   6. Generate MDS and PCA plots (coloured by first biological variable)
#   7. Generate heatmaps per biological variable
#
# Expects in environment:
#   counts_pb, sample_metadata,
#   BIO_VARS  — named list: list(age=list(type="numeric"), sex=list(type="categorical", reference="F"))
#   COVARIATES, USE_SVA, SVA_METHOD,
#   FDR_THRESHOLD, MIN_SAMPLES_EXPRESSED,
#   OUTPUT_DIR, PROJECT_NAME,
#   stamp(), stamp_pdf(), make_sample_colors(), build_design()
#
# Returns: list(vst_counts, deseq_results)
#   deseq_results — nested list: [[cell_type]][[bio_var]] = DESeqResults
################################################################################

run_deseq2_analysis <- function(counts_pb, sample_metadata) {

  message("\n── Step 2: DESeq2 analysis ──────────────────────────────────────────────────")

  bio_var_names <- names(BIO_VARS)

  # Results nested by cell type then by biological variable
  vst_counts    <- vector("list", length(counts_pb))
  deseq_results <- vector("list", length(counts_pb))
  names(vst_counts) <- names(deseq_results) <- names(counts_pb)

  for (cell_type in names(counts_pb)) {

    message("\n  ── Cell type: ", cell_type, " ──")

    # ── 1. Sample and gene filtering ──────────────────────────────────────────

    counts <- counts_pb[[cell_type]]

    # Drop samples with all-zero counts (muscat keeps zero columns for samples
    # that fall below MIN_CELLS_PER_SAMPLE; they break size factor estimation)
    sample_totals <- colSums(counts)
    if (any(sample_totals == 0)) {
      zero_samps <- names(sample_totals[sample_totals == 0])
      message("    Removing ", length(zero_samps),
              " all-zero sample(s): ", paste(zero_samps, collapse = ", "))
      counts <- counts[, sample_totals > 0, drop = FALSE]
    }
    if (ncol(counts) < MIN_SAMPLES_PASSING) {
      message("    Skipping — only ", ncol(counts), " non-zero sample(s) remain",
              " (< MIN_SAMPLES_PASSING = ", MIN_SAMPLES_PASSING, ").")
      next
    }

    good_genes <- rowSums(counts > 0) >= MIN_SAMPLES_EXPRESSED
    counts     <- counts[good_genes, ]

    # Coerce to integer (handles non-integer source assays, e.g. Tabula Sapiens)
    storage.mode(counts) <- "integer"

    if (nrow(counts) == 0) {
      message("    Skipping — no genes passed expression filter.")
      next
    }

    # ── 2. Build colData ───────────────────────────────────────────────────────

    missing_meta <- setdiff(colnames(counts), rownames(sample_metadata))
    if (length(missing_meta) > 0) {
      stop(
        "Sample ID mismatch — samples in count matrix not found in metadata.\n",
        "  In count matrix : ", paste(sort(colnames(counts)),         collapse = ", "), "\n",
        "  In metadata     : ", paste(sort(rownames(sample_metadata)), collapse = ", "), "\n",
        "  Missing         : ", paste(missing_meta, collapse = ", "), "\n",
        "Fix: update the sample_id column in your metadata file to exactly match\n",
        "the values in the '", SAMPLE_COL, "' column of your Seurat object."
      )
    }

    col_data <- sample_metadata[colnames(counts), , drop = FALSE]

    # Cast each biological variable to its declared type.
    # Numeric variables are centered and scaled for DESeq2.
    # Original (unscaled) values are stored separately for visualization.
    scale_factors  <- list()
    original_vals  <- list()   # unscaled values used for plot colors/ordering

    skip_ct <- FALSE
    for (bv in bio_var_names) {
      spec <- BIO_VARS[[bv]]
      if (spec$type == "numeric") {
        raw_vals           <- as.numeric(col_data[[bv]])
        original_vals[[bv]] <- raw_vals          # keep raw for plots
        bv_mean            <- mean(raw_vals, na.rm = TRUE)
        bv_sd              <- stats::sd(raw_vals, na.rm = TRUE)
        if (is.na(bv_sd) || bv_sd == 0) {
          message("    [", bv, "] Zero variance after sample filtering",
                  " (all remaining samples share the same value) — skipping cell type.")
          skip_ct <- TRUE
          break
        }
        col_data[[bv]]     <- (raw_vals - bv_mean) / bv_sd
        scale_factors[[bv]] <- list(mean = bv_mean, sd = bv_sd)
        message("    [", bv, "] Centered and scaled (mean=", round(bv_mean, 2),
                ", sd=", round(bv_sd, 2), ") — log2FC per SD of ", bv)
      } else {
        ref               <- if (!is.null(spec$reference)) spec$reference else
                                levels(factor(col_data[[bv]]))[1]
        col_data[[bv]]    <- relevel(factor(col_data[[bv]]), ref = ref)
        original_vals[[bv]] <- as.character(col_data[[bv]])   # factor levels for plots
      }
    }
    if (skip_ct) next

    # Validate and drop technical covariates with < 2 levels
    valid_covs <- COVARIATES[vapply(COVARIATES, function(v) {
      v %in% colnames(col_data) && length(unique(col_data[[v]])) >= 2
    }, logical(1))]

    dropped_covs <- setdiff(COVARIATES, valid_covs)
    if (length(dropped_covs) > 0)
      message("    Dropping covariate(s) with < 2 levels: ",
              paste(dropped_covs, collapse = ", "))

    # ── 3. Optional SVA ────────────────────────────────────────────────────────
    # mod1 includes ALL biological variables so SVA captures only unexplained
    # variation, not variation due to age, sex, treatment, etc.

    sv_names      <- character(0)
    vst_for_plots <- NULL

    if (USE_SVA) {
      message("    Running SVA (method = '", SVA_METHOD, "')...")

      lib_sizes <- colSums(counts)
      log_cpm   <- log2(sweep(counts, 2, lib_sizes, "/") * 1e6 + 1)

      mod1 <- model.matrix(build_design(bio_var_names, valid_covs), data = col_data)
      mod0 <- model.matrix(~ 1, data = col_data)

      n_sv <- tryCatch(
        sva::num.sv(as.matrix(log_cpm), mod1, method = SVA_METHOD),
        error = function(e) {
          message("    num.sv failed: ", conditionMessage(e), " — skipping SVA.")
          0L
        }
      )

      if (n_sv > 0) {
        message("    Detected ", n_sv, " surrogate variable(s).")
        svseq    <- sva::svaseq(as.matrix(counts), mod1, mod0, n.sv = n_sv)
        sv_df    <- as.data.frame(svseq$sv)
        sv_names <- paste0("SV", seq_len(n_sv))
        colnames(sv_df) <- sv_names
        col_data <- cbind(col_data, sv_df)

        # Corrected values for VISUALIZATION ONLY
        vis_covariates <- cbind(
          svseq$sv,
          if (length(valid_covs) > 0)
            as.matrix(col_data[, valid_covs, drop = FALSE])
          else NULL
        )
        vis_design    <- model.matrix(build_design(bio_var_names), data = col_data)
        vst_for_plots <- limma::removeBatchEffect(log2(counts + 0.5),
                                                   covariates = vis_covariates,
                                                   design     = vis_design)
      } else {
        message("    No surrogate variables detected — SVA not applied.")
      }
    }

    # ── 4. DESeq2 on raw counts ────────────────────────────────────────────────

    design_formula <- build_design(bio_var_names, valid_covs, sv_names)
    message("    DESeq2 design: ", deparse(design_formula))

    dds      <- DESeqDataSetFromMatrix(countData = counts,
                                        colData   = col_data,
                                        design    = design_formula)

    # Pre-estimate size factors; fall back to poscounts for sparse cell types
    # where every gene has at least one zero across samples.
    dds <- tryCatch(
      estimateSizeFactors(dds),
      error = function(e) {
        if (grepl("every gene contains at least one zero", conditionMessage(e))) {
          message("    Sparse data — retrying with poscounts size factor estimation.")
          estimateSizeFactors(dds, type = "poscounts")
        } else stop(e)
      }
    )
    dds      <- DESeq(dds, quiet = TRUE)
    vst_data <- getVarianceStabilizedData(dds)

    vst_counts[[cell_type]] <- vst_data

    # ── 5. Extract results for each biological variable ────────────────────────

    cell_res <- vector("list", length(bio_var_names))
    names(cell_res) <- bio_var_names

    for (bv in bio_var_names) {
      spec         <- BIO_VARS[[bv]]
      result_names <- resultsNames(dds)

      res <- if (spec$type == "numeric") {
        coef_name <- if (bv %in% result_names) {
          bv
        } else {
          matched <- result_names[grep(paste0("^", bv), result_names)]
          if (length(matched) == 0)
            stop("Cannot find DESeq2 coefficient for numeric variable '", bv,
                 "'. resultsNames: ", paste(result_names, collapse = ", "))
          matched[1]
        }
        results(dds, name = coef_name)

      } else {
        bv_levels   <- levels(col_data[[bv]])
        ref_level   <- bv_levels[1]
        test_levels <- bv_levels[-1]

        if (length(test_levels) == 1) {
          results(dds, contrast = c(bv, test_levels[1], ref_level))
        } else {
          message("    [", bv, "] Multiple non-reference levels: ",
                  paste(test_levels, collapse = ", "),
                  " — extracting '", test_levels[1], "' vs '", ref_level, "'")
          results(dds, contrast = c(bv, test_levels[1], ref_level))
        }
      }

      res             <- res[!is.na(res$padj), ]
      cell_res[[bv]]  <- res

      n_sig <- sum(res$padj < FDR_THRESHOLD, na.rm = TRUE)
      sf    <- scale_factors[[bv]]
      if (!is.null(sf))
        message("    [", bv, "] Significant genes (FDR < ", FDR_THRESHOLD, "): ", n_sig,
                "  (log2FC per SD = ", round(sf$sd, 2), " ", bv, ")")
      else
        message("    [", bv, "] Significant genes (FDR < ", FDR_THRESHOLD, "): ", n_sig)
    }

    deseq_results[[cell_type]] <- cell_res

    # ── 6. MDS and PCA — coloured by first biological variable ─────────────────
    # Use original (unscaled) values for colors and labels

    plot_data   <- if (!is.null(vst_for_plots)) vst_for_plots else vst_data
    color_bv    <- bio_var_names[1]
    sample_vals <- original_vals[[color_bv]]   # unscaled for display
    col_vec     <- make_sample_colors(sample_vals)
    outprefix   <- paste0(Sys.Date(), "_", PROJECT_NAME, "_PB_DESeq2_", cell_type)

    # MDS
    mds_result <- cmdscale(1 - cor(plot_data, method = "spearman"), k = 2)
    pdf(stamp_pdf("MDS_Plots", outprefix, "MDS_plot"))
    plot(mds_result[, 1], mds_result[, 2],
         xlab = "MDS dimension 1", ylab = "MDS dimension 2",
         main = paste0(cell_type, " MDS (coloured by ", color_bv, ")"),
         cex = 3, pch = 16, col = col_vec, cex.lab = 1.3, cex.axis = 1.2)
    legend("topright", legend = unique(as.character(sample_vals)),
           col = unique(col_vec), pch = 16, cex = 0.8, title = color_bv)
    dev.off()

    # PCA
    gene_vars <- apply(plot_data, 1, stats::var)
    pos_var   <- gene_vars > 0
    if (sum(pos_var) >= 2) {
      pca_res <- prcomp(t(plot_data[pos_var, ]), scale. = TRUE)
      pca_sum <- summary(pca_res)
      pc1_pct <- round(100 * pca_sum$importance[2, 1], 1)
      pc2_pct <- round(100 * pca_sum$importance[2, 2], 1)
      x_range <- range(pca_res$x[, 1]); x_pad <- diff(x_range) * 0.15
      y_range <- range(pca_res$x[, 2]); y_pad <- diff(y_range) * 0.15

      pdf(stamp_pdf("PCA_Plots", outprefix, "PCA_plot"))
      plot(pca_res$x[, 1], pca_res$x[, 2],
           xlab  = paste0("PC1 (", pc1_pct, "%)"),
           ylab  = paste0("PC2 (", pc2_pct, "%)"),
           main  = paste0(cell_type, " PCA (coloured by ", color_bv, ")"),
           xlim  = c(x_range[1] - x_pad, x_range[2] + x_pad),
           ylim  = c(y_range[1] - y_pad, y_range[2] + y_pad),
           cex   = 3, pch = 16, col = col_vec, cex.lab = 1.3, cex.axis = 1.2)
      text(pca_res$x[, 1], pca_res$x[, 2],
           labels = colnames(plot_data), cex = 0.6, pos = 3)
      legend("topright", legend = unique(as.character(sample_vals)),
             col = unique(col_vec), pch = 16, cex = 0.8, title = color_bv)
      dev.off()
    }

    # ── 7. Heatmap per biological variable ────────────────────────────────────
    # Order columns by original (unscaled) values

    for (bv in bio_var_names) {
      res       <- cell_res[[bv]]
      sig_genes <- rownames(res)[res$padj < FDR_THRESHOLD]

      if (length(sig_genes) > 2) {
        ordered_cols <- colnames(plot_data)[order(original_vals[[bv]])]
        heatmap_data <- plot_data[sig_genes, ordered_cols]

        pdf(stamp_pdf("Heatmaps", paste0(cell_type, "_", bv), "heatmap_FDR5"))
        pheatmap::pheatmap(
          heatmap_data,
          cluster_cols  = FALSE,
          cluster_rows  = TRUE,
          show_rownames = FALSE,
          scale         = "row",
          color         = colorRampPalette(rev(c(
            "#CC3333", "#FF9999", "#FFCCCC", "white",
            "#CCCCFF", "#9999FF", "#333399"
          )))(50),
          main         = paste0(cell_type, " | ", bv,
                                " — FDR < ", FDR_THRESHOLD,
                                " (", length(sig_genes), " genes)"),
          cellwidth    = 20,
          cellheight   = 1,
          border_color = NA
        )
        dev.off()
      }
    }
  }

  invisible(list(vst_counts = vst_counts, deseq_results = deseq_results))
}
