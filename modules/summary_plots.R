################################################################################
# modules/summary_plots.R — Summary strip plot and result tables
#
# Generates one strip plot per biological variable and saves per-cell-type
# result tables.
#
# Expects in environment:
#   deseq_results — nested list: [[cell_type]][[bio_var]] = DESeqResults
#   BIO_VARS, FDR_THRESHOLD, CELL_TYPE_ORDER,
#   OUTPUT_DIR, PROJECT_NAME,
#   stamp(), stamp_pdf()
################################################################################

run_summary_plots <- function(deseq_results) {

  message("\n── Step 3: Summary plots ────────────────────────────────────────────────────")

  bio_var_names <- names(BIO_VARS)

  # ── Order cell types ──────────────────────────────────────────────────────────

  available <- if (!is.null(CELL_TYPE_ORDER))
                 c(intersect(CELL_TYPE_ORDER, names(deseq_results)),
                   setdiff(names(deseq_results), CELL_TYPE_ORDER))
               else
                 sort(names(deseq_results))

  res_ordered <- deseq_results[available]

  # ── Save per-cell-type, per-variable result tables ────────────────────────────

  for (ct in names(res_ordered)) {
    for (var in bio_var_names) {
      res <- res_ordered[[ct]][[var]]
      if (is.null(res) || nrow(res) == 0) next
      write.table(
        as.data.frame(res),
        file  = stamp_dir("DESeq2_results", ct, var, "DESeq2_results", ext = "txt"),
        sep   = "\t", quote = FALSE, row.names = TRUE
      )
    }
  }
  message("  Per-cell-type result tables saved.")

  # ── One strip plot per biological variable ─────────────────────────────────────

  for (var in bio_var_names) {

    # Collect results for this variable across all cell types
    var_results <- lapply(res_ordered, function(ct_res) ct_res[[var]])
    var_results <- var_results[!sapply(var_results, is.null)]

    if (length(var_results) == 0) next

    # Significant gene counts
    sig_counts <- sapply(var_results, function(r)
      sum(!is.na(r$padj) & r$padj < FDR_THRESHOLD))

    # Save summary table
    summary_df <- data.frame(
      cell_type   = names(sig_counts),
      bio_var     = var,
      n_sig_genes = as.integer(sig_counts),
      stringsAsFactors = FALSE
    )
    write.table(summary_df,
                file = stamp_dir("DESeq2_results", var, "significant_gene_summary", ext = "txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # Strip plot
    xlab_vec <- paste0(names(sig_counts), "\n(", sig_counts, " sig.)")
    all_lfc  <- unlist(lapply(var_results, function(r) r$log2FoldChange))
    # y-axis: use STRIP_YLIM from config if set, otherwise derive from data
    # with 5% padding. No clamping — all points are visible.
    if (!exists("STRIP_YLIM") || is.null(STRIP_YLIM)) {
      y_range <- range(all_lfc, na.rm = TRUE)
      y_pad   <- diff(y_range) * 0.05
      y_lim   <- c(y_range[1] - y_pad, y_range[2] + y_pad)
    } else {
      y_lim   <- STRIP_YLIM
    }
    y_label  <- if (BIO_VARS[[var]]$type == "numeric")
                  paste0("Log2FC per unit of ", var)
                else
                  paste0("Log2FC (", var, ")")

    set.seed(123123123)

    pdf(stamp_pdf("Summary", paste0("stripplot_", var, "_significant_genes")),
        width = 8, height = 5)
    par(mar = c(3.1, 4.1, 1.5, 1), oma = c(7, 2, 1, 1))

    plot(x = 1, y = 1, type = "n",
         xlim = c(0.5, length(var_results) + 0.5),
         ylim = y_lim, axes = FALSE,
         xlab = "", ylab = y_label,
         main = paste0("DEG summary — ", var))

    abline(h = 0)
    abline(h = pretty(y_lim), lty = "dotted", col = "gray80")

    for (i in seq_along(var_results)) {
      r      <- var_results[[i]]
      is_sig <- !is.na(r$padj) & r$padj < FDR_THRESHOLD
      cols   <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), nrow(r))
      cols[is_sig & r$log2FoldChange > 0] <- "deeppink4"
      cols[is_sig & r$log2FoldChange < 0] <- "deeppink1"
      points(x = jitter(rep(i, nrow(r)), amount = 0.2),
             y = r$log2FoldChange,
             pch = 16, col = cols, cex = 0.5)
    }

    axis(1, at = seq_along(var_results), tick = FALSE, las = 2, lwd = 0,
         labels = xlab_vec, cex.axis = 0.7)
    axis(2, las = 1, at = pretty(y_lim))
    box()
    dev.off()
    message("  Strip plot saved: stripplot_", var)
  }

  message("Summary complete.\n")
  invisible(NULL)
}

