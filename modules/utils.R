################################################################################
# modules/utils.R — Shared helpers
################################################################################

#' Build a dated, project-stamped output path
stamp <- function(..., ext) {
  fname <- paste(c(as.character(Sys.Date()), PROJECT_NAME, ...), collapse = "_")
  file.path(OUTPUT_DIR, paste0(fname, ".", ext))
}

#' Build a dated, project-stamped path in a named subdirectory
stamp_dir <- function(subdir, ..., ext) {
  fname <- paste(c(as.character(Sys.Date()), PROJECT_NAME, ...), collapse = "_")
  file.path(OUTPUT_DIR, subdir, paste0(fname, ".", ext))
}

#' Build a dated, project-stamped path for PDF files in a subdirectory
stamp_pdf <- function(subdir, ...) {
  fname <- paste(c(as.character(Sys.Date()), PROJECT_NAME, ...), collapse = "_")
  file.path(OUTPUT_DIR, subdir, paste0(fname, ".pdf"))
}

#' Assign colors to samples based on a grouping variable's values.
#' Uses SAMPLE_COLORS from config if provided, otherwise assigns automatically.
make_sample_colors <- function(values) {
  if (!is.null(SAMPLE_COLORS)) {
    cols <- SAMPLE_COLORS[as.character(values)]
    cols[is.na(cols)] <- "grey60"
    return(unname(cols))
  }
  # Auto-assign from a palette
  unique_vals <- unique(values)
  palette     <- if (length(unique_vals) <= 8)
                   RColorBrewer::brewer.pal(max(3, length(unique_vals)), "Set1")[seq_along(unique_vals)]
                 else
                   colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(unique_vals))
  color_map   <- setNames(palette, unique_vals)
  unname(color_map[as.character(values)])
}

#' Build DESeq2 design formula from biological variables, technical covariates,
#' and optional surrogate variable names.
#' Terms are ordered: covariates first, then SVs, then biological variables
#' (last term is the one DESeq2 uses for the default results() call).
build_design <- function(bio_var_names, covariates = character(0),
                          sv_names = character(0)) {
  terms <- c(covariates, sv_names, bio_var_names)
  as.formula(paste("~", paste(terms, collapse = " + ")))
}
