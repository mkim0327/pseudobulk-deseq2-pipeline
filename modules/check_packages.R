################################################################################
# modules/check_packages.R — Package installation check
################################################################################

check_and_install_packages <- function() {

  message("\n── Checking required packages ───────────────────────────────────────────────")

  packages <- list(
    # CRAN
    Seurat        = list(source = "CRAN",         repo = "Seurat"),
    ggplot2       = list(source = "CRAN",         repo = "ggplot2"),
    dplyr         = list(source = "CRAN",         repo = "dplyr"),
    tibble        = list(source = "CRAN",         repo = "tibble"),
    RColorBrewer  = list(source = "CRAN",         repo = "RColorBrewer"),
    pheatmap      = list(source = "CRAN",         repo = "pheatmap"),
    BiocManager   = list(source = "CRAN",         repo = "BiocManager"),
    # Bioconductor
    muscat            = list(source = "Bioconductor", repo = "muscat"),
    DESeq2            = list(source = "Bioconductor", repo = "DESeq2"),
    SingleCellExperiment = list(source = "Bioconductor", repo = "SingleCellExperiment"),
    sva               = list(source = "Bioconductor", repo = "sva"),
    limma             = list(source = "Bioconductor", repo = "limma"),
    ComplexHeatmap    = list(source = "Bioconductor", repo = "ComplexHeatmap"),
    DEGreport         = list(source = "Bioconductor", repo = "DEGreport")
  )

  pkg_names      <- names(packages)
  installed      <- rownames(installed.packages())
  missing        <- pkg_names[!pkg_names %in% installed]

  if (length(missing) == 0) {
    message("All packages already installed.")
    message("Loading packages...")
    for (pkg in pkg_names)
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    message("All packages loaded successfully.\n")
    return(invisible(TRUE))
  }

  message("Missing packages: ", paste(missing, collapse = ", "))

  missing_cran <- missing[vapply(missing, function(p) packages[[p]]$source == "CRAN",         logical(1))]
  missing_bioc <- missing[vapply(missing, function(p) packages[[p]]$source == "Bioconductor", logical(1))]

  if (length(missing_cran) > 0) {
    message("Installing from CRAN: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran, repos = "https://cloud.r-project.org", quiet = TRUE)
  }

  if (length(missing_bioc) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
    message("Installing from Bioconductor: ", paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }

  still_missing <- pkg_names[!pkg_names %in% rownames(installed.packages())]
  if (length(still_missing) > 0)
    stop("Could not install: ", paste(still_missing, collapse = ", "),
         "\nPlease install manually and re-run.")

  message("Loading packages...")
  for (pkg in pkg_names)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  message("All packages loaded successfully.\n")
  invisible(TRUE)
}
