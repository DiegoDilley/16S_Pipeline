# ---------------------------------------------------------------------------
# 00 — Dependencies
#
# Run once. Grouped by where the package actually comes from, because that is
# the thing that breaks for someone cloning a repo.
# ---------------------------------------------------------------------------

# --- CRAN ------------------------------------------------------------------
cran_pkgs <- c(
  "tidyverse",    # dplyr, ggplot2, tidyr, readr, tibble
  "vegan",        # distances, NMDS, adonis2, betadisper
  "remotes",      # GitHub installs
  "BiocManager",  # Bioconductor installs
  "patchwork",    # multi-panel figures
  "RColorBrewer", # palettes
  "ggrepel"       # non-overlapping labels on ordinations
)

missing <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing)

# --- Bioconductor ----------------------------------------------------------
# phyloseq is not on CRAN. This is the most common first stumble.
if (!requireNamespace("phyloseq", quietly = TRUE)) {
  BiocManager::install("phyloseq")
}

# --- GitHub ----------------------------------------------------------------
# None of these three are on CRAN. Install lines are given explicitly rather
# than left as "install the required packages".
if (!requireNamespace("FEAST", quietly = TRUE)) {
  remotes::install_github("cozygene/FEAST")
}
if (!requireNamespace("microeco", quietly = TRUE)) {
  remotes::install_github("ChiLiubio/microeco")
}
if (!requireNamespace("pairwiseAdonis", quietly = TRUE)) {
  remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
}

# --- External, non-R -------------------------------------------------------
# FAPROTAX (script 07) is a Python script, not an R package. Download it from
# pages.uoregon.edu/slouca/LoucaLab and note the version you used.

# --- Record the environment ------------------------------------------------
# Capture this once and commit it. It answers "which version of vegan?" two
# years later, when the answer matters and nobody remembers.
writeLines(capture.output(sessionInfo()), "outputs/sessionInfo.txt")
