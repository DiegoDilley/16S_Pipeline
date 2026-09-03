# ---------------------------------------------------------------------------
# 06 — Source tracking with FEAST
#
# Input : outputs/ps_filtered.rds, data/metadata.tsv
# Output: outputs/feast_results.csv, outputs/plot_sources.rds
#
# FEAST estimates what fraction of each sink community came from each
# candidate source, plus an "Unknown" fraction for everything unexplained.
# The Unknown fraction is a result, not an error term — in many designs it is
# the most interesting number in the table.
#
# The method is only as good as the source set. Anything real that you didn't
# sample lands in Unknown.
# ---------------------------------------------------------------------------

library(phyloseq)
library(FEAST)
library(tidyverse)

set.seed(42)
ps <- readRDS("outputs/ps_filtered.rds")

# ---------------------------------------------------------------------------
# THE WORKING-DIRECTORY GOTCHA
#
# FEAST() calls setwd() internally and does not restore it. Every relative
# path after the first call resolves somewhere unexpected, and in a loop the
# damage compounds. Capture the project root once and restore it after every
# call.
# ---------------------------------------------------------------------------
PROJECT_ROOT <- getwd()

# --- 1. Count matrix -------------------------------------------------------
# FEAST wants SAMPLES AS ROWS — the transpose of phyloseq's default layout.
# It also wants raw integer counts, not proportions.
counts <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) counts <- t(counts)
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

# --- 2. Metadata -----------------------------------------------------------
# Three columns, named exactly:
#   Env         label for the source or sink type
#   SourceSink  literally "Source" or "Sink"
#   id          integer linking each sink to its source set; NA for sources
#               when all sinks share one source pool
#
# Rownames must match the count matrix rownames, in the same order.
feast_meta <- read.table("data/metadata.tsv", sep = "\t",
                         header = TRUE, row.names = 1) %>%
  mutate(
    Env        = source_type,                       # e.g. "soil", "water", "sink_community"
    SourceSink = if_else(source_type == "sink_community", "Sink", "Source"),
    id         = if_else(SourceSink == "Sink", row_number(), NA_integer_)
  ) %>%
  select(Env, SourceSink, id)

# Align. A mismatch here produces a cryptic subscript error deep inside FEAST.
common <- intersect(rownames(counts), rownames(feast_meta))
counts     <- counts[common, , drop = FALSE]
feast_meta <- feast_meta[common, , drop = FALSE]

stopifnot(identical(rownames(counts), rownames(feast_meta)))
stopifnot(all(feast_meta$SourceSink %in% c("Source", "Sink")))

# --- 3. Run ----------------------------------------------------------------
# different_sources_flag:
#   0  all sinks share one source pool (the common case)
#   1  each sink has its own source set, matched via `id`
#
# EM_iterations: 1000 is the default and usually adequate. Raise it if
# repeated runs with different seeds give unstable proportions.
dir.create("outputs/feast", showWarnings = FALSE, recursive = TRUE)

feast_out <- FEAST(
  C = counts,
  metadata = feast_meta,
  different_sources_flag = 0,
  dir_path = file.path(PROJECT_ROOT, "outputs/feast"),
  outfile = "feast_run",
  EM_iterations = 1000
)

setwd(PROJECT_ROOT)   # <- restore. Do this after EVERY FEAST() call.

# --- 4. Tidy ---------------------------------------------------------------
# Output is a matrix: sinks as rows, sources plus Unknown as columns.
results <- as.data.frame(feast_out) %>%
  rownames_to_column("sink") %>%
  pivot_longer(-sink, names_to = "source", values_to = "proportion") %>%
  # Column names come back as "source_Env"; strip the suffix.
  mutate(source = str_remove(source, "_.*$")) %>%
  filter(!is.na(proportion))

# Sanity check: proportions per sink should sum to ~1.
results %>%
  group_by(sink) %>%
  summarise(total = sum(proportion)) %>%
  filter(abs(total - 1) > 0.01)

write_csv(results, "outputs/feast_results.csv")

# --- 5. Summarise ----------------------------------------------------------
results %>%
  group_by(source) %>%
  summarise(mean_contribution = mean(proportion),
            sd = sd(proportion), .groups = "drop") %>%
  arrange(desc(mean_contribution))

# --- 6. Plot ---------------------------------------------------------------
# Unknown pinned last so it reads as the residual it is.
source_levels <- setdiff(unique(results$source), "Unknown")

p_sources <- results %>%
  mutate(source = factor(source, levels = c(source_levels, "Unknown"))) %>%
  ggplot(aes(x = sink, y = proportion, fill = source)) +
  geom_col(width = 0.9) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Estimated source contribution", fill = "Source")

saveRDS(p_sources, "outputs/plot_sources.rds")

# --- 7. Following up on Unknown --------------------------------------------
# A large Unknown fraction is worth characterising rather than reporting as a
# bare number. Script 05's subset workflow applies directly: pull the features
# that are abundant in sinks and near-absent from every source, and describe
# them. That set is either an unsampled source or a genuinely resident
# community, and the taxonomy usually tells you which.
