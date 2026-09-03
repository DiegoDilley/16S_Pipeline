# ---------------------------------------------------------------------------
# 01 — From exported tables to a filtered phyloseq object
#
# Input : data/feature_table.tsv, data/taxonomy.tsv, data/metadata.tsv
# Output: outputs/ps_raw.rds, outputs/ps_filtered.rds
#
# This is the handoff point from QIIME2. Everything upstream (import,
# demultiplexing, denoising, classification) happened outside R.
# ---------------------------------------------------------------------------

library(phyloseq)
library(tidyverse)

set.seed(42)

# --- 1. Feature table ------------------------------------------------------
# A QIIME2 biom export starts with a "# Constructed from biom file" comment
# line. skip = 1 drops it. check.names = FALSE stops R prefixing sample IDs
# that begin with a digit with an "X" — a silent way to break the join to
# metadata later.
counts <- read.table(
  "data/feature_table.tsv",
  sep = "\t", header = TRUE, row.names = 1,
  skip = 1, comment.char = "", check.names = FALSE
)
counts <- as.matrix(counts)

# --- 2. Taxonomy -----------------------------------------------------------
# One column of semicolon-delimited lineage, split into ranks. fill = "right"
# handles features only classified to family or above.
taxonomy <- read.table(
  "data/taxonomy.tsv",
  sep = "\t", header = TRUE, row.names = 1,
  quote = "", comment.char = ""
) %>%
  separate(
    Taxon,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = ";\\s*", fill = "right", remove = TRUE
  ) %>%
  # Strip the "d__" / "p__" rank prefixes SILVA and GreenGenes add.
  mutate(across(Kingdom:Species, ~ str_remove(.x, "^[a-z]__"))) %>%
  # Empty strings are not the same as NA to subset_taxa(). Make them NA now.
  mutate(across(Kingdom:Species, ~ na_if(str_trim(.x), "")))

taxonomy <- as.matrix(taxonomy[, c("Kingdom", "Phylum", "Class",
                                   "Order", "Family", "Genus", "Species")])

# --- 3. Metadata -----------------------------------------------------------
metadata <- read.table(
  "data/metadata.tsv",
  sep = "\t", header = TRUE, row.names = 1, comment.char = ""
)

# Set factor levels explicitly. Otherwise R orders them alphabetically and
# every plot and model contrast silently inherits that order.
metadata$group <- factor(metadata$group, levels = c("group_a", "group_b"))

# --- 4. Assemble -----------------------------------------------------------
ps_raw <- phyloseq(
  otu_table(counts, taxa_are_rows = TRUE),
  tax_table(taxonomy),
  sample_data(metadata)
)

# Check this. phyloseq intersects the three tables silently, so a sample-ID
# mismatch shows up as samples quietly vanishing rather than as an error.
message("Samples in metadata: ", nrow(metadata))
message("Samples in object:   ", nsamples(ps_raw))
message("Features:            ", ntaxa(ps_raw))

saveRDS(ps_raw, "outputs/ps_raw.rds")

# --- 5. Filter -------------------------------------------------------------
# Order matters: drop non-target lineages first, then filter on prevalence,
# so prevalence is computed on the features you're actually keeping.

# 5a. Non-target lineages. The is.na() guard is essential — comparing NA to a
# string yields NA, and subset_taxa drops NA rows, so `Order != "Chloroplast"`
# alone would also throw away everything unclassified at order level.
ps <- subset_taxa(
  ps_raw,
  !is.na(Kingdom) & Kingdom != "Unassigned" &
    (is.na(Order)  | Order  != "Chloroplast") &
    (is.na(Family) | Family != "Mitochondria")
)

# 5b. Prevalence. Removes features seen in fewer than `min_samples` samples
# with at least `min_count` reads. Singletons dominate feature counts and add
# noise to distances without adding signal.
min_count   <- 2
min_samples <- 3
keep <- genefilter_sample(ps, filterfun_sample(function(x) x >= min_count),
                          A = min_samples)
ps <- prune_taxa(keep, ps)

# 5c. Drop samples with too few reads to be informative. Inspect the
# distribution before picking a threshold rather than reaching for a round one.
sort(sample_sums(ps))
ps <- prune_samples(sample_sums(ps) >= 1000, ps)

# 5d. Features left at zero after dropping samples.
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

message("After filtering: ", ntaxa(ps), " features, ", nsamples(ps), " samples")
saveRDS(ps, "outputs/ps_filtered.rds")

# --- 6. Rarefaction --------------------------------------------------------
# Contested. Rarefying discards data and is a poor fit for differential
# abundance, but it remains the common choice for richness estimates, which
# are sensitive to sequencing depth in a way Shannon largely isn't.
#
# The pragmatic split used in these scripts:
#   - alpha diversity  -> rarefied counts (script 02)
#   - beta diversity   -> relative abundance, not rarefied (script 03)
#   - composition      -> relative abundance (script 04)
#
# Whatever you choose, state it in the methods and keep it consistent.
depth <- min(sample_sums(ps))
message("Rarefying to ", depth)

ps_rare <- rarefy_even_depth(
  ps,
  sample.size = depth,
  rngseed = 42,          # required for reproducibility; without it, each run differs
  replace = FALSE,
  verbose = FALSE
)

saveRDS(ps_rare, "outputs/ps_rarefied.rds")
