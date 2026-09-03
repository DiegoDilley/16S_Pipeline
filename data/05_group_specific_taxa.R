# ---------------------------------------------------------------------------
# 05 — Group-specific taxa, and describing a subset
#
# Input : outputs/ps_filtered.rds
# Output: outputs/unique_taxa.csv, outputs/subset_summary.csv
#
# Two related jobs:
#   A. which features occur in one group and not another
#   B. take an arbitrary slice of the community and characterise it
# ---------------------------------------------------------------------------

library(phyloseq)
library(tidyverse)

ps <- readRDS("outputs/ps_filtered.rds")

# ===========================================================================
# A. Taxa unique to a group
# ===========================================================================

# --- 1. Prevalence per group ----------------------------------------------
# "Unique" should mean reliably present in one group and absent from the
# other. Presence in a single sample is usually noise — cross-contamination,
# index hopping, or a stray read — so require a minimum prevalence.
MIN_PREV <- 2   # samples within a group

prevalence <- psmelt(ps) %>%
  group_by(OTU, group) %>%
  summarise(n_present = sum(Abundance > 0), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = n_present, values_fill = 0)

# --- 2. Set logic ----------------------------------------------------------
unique_a <- prevalence %>% filter(group_a >= MIN_PREV, group_b == 0) %>% pull(OTU)
unique_b <- prevalence %>% filter(group_b >= MIN_PREV, group_a == 0) %>% pull(OTU)
shared   <- prevalence %>% filter(group_a > 0, group_b > 0) %>% pull(OTU)

message(sprintf("Unique to A: %d | Unique to B: %d | Shared: %d",
                length(unique_a), length(unique_b), length(shared)))

# Caveat worth stating in any results section: absence here is absence from
# *this sample set at this sequencing depth*. Unequal group sizes or depths
# inflate the "unique" count for the better-sampled group. If the groups are
# uneven, rarefy or subsample to equal n before making this claim.

# --- 3. Who are they? ------------------------------------------------------
# A count of unique features is not a finding. Their identity is.
unique_tax <- as.data.frame(tax_table(ps)) %>%
  rownames_to_column("OTU") %>%
  filter(OTU %in% c(unique_a, unique_b)) %>%
  mutate(exclusive_to = if_else(OTU %in% unique_a, "group_a", "group_b")) %>%
  left_join(
    psmelt(transform_sample_counts(ps, function(x) x / sum(x))) %>%
      group_by(OTU) %>%
      summarise(mean_rel_abund = mean(Abundance), .groups = "drop"),
    by = "OTU"
  ) %>%
  arrange(exclusive_to, desc(mean_rel_abund))

write_csv(unique_tax, "outputs/unique_taxa.csv")

# Are the unique features rare, or genuinely abundant where they occur?
unique_tax %>%
  group_by(exclusive_to) %>%
  summarise(n = n(),
            median_abund = median(mean_rel_abund),
            max_abund    = max(mean_rel_abund), .groups = "drop")

# --- 4. Which lineages are they drawn from? --------------------------------
unique_tax %>%
  count(exclusive_to, Phylum, sort = TRUE) %>%
  group_by(exclusive_to) %>%
  slice_max(n, n = 10)

# ===========================================================================
# B. Characterising an arbitrary subset
# ===========================================================================
# Same pattern for any slice: taxa unclassified below family, taxa above an
# abundance threshold, taxa in a target genus. Define the subset, then ask how
# large it is, where it sits, and what it contains.

# --- 1. Define -------------------------------------------------------------
ps_subset <- subset_taxa(ps, is.na(Genus))   # e.g. the unclassified fraction

message(sprintf("Subset: %d of %d features (%.1f%%)",
                ntaxa(ps_subset), ntaxa(ps),
                100 * ntaxa(ps_subset) / ntaxa(ps)))

# --- 2. How much of each sample does it account for? -----------------------
# Feature count and read share are different numbers, and the second is the
# one that matters. A subset can be half the features and 2% of the reads.
subset_share <- tibble(
  sample_id     = sample_names(ps),
  subset_reads  = sample_sums(ps_subset),
  total_reads   = sample_sums(ps),
  proportion    = subset_reads / total_reads
) %>%
  left_join(as(sample_data(ps), "data.frame") %>% rownames_to_column("sample_id"),
            by = "sample_id")

write_csv(subset_share, "outputs/subset_summary.csv")

# --- 3. Does the subset differ between groups? -----------------------------
wilcox.test(proportion ~ group, data = subset_share, exact = FALSE)

p_subset <- ggplot(subset_share, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.15, size = 2) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion of reads in subset")

saveRDS(p_subset, "outputs/plot_subset.rds")

# --- 4. What is in it, as far as classification goes? ----------------------
# Even an "unclassified" subset is usually resolved at a higher rank.
as.data.frame(tax_table(ps_subset)) %>%
  count(Phylum, Family, sort = TRUE) %>%
  head(20)
