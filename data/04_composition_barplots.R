# ---------------------------------------------------------------------------
# 04 — Taxonomic composition
#
# Input : outputs/ps_filtered.rds
# Output: outputs/composition_table.csv, outputs/plot_composition.rds
#
# The pattern generalises to any rank: agglomerate, convert to proportions,
# keep the top N, lump the rest, then plot.
# ---------------------------------------------------------------------------

library(phyloseq)
library(tidyverse)

ps <- readRDS("outputs/ps_filtered.rds")

RANK  <- "Phylum"   # swap for Class, Family, Genus...
TOP_N <- 10         # beyond ~12, colours stop being distinguishable

# --- 1. Agglomerate --------------------------------------------------------
# NArm = FALSE keeps features unclassified at this rank. Dropping them makes
# the remaining proportions add to 1 while quietly hiding part of the
# community, which is how "100% assigned" plots happen.
ps_rank <- tax_glom(ps, taxrank = RANK, NArm = FALSE)

# --- 2. Relative abundance -------------------------------------------------
# Per sample, after agglomeration.
ps_rel <- transform_sample_counts(ps_rank, function(x) x / sum(x))

# --- 3. Long format --------------------------------------------------------
# psmelt joins the OTU table, taxonomy and metadata into one long data frame.
# It is slow on large objects — agglomerate first, as above.
melted <- psmelt(ps_rel) %>%
  mutate(taxon = replace_na(as.character(.data[[RANK]]), "Unclassified"))

# --- 4. Top N + Other ------------------------------------------------------
# Rank by mean relative abundance across samples, not by total, so that
# deeply sequenced samples don't decide the legend.
top_taxa <- melted %>%
  group_by(taxon) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  slice_max(mean_abund, n = TOP_N) %>%
  pull(taxon)

melted <- melted %>%
  mutate(taxon_lumped = if_else(taxon %in% top_taxa, taxon, "Other"))

# Order the stack by abundance, with Other and Unclassified pinned last so
# they don't float around mid-stack between panels.
taxon_order <- melted %>%
  filter(!taxon_lumped %in% c("Other", "Unclassified")) %>%
  group_by(taxon_lumped) %>%
  summarise(m = mean(Abundance), .groups = "drop") %>%
  arrange(desc(m)) %>%
  pull(taxon_lumped)

melted <- melted %>%
  mutate(taxon_lumped = factor(taxon_lumped,
                               levels = c(taxon_order, "Unclassified", "Other")))

# --- 5. Per-sample bars ----------------------------------------------------
p_samples <- melted %>%
  ggplot(aes(x = Sample, y = Abundance, fill = taxon_lumped)) +
  geom_col(width = 0.9) +
  # space = "free_x" keeps bar widths equal when groups have different n.
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Relative abundance", fill = RANK)

# --- 6. Per-group bars -----------------------------------------------------
# Mean of the per-sample proportions, not the proportion of pooled counts.
# Pooling lets the deepest-sequenced sample dominate its whole group.
group_summary <- melted %>%
  group_by(group, taxon_lumped) %>%
  summarise(mean_abund = mean(Abundance),
            sd_abund   = sd(Abundance), .groups = "drop")

write_csv(group_summary, "outputs/composition_table.csv")

p_groups <- ggplot(group_summary,
                   aes(x = group, y = mean_abund, fill = taxon_lumped)) +
  geom_col(width = 0.6) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Mean relative abundance", fill = RANK)

saveRDS(p_samples, "outputs/plot_composition.rds")
saveRDS(p_groups,  "outputs/plot_composition_grouped.rds")

# --- Note on stacked bars --------------------------------------------------
# Only the bottom segment sits on a common baseline, so everything above it is
# hard to compare across bars. For a handful of taxa you care about, a faceted
# boxplot of relative abundance per group is easier to read and shows the
# spread that a stacked bar hides.
