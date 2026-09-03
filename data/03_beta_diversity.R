# ---------------------------------------------------------------------------
# 03 — Beta diversity: ordination and PERMANOVA
#
# Input : outputs/ps_filtered.rds
# Output: outputs/permanova.csv, outputs/plot_nmds.rds
#
# Between-sample differences. The sequence is: transform, distance, ordinate,
# test, then check that the test means what you think it means.
# ---------------------------------------------------------------------------

library(phyloseq)
library(vegan)
library(tidyverse)
library(pairwiseAdonis)

set.seed(42)
ps <- readRDS("outputs/ps_filtered.rds")

# --- 1. Transform ----------------------------------------------------------
# Bray-Curtis is abundance-weighted, so uneven sequencing depth alone will
# separate samples. Convert to proportions first.
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# --- 2. Distance -----------------------------------------------------------
# bray     abundance-weighted, no phylogeny — the default choice
# jaccard  presence/absence; use binary = TRUE, or rare taxa are downweighted
# unifrac  phylogeny-aware; requires a tree in the phyloseq object
dist_bc <- phyloseq::distance(ps_rel, method = "bray")

meta <- as(sample_data(ps_rel), "data.frame")

# --- 3. Ordinate -----------------------------------------------------------
# NMDS is rank-based and does not assume a linear response, which suits
# ecological distances. It is iterative, so set a seed and report stress.
ord_nmds <- metaMDS(dist_bc, k = 2, trymax = 100, trace = FALSE)

# Stress is the headline diagnostic. Rules of thumb:
#   < 0.05 excellent | < 0.10 good | < 0.20 usable | > 0.20 do not interpret
# Above 0.20, try k = 3 or switch to PCoA.
message("NMDS stress: ", round(ord_nmds$stress, 3))
stressplot(ord_nmds)   # Shepard plot: observed vs ordination distance

# PCoA alternative — eigenvalue-based, so axes carry a % variance explained
# that NMDS axes do not have.
ord_pcoa <- ordinate(ps_rel, method = "PCoA", distance = "bray")
head(ord_pcoa$values$Relative_eig)

# --- 4. PERMANOVA ----------------------------------------------------------
# Tests whether group centroids differ in multivariate space. by = "terms"
# gives sequential sums of squares, so with more than one predictor the order
# of terms in the formula changes the result. State the order you used.
perm <- adonis2(
  dist_bc ~ group,
  data = meta,
  permutations = 999,
  by = "terms"
)
print(perm)

# R2 is the effect size and matters more than the p-value. Enough samples will
# make almost anything significant; R2 says how much variation it explains.
write.csv(as.data.frame(perm), "outputs/permanova.csv")

# Nested or repeated-measures designs: restrict permutations rather than
# treating samples as independent.
# adonis2(dist_bc ~ group, data = meta, permutations = 999,
#         strata = meta$site)

# --- 5. Dispersion check ---------------------------------------------------
# The step people skip, and the reason PERMANOVA results get over-read.
# A significant adonis2 means "centroids differ, or spread differs, or both".
# If dispersion also differs, you cannot attribute the result to location.
disp <- betadisper(dist_bc, meta$group)
permutest(disp, permutations = 999)
plot(disp)

# Reporting: a significant PERMANOVA with a non-significant betadisper is a
# clean location effect. Both significant means you say so.

# --- 6. Pairwise, for >2 groups --------------------------------------------
if (nlevels(meta$group) > 2) {
  pw <- pairwise.adonis2(dist_bc ~ group, data = meta, permutations = 999)
  print(pw)
}

# --- 7. Plot ---------------------------------------------------------------
scores_df <- as.data.frame(vegan::scores(ord_nmds, display = "sites")) %>%
  rownames_to_column("sample_id") %>%
  left_join(meta %>% rownames_to_column("sample_id"), by = "sample_id")

p_nmds <- ggplot(scores_df, aes(x = NMDS1, y = NMDS2, colour = group)) +
  geom_point(size = 3, alpha = 0.9) +
  # 95% confidence ellipse around the group mean. Needs n >= 4 per group;
  # below that it is drawn from too little information to mean much.
  stat_ellipse(aes(group = group), type = "t", level = 0.95, linewidth = 0.5) +
  # Put the diagnostics on the figure, not only in the text.
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.8, size = 3,
           label = sprintf("Stress = %.3f\nPERMANOVA R2 = %.3f, p = %.3f",
                           ord_nmds$stress, perm$R2[1], perm$`Pr(>F)`[1])) +
  coord_equal() +   # NMDS axes are in the same units; do not distort them
  labs(x = "NMDS1", y = "NMDS2")

saveRDS(p_nmds, "outputs/plot_nmds.rds")
