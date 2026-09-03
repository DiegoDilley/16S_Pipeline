# ---------------------------------------------------------------------------
# 02 — Alpha diversity
#
# Input : outputs/ps_rarefied.rds
# Output: outputs/alpha_table.csv, outputs/alpha_tests.csv
#
# Within-sample diversity. Three metrics answering different questions:
#   Observed  how many features (richness) — most depth-sensitive
#   Shannon   richness weighted by evenness
#   Simpson   dominance; ~ probability two random reads differ
# ---------------------------------------------------------------------------

library(phyloseq)
library(tidyverse)

ps_rare <- readRDS("outputs/ps_rarefied.rds")

# --- 1. Estimate -----------------------------------------------------------
# estimate_richness needs integer counts. Passing relative abundance returns
# nonsense without warning — Observed becomes "number of non-zero proportions",
# which is not richness at any defined depth.
alpha <- estimate_richness(ps_rare, measures = c("Observed", "Shannon", "Simpson"))

# Rejoin metadata. estimate_richness returns rownames only, and mangles IDs
# containing hyphens, so match on cleaned names rather than assuming row order.
alpha <- alpha %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-")) %>%
  left_join(
    as(sample_data(ps_rare), "data.frame") %>% rownames_to_column("sample_id"),
    by = "sample_id"
  )

stopifnot(!any(is.na(alpha$group)))   # a failure here means the join missed
write_csv(alpha, "outputs/alpha_table.csv")

# --- 2. Choose a test ------------------------------------------------------
# Check normality per group rather than assuming it. Diversity indices are
# often skewed, and amplicon studies are usually small enough that a
# non-parametric test costs little power.
alpha %>%
  group_by(group) %>%
  summarise(
    n = n(),
    shapiro_p = shapiro.test(Shannon)$p.value,
    .groups = "drop"
  )

# Also check variance homogeneity if you were considering a t-test / ANOVA.
bartlett.test(Shannon ~ group, data = alpha)

# --- 3. Test ---------------------------------------------------------------
# Two groups -> Wilcoxon rank-sum. More than two -> Kruskal-Wallis, then
# pairwise Wilcoxon with multiple-testing correction.
metrics <- c("Observed", "Shannon", "Simpson")

tests <- map_dfr(metrics, function(m) {
  f <- as.formula(paste(m, "~ group"))
  if (nlevels(alpha$group) == 2) {
    tt <- wilcox.test(f, data = alpha, exact = FALSE)
    tibble(metric = m, test = "Wilcoxon", statistic = tt$statistic, p = tt$p.value)
  } else {
    tt <- kruskal.test(f, data = alpha)
    tibble(metric = m, test = "Kruskal-Wallis", statistic = tt$statistic, p = tt$p.value)
  }
})

# Correct across metrics — three tests on the same samples is three chances at
# a false positive.
tests <- tests %>% mutate(p_adj = p.adjust(p, method = "BH"))
write_csv(tests, "outputs/alpha_tests.csv")
print(tests)

# Pairwise follow-up, only if there are more than two groups and the omnibus
# test was significant.
if (nlevels(alpha$group) > 2) {
  pairwise.wilcox.test(alpha$Shannon, alpha$group, p.adjust.method = "BH")
}

# --- 4. Plot ---------------------------------------------------------------
# Boxplot plus the underlying points. With n < 10 per group a boxplot alone
# implies more precision than the data supports.
alpha_long <- alpha %>%
  pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = metrics))

p_alpha <- ggplot(alpha_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~ metric, scales = "free_y") +   # free_y: the metrics share no scale
  labs(x = NULL, y = "Diversity estimate")

saveRDS(p_alpha, "outputs/plot_alpha.rds")   # themed and written out in script 08
