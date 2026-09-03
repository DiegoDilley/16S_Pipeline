# ---------------------------------------------------------------------------
# 07 — Functional prediction with FAPROTAX
#
# Input : outputs/ps_filtered.rds
# Output: outputs/functional_table.csv, outputs/plot_functions.rds
#
# FAPROTAX maps taxonomic assignments to metabolic functions using a curated
# database of cultured representatives. It is a lookup table, not a genomic
# prediction: if a genus is in the database, its known functions are assigned
# to every feature classified to it.
#
# What that means in practice:
#   - only well-studied lineages get annotated; expect 20-40% of reads assigned
#   - it assumes an uncultured relative behaves like its cultured cousin
#   - results are hypotheses to test, not measured activity
# Say this in the methods. Reviewers ask.
# ---------------------------------------------------------------------------

library(phyloseq)
library(tidyverse)

ps <- readRDS("outputs/ps_filtered.rds")

# ===========================================================================
# Route A — the standalone Python script
# ===========================================================================
# FAPROTAX ships as collapse_table.py plus a database file. This route is
# closest to the published method and pins an explicit database version.

dir.create("outputs/faprotax", showWarnings = FALSE, recursive = TRUE)

# --- 1. Export in the expected format --------------------------------------
# A classic-format table: features as rows, samples as columns, with a final
# "taxonomy" column holding the semicolon-delimited lineage.
tax_string <- as.data.frame(tax_table(ps)) %>%
  mutate(across(everything(), ~ replace_na(.x, ""))) %>%
  transmute(taxonomy = paste(Kingdom, Phylum, Class, Order,
                             Family, Genus, Species, sep = ";"))

export <- as.data.frame(as(otu_table(ps), "matrix"))
if (!taxa_are_rows(ps)) export <- as.data.frame(t(export))
export$taxonomy <- tax_string$taxonomy[match(rownames(export), rownames(tax_string))]

out_path <- "outputs/faprotax/feature_table_tax.tsv"
writeLines("# Constructed from biom file", out_path)
write.table(export, out_path, sep = "\t", quote = FALSE,
            col.names = NA, append = TRUE)

# --- 2. Run it -------------------------------------------------------------
# system2 rather than system(): arguments are passed as a vector, so paths
# with spaces don't need manual quoting, and the exit status comes back.
status <- system2(
  "python",
  args = c(
    "tools/FAPROTAX/collapse_table.py",
    "-i", out_path,
    "-o", "outputs/faprotax/functional_table.tsv",
    "-g", "tools/FAPROTAX/FAPROTAX.txt",
    "--collapse_by_metadata", "taxonomy",
    "--group_leftovers_as", "Unassigned",
    "-r", "outputs/faprotax/report.txt",   # keep this: it logs what matched
    "-v"
  )
)
if (status != 0) stop("FAPROTAX failed — check the report file")

# The report gives the assignment rate. Quote it in the results.
# readLines("outputs/faprotax/report.txt") %>% head(20)

func <- read.table("outputs/faprotax/functional_table.tsv",
                   sep = "\t", header = TRUE, row.names = 1,
                   skip = 1, comment.char = "", check.names = FALSE)

# ===========================================================================
# Route B — microeco, entirely in R
# ===========================================================================
# Bundles the same database. Convenient, but the version is whatever the
# installed package ships; record it if you use this route.
#
# library(microeco)
# me <- microeco::microtable$new(
#   otu_table  = as.data.frame(as(otu_table(ps), "matrix")),
#   tax_table  = as.data.frame(tax_table(ps)),
#   sample_table = as(sample_data(ps), "data.frame")
# )
# me$tidy_dataset()
# tf <- trans_func$new(me)
# tf$cal_spe_func(prok_database = "FAPROTAX")
# tf$cal_spe_func_perc(abundance_weighted = TRUE)
# func <- t(tf$res_spe_func_perc)

# ===========================================================================
# Downstream — identical either way
# ===========================================================================

# --- 3. Relative abundance and tidy ----------------------------------------
func_rel <- sweep(as.matrix(func), 2, colSums(as.matrix(func)), "/")

meta <- as(sample_data(ps), "data.frame") %>% rownames_to_column("sample_id")

func_long <- as.data.frame(func_rel) %>%
  rownames_to_column("function_group") %>%
  pivot_longer(-function_group, names_to = "sample_id", values_to = "abundance") %>%
  left_join(meta, by = "sample_id")

write_csv(func_long, "outputs/functional_table.csv")

# --- 4. Which functions are worth looking at? ------------------------------
# Most rows will be near zero. Filter before testing, so multiple-testing
# correction isn't spent on functions that were never detected.
keep_funcs <- func_long %>%
  group_by(function_group) %>%
  summarise(mean_abund = mean(abundance),
            prevalence = mean(abundance > 0), .groups = "drop") %>%
  filter(mean_abund > 0.001, prevalence > 0.25) %>%
  pull(function_group)

# --- 5. Test between groups ------------------------------------------------
func_tests <- func_long %>%
  filter(function_group %in% keep_funcs) %>%
  group_by(function_group) %>%
  summarise(
    p = wilcox.test(abundance ~ group, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(p_adj)

print(head(func_tests, 15))

# --- 6. Plot ---------------------------------------------------------------
top_funcs <- func_long %>%
  filter(function_group %in% keep_funcs) %>%
  group_by(function_group) %>%
  summarise(m = mean(abundance), .groups = "drop") %>%
  slice_max(m, n = 15) %>%
  pull(function_group)

p_func <- func_long %>%
  filter(function_group %in% top_funcs) %>%
  group_by(function_group, group) %>%
  summarise(mean_abund = mean(abundance), .groups = "drop") %>%
  ggplot(aes(x = mean_abund,
             y = reorder(function_group, mean_abund),
             fill = group)) +
  geom_col(position = "dodge") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Mean relative abundance of assigned reads",
       y = NULL, fill = NULL)

saveRDS(p_func, "outputs/plot_functions.rds")
