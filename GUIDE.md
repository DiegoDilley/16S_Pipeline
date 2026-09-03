# Amplicon Analysis in R: A Walkthrough

This guide covers what happens to 16S/18S amplicon data after denoising, from a
raw feature table through to publication figures.

This is for refrence not runnable code.

**Contents**

1. [Where the data comes from](#1-where-the-data-comes-from)
2. [Setting up](#2-setting-up)
3. [Building a phyloseq object](#3-building-a-phyloseq-object)
4. [Filtering](#4-filtering)
5. [The rarefaction decision](#5-the-rarefaction-decision)
6. [Alpha diversity](#6-alpha-diversity)
7. [Beta diversity](#7-beta-diversity)
8. [Taxonomic composition](#8-taxonomic-composition)
9. [Group-specific taxa](#9-group-specific-taxa)
10. [Source tracking](#10-source-tracking)
11. [Functional prediction](#11-functional-prediction)
12. [Figures](#12-figures)
13. [Things that go wrong](#13-things-that-go-wrong)

Throughout, the example study compares two groups, `group_a` and `group_b`,
stored in a metadata column called `group`.

---

## 1. Where the data comes from

R is not the start of this pipeline. Before any of it, you've run something like
QIIME2 or DADA2 in a shell (i didnt record my code for that so its not here) : 

- imported and demultiplexed raw reads
- trimmed and quality-filtered them
- denoised into ASVs with Deblur or DADA2
- classified those ASVs against a reference database like SILVA

That produces three things, which are the entire input to everything below:

| File | What it is |
|---|---|
| `feature_table.tsv` | Features (ASVs) by samples, containing read counts |
| `taxonomy.tsv` | Each feature ID mapped to a semicolon-delimited lineage |
| `metadata.tsv` | Each sample ID mapped to its group and any covariates |

If your feature table came out of QIIME2 as a `.qza`, you export it to biom and
then to TSV. Worth knowing: a `.qza` is just a zip file carrying its own
provenance, so if you've lost the script that produced it, the parameters are
still recoverable:

```bash
unzip -o table.qza -d table_prov
cat table_prov/*/provenance/action/action.yaml
```

That file records the plugin, the action, and every parameter. It's the
authoritative answer to "what trim length did I use", which is a question you
will eventually be asked.

---

## 2. Setting up

Three of the packages used here don't come from CRAN, and a README that says
"install the required packages" will fail for anyone who clones your repo. List
the install lines explicitly.

```r
# CRAN
install.packages(c("tidyverse", "vegan", "remotes", "BiocManager",
                   "patchwork", "RColorBrewer", "ggrepel"))

# Bioconductor — phyloseq is not on CRAN, which is the usual first stumble
BiocManager::install("phyloseq")

# GitHub
remotes::install_github("cozygene/FEAST")
remotes::install_github("ChiLiubio/microeco")
remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
```

FAPROTAX (section 11) isn't an R package at all. It's a Python script you
download separately.

Capture your environment once and commit the result. Two years later it answers
"which version of vegan produced this number".

```r
writeLines(capture.output(sessionInfo()), "outputs/sessionInfo.txt")
```

And set a seed at the top of anything involving randomness: NMDS, permutation
tests, rarefaction, and the EM algorithm in FEAST all give different answers
between runs otherwise.

```r
set.seed(42)
```

---

## 3. Building a phyloseq object

`phyloseq` bundles the count table, taxonomy and metadata into a single object
so that subsetting a sample subsets all three together. Getting the three tables
in cleanly is most of the work.

### Reading the feature table

```r
counts <- read.table(
  "data/feature_table.tsv",
  sep = "\t", header = TRUE, row.names = 1,
  skip = 1, comment.char = "", check.names = FALSE
)
counts <- as.matrix(counts)
```

Two arguments matter more than they look:

`skip = 1` drops the `# Constructed from biom file` line that QIIME2 puts at the
top of an exported table. Without it, that line becomes your header.

`check.names = FALSE` stops R prepending an `X` to sample IDs that start with a
digit. R does this silently, and the result is that your sample IDs no longer
match your metadata, so phyloseq quietly keeps zero samples.

### Reading the taxonomy

The taxonomy file has one column of lineage strings that need splitting into
ranks.

```r
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
  mutate(across(Kingdom:Species, ~ str_remove(.x, "^[a-z]__"))) %>%
  mutate(across(Kingdom:Species, ~ na_if(str_trim(.x), "")))

taxonomy <- as.matrix(taxonomy[, c("Kingdom", "Phylum", "Class",
                                   "Order", "Family", "Genus", "Species")])
```

`fill = "right"` handles features that were only classified to family or above,
padding the remaining ranks rather than erroring.

The `str_remove` strips the `d__` / `p__` rank prefixes that SILVA and
GreenGenes attach. The `na_if` converts empty strings to real `NA`, which
matters later: `subset_taxa()` treats `""` and `NA` differently, and an empty
string will not behave like missing data when you filter.

### Reading the metadata

```r
metadata <- read.table("data/metadata.tsv", sep = "\t",
                       header = TRUE, row.names = 1, comment.char = "")

metadata$group <- factor(metadata$group, levels = c("group_a", "group_b"))
```

Set factor levels explicitly. If you don't, R orders them alphabetically, and
that order silently propagates into every plot legend and every model contrast
for the rest of the analysis.

### Assembling and checking

```r
ps_raw <- phyloseq(
  otu_table(counts, taxa_are_rows = TRUE),
  tax_table(taxonomy),
  sample_data(metadata)
)

message("Samples in metadata: ", nrow(metadata))
message("Samples in object:   ", nsamples(ps_raw))
message("Features:            ", ntaxa(ps_raw))
```

Do check those numbers. phyloseq intersects the three tables without
complaining, so a sample-ID mismatch doesn't raise an error. It just gives you
an object with fewer samples than you expected, and you find out at the figure
stage.

---

## 4. Filtering

Order matters here. Drop non-target lineages first, then filter on prevalence,
so prevalence is calculated on the features you're actually keeping.

### Non-target lineages

Chloroplast and mitochondrial sequences amplify with 16S primers but aren't
bacteria you're interested in. Unassigned features at kingdom level are usually
noise.

```r
ps <- subset_taxa(
  ps_raw,
  !is.na(Kingdom) & Kingdom != "Unassigned" &
    (is.na(Order)  | Order  != "Chloroplast") &
    (is.na(Family) | Family != "Mitochondria")
)
```

The `is.na()` guards are essential and easy to leave out. Comparing `NA` to a
string yields `NA`, and `subset_taxa` drops `NA` rows. So writing
`Order != "Chloroplast"` on its own removes chloroplasts *and* every feature
that was never classified to order level, which can be a large fraction of your
data.

### Prevalence

Features seen once, in one sample, are mostly index hopping and sequencing
error. They inflate feature counts and add noise to distance matrices without
adding signal.

```r
min_count   <- 2
min_samples <- 3
keep <- genefilter_sample(ps, filterfun_sample(function(x) x >= min_count),
                          A = min_samples)
ps <- prune_taxa(keep, ps)
```

The thresholds are a judgement call. Justify them in your methods rather than
presenting them as standard, because they aren't.

### Sample depth

```r
sort(sample_sums(ps))
ps <- prune_samples(sample_sums(ps) >= 1000, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
```

Look at the sorted distribution before choosing a cutoff. Usually there's a
visible gap between a handful of failed libraries and everything else, and the
cutoff belongs in that gap rather than at a round number you picked in advance.

The final line removes features left at zero once those samples are gone.

---

## 5. The rarefaction decision

Rarefying means subsampling every library down to the same depth. It's
contested, because it throws away real data, and it's a poor choice for
differential abundance testing. But sequencing depth strongly affects richness
estimates, and rarefying is still the common way to handle that.

A workable split, and the one used in the rest of this guide:

| Analysis | Input |
|---|---|
| Alpha diversity | Rarefied counts |
| Beta diversity | Relative abundance, not rarefied |
| Composition | Relative abundance |
| Differential abundance | Raw counts, into a method that models depth |

```r
depth <- min(sample_sums(ps))

ps_rare <- rarefy_even_depth(
  ps,
  sample.size = depth,
  rngseed = 42,
  replace = FALSE,
  verbose = FALSE
)
```

`rngseed` is not optional if you want the same numbers twice.

Whatever you choose, state it in the methods and apply it consistently. The
choice is defensible either way; being inconsistent about it isn't.

---

## 6. Alpha diversity

Alpha diversity is within-sample diversity. The three standard metrics answer
different questions:

- **Observed** counts features. It's the most sensitive to sequencing depth,
  which is the reason for rarefying.
- **Shannon** weights richness by evenness. A community with one dominant taxon
  scores lower than an even one with the same number of taxa.
- **Simpson** measures dominance, roughly the probability that two randomly
  drawn reads belong to different taxa.

### Estimating

```r
alpha <- estimate_richness(ps_rare, measures = c("Observed", "Shannon", "Simpson"))
```

This needs integer counts. Passing relative abundance returns numbers without
warning, but "Observed" then means "count of non-zero proportions", which isn't
richness at any defined depth.

Rejoining the metadata is fiddlier than it should be, because
`estimate_richness` mangles sample IDs containing hyphens:

```r
alpha <- alpha %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_replace_all(sample_id, "\\.", "-")) %>%
  left_join(
    as(sample_data(ps_rare), "data.frame") %>% rownames_to_column("sample_id"),
    by = "sample_id"
  )

stopifnot(!any(is.na(alpha$group)))
```

That `stopifnot` catches a failed join immediately rather than at the plotting
stage.

### Choosing a test

Check rather than assume. Diversity indices are often skewed, and amplicon
studies are usually small enough that a non-parametric test costs little power.

```r
alpha %>%
  group_by(group) %>%
  summarise(n = n(), shapiro_p = shapiro.test(Shannon)$p.value, .groups = "drop")

bartlett.test(Shannon ~ group, data = alpha)
```

With two groups, use a Wilcoxon rank-sum test. With more, use Kruskal-Wallis
followed by pairwise Wilcoxon tests if the omnibus test is significant.

```r
wilcox.test(Shannon ~ group, data = alpha, exact = FALSE)

# more than two groups
kruskal.test(Shannon ~ group, data = alpha)
pairwise.wilcox.test(alpha$Shannon, alpha$group, p.adjust.method = "BH")
```

Correct across metrics too. Testing three indices on the same samples is three
chances at a false positive.

```r
tests <- tests %>% mutate(p_adj = p.adjust(p, method = "BH"))
```

### Plotting

```r
p_alpha <- alpha_long %>%
  ggplot(aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = NULL, y = "Diversity estimate")
```

Show the points. With fewer than about ten samples per group, a boxplot on its
own implies more precision than the data supports. `scales = "free_y"` because
the three metrics share no common scale.

---

## 7. Beta diversity

Beta diversity is between-sample difference. The sequence is: transform,
calculate distances, ordinate, test, then check that the test means what you
think it means.

### Transform, then distance

```r
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
dist_bc <- phyloseq::distance(ps_rel, method = "bray")
```

Bray-Curtis is abundance-weighted, so if you feed it raw counts, uneven
sequencing depth alone will separate your samples. Convert to proportions first.

Alternatives: Jaccard for presence/absence (pass `binary = TRUE`, otherwise rare
taxa get downweighted anyway), or UniFrac if you have a phylogenetic tree in the
object.

### Ordination

```r
ord_nmds <- metaMDS(dist_bc, k = 2, trymax = 100, trace = FALSE)
message("NMDS stress: ", round(ord_nmds$stress, 3))
stressplot(ord_nmds)
```

NMDS is rank-based, so it doesn't assume a linear response, which suits
ecological distances. It's iterative, hence the seed and the `trymax`.

Stress is the diagnostic that decides whether the plot is interpretable at all:

| Stress | Reading |
|---|---|
| < 0.05 | Excellent |
| < 0.10 | Good |
| < 0.20 | Usable |
| > 0.20 | Don't interpret it |

Above 0.20, try three dimensions or switch to PCoA. PCoA is eigenvalue-based,
so its axes carry a percentage of variance explained, which NMDS axes do not.

```r
ord_pcoa <- ordinate(ps_rel, method = "PCoA", distance = "bray")
head(ord_pcoa$values$Relative_eig)
```

### PERMANOVA

An ordination is a picture. PERMANOVA is the test.

```r
perm <- adonis2(dist_bc ~ group, data = meta, permutations = 999, by = "terms")
```

Two things to be careful about.

`by = "terms"` gives sequential sums of squares, so with more than one predictor
the order of terms in the formula changes the result. Report the order you used.

`R2` is the effect size, and it matters more than the p-value. With enough
samples almost anything reaches significance; R2 tells you how much of the
variation the grouping actually explains.

For nested or repeated-measures designs, restrict the permutations instead of
treating samples as independent:

```r
adonis2(dist_bc ~ group, data = meta, permutations = 999, strata = meta$site)
```

### The dispersion check

This is the step most often skipped, and it's why PERMANOVA results get
over-interpreted.

A significant `adonis2` means the groups differ in centroid position, *or* in
spread, or both. If the groups have different within-group variability, you
can't attribute the result to a location difference.

```r
disp <- betadisper(dist_bc, meta$group)
permutest(disp, permutations = 999)
plot(disp)
```

A significant PERMANOVA with a non-significant `betadisper` is a clean location
effect. If both are significant, say so, and describe the result as a difference
in composition and heterogeneity.

### Plotting

```r
p_nmds <- ggplot(scores_df, aes(x = NMDS1, y = NMDS2, colour = group)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = group), type = "t", level = 0.95, linewidth = 0.5) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.8, size = 3,
           label = sprintf("Stress = %.3f\nPERMANOVA R2 = %.3f, p = %.3f",
                           ord_nmds$stress, perm$R2[1], perm$`Pr(>F)`[1])) +
  coord_equal()
```

`coord_equal()` matters: NMDS axes are in the same units, and stretching one
distorts the distances the whole plot represents.

Put the stress value and the PERMANOVA result on the figure. Confidence
ellipses need at least four samples per group to mean much.

---

## 8. Taxonomic composition

The pattern is the same at any rank: agglomerate, convert to proportions, keep
the top N, lump the rest.

```r
RANK  <- "Phylum"
TOP_N <- 10

ps_rank <- tax_glom(ps, taxrank = RANK, NArm = FALSE)
ps_rel  <- transform_sample_counts(ps_rank, function(x) x / sum(x))
melted  <- psmelt(ps_rel) %>%
  mutate(taxon = replace_na(as.character(.data[[RANK]]), "Unclassified"))
```

`NArm = FALSE` keeps features that weren't classified at this rank. Dropping
them makes the remaining proportions sum neatly to 1 while hiding part of the
community, and that's how plots claiming 100% classification happen.

`psmelt` joins everything into one long data frame. It's slow on large objects,
so agglomerate first.

### Top N and Other

```r
top_taxa <- melted %>%
  group_by(taxon) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  slice_max(mean_abund, n = TOP_N) %>%
  pull(taxon)

melted <- melted %>%
  mutate(taxon_lumped = if_else(taxon %in% top_taxa, taxon, "Other"))
```

Rank by *mean* relative abundance across samples rather than by total. Using the
total lets your most deeply sequenced sample decide the legend.

Then fix the stacking order, pinning the residual categories last so they don't
drift around mid-stack between panels:

```r
melted <- melted %>%
  mutate(taxon_lumped = factor(taxon_lumped,
                               levels = c(taxon_order, "Unclassified", "Other")))
```

### The plot

```r
p_samples <- ggplot(melted, aes(x = Sample, y = Abundance, fill = taxon_lumped)) +
  geom_col(width = 0.9) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Relative abundance", fill = RANK)
```

`space = "free_x"` keeps bar widths consistent when the groups have different
sample counts.

For a per-group version, average the per-sample proportions rather than pooling
counts and taking proportions of the pool. Pooling gives the deepest-sequenced
sample control of its group's bar.

```r
group_summary <- melted %>%
  group_by(group, taxon_lumped) %>%
  summarise(mean_abund = mean(Abundance), sd_abund = sd(Abundance), .groups = "drop")
```

One caveat on stacked bars generally: only the bottom segment sits on a common
baseline, so everything above it is hard to compare across bars. If you care
about a handful of specific taxa, a faceted boxplot of their relative abundance
per group is easier to read and shows the spread that stacking hides.

---

## 9. Group-specific taxa

Two related questions: which features occur in one group and not the other, and
what a given slice of the community looks like.

### Presence and absence

```r
MIN_PREV <- 2

prevalence <- psmelt(ps) %>%
  group_by(OTU, group) %>%
  summarise(n_present = sum(Abundance > 0), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = n_present, values_fill = 0)

unique_a <- prevalence %>% filter(group_a >= MIN_PREV, group_b == 0) %>% pull(OTU)
unique_b <- prevalence %>% filter(group_b >= MIN_PREV, group_a == 0) %>% pull(OTU)
shared   <- prevalence %>% filter(group_a > 0, group_b > 0) %>% pull(OTU)
```

Requiring a minimum prevalence matters. A feature present in exactly one sample
of one group is usually cross-contamination or index hopping, not biology.

Be explicit about what "unique" means when you write it up: absent from this
sample set at this sequencing depth. Unequal group sizes or depths will inflate
the unique count for whichever group was better sampled. If the groups are
uneven, subsample to equal n before making the claim.

### Who they are

A count of unique features isn't a finding. Their identity is.

```r
unique_tax <- as.data.frame(tax_table(ps)) %>%
  rownames_to_column("OTU") %>%
  filter(OTU %in% c(unique_a, unique_b)) %>%
  mutate(exclusive_to = if_else(OTU %in% unique_a, "group_a", "group_b")) %>%
  left_join(mean_abundances, by = "OTU") %>%
  arrange(exclusive_to, desc(mean_rel_abund))
```

Then ask whether they're rare or actually abundant where they occur, and which
lineages they come from:

```r
unique_tax %>%
  group_by(exclusive_to) %>%
  summarise(n = n(),
            median_abund = median(mean_rel_abund),
            max_abund    = max(mean_rel_abund), .groups = "drop")

unique_tax %>% count(exclusive_to, Phylum, sort = TRUE)
```

### Describing any subset

The same pattern works for any slice you can define: taxa unclassified below
family, taxa above an abundance threshold, taxa in a target genus.

```r
ps_subset <- subset_taxa(ps, is.na(Genus))

subset_share <- tibble(
  sample_id    = sample_names(ps),
  subset_reads = sample_sums(ps_subset),
  total_reads  = sample_sums(ps),
  proportion   = subset_reads / total_reads
)
```

Feature count and read share are different numbers, and the second is usually
the one that matters. A subset can easily be half your features and 2% of your
reads.

From there, test whether the subset differs between groups, and look at what
it's composed of. Even an "unclassified" subset is normally resolved at a higher
rank.

---

## 10. Source tracking

Source tracking estimates what fraction of a community came from each candidate
source. FEAST does this with a Bayesian mixture model, and reports an
**Unknown** fraction covering everything the sources don't explain.

The Unknown fraction is a result, not an error term. In many designs it's the
most informative number in the table. It's also a direct measure of how well you
sampled: anything real that you didn't collect ends up there.

### The working-directory problem

Deal with this first, because it will otherwise cost you an afternoon.

`FEAST()` calls `setwd()` internally and never restores it. Every relative path
after the first call resolves somewhere unexpected, and in a loop it compounds.

```r
PROJECT_ROOT <- getwd()

# ... FEAST() call ...

setwd(PROJECT_ROOT)   # after every single call
```

### Input format

FEAST wants **samples as rows**, which is the transpose of phyloseq's default,
and it wants raw integer counts.

```r
counts <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) counts <- t(counts)
storage.mode(counts) <- "integer"
```

The metadata needs exactly three columns, named exactly:

| Column | Contents |
|---|---|
| `Env` | Label for the source or sink type |
| `SourceSink` | Literally `"Source"` or `"Sink"` |
| `id` | Integer linking each sink to its source set; `NA` for sources when all sinks share one pool |

```r
feast_meta <- metadata %>%
  mutate(
    Env        = source_type,
    SourceSink = if_else(source_type == "sink_community", "Sink", "Source"),
    id         = if_else(SourceSink == "Sink", row_number(), NA_integer_)
  ) %>%
  select(Env, SourceSink, id)

common     <- intersect(rownames(counts), rownames(feast_meta))
counts     <- counts[common, , drop = FALSE]
feast_meta <- feast_meta[common, , drop = FALSE]

stopifnot(identical(rownames(counts), rownames(feast_meta)))
```

Rownames must match and be in the same order. A mismatch surfaces as a cryptic
subscript error from deep inside the package rather than as a useful message.

### Running

```r
feast_out <- FEAST(
  C = counts,
  metadata = feast_meta,
  different_sources_flag = 0,
  dir_path = file.path(PROJECT_ROOT, "outputs/feast"),
  outfile = "feast_run",
  EM_iterations = 1000
)

setwd(PROJECT_ROOT)
```

`different_sources_flag = 0` means all sinks draw from one shared source pool,
which is the common case. Set it to `1` when each sink has its own source set,
matched through the `id` column.

Raise `EM_iterations` if repeated runs with different seeds give unstable
proportions.

### Checking and plotting

Proportions per sink should sum to about 1. Check it:

```r
results %>%
  group_by(sink) %>%
  summarise(total = sum(proportion)) %>%
  filter(abs(total - 1) > 0.01)
```

Plot with Unknown pinned last in the factor levels, so it reads as the residual
it is.

If the Unknown fraction is large, it's worth characterising rather than
reporting as a bare number. The subset workflow from section 9 applies directly:
pull the features abundant in sinks and near-absent from every source, and
describe them. That set is either an unsampled source or a genuinely resident
community, and the taxonomy usually distinguishes the two.

---

## 11. Functional prediction

FAPROTAX maps taxonomic assignments to metabolic functions using a curated
database of cultured organisms. It's a lookup table. If a genus appears in the
database, its documented functions are assigned to every feature classified to
that genus.

That has three consequences worth stating in your methods, because reviewers ask:

- Only well-studied lineages get annotated. An assignment rate of 20-40% of
  reads is normal.
- It assumes an uncultured organism behaves like its cultured relative.
- The output is a hypothesis about potential function, not a measurement of
  activity.

### Running it

FAPROTAX is a Python script. Export a classic-format table with a `taxonomy`
column, then call it:

```r
status <- system2(
  "python",
  args = c(
    "tools/FAPROTAX/collapse_table.py",
    "-i", out_path,
    "-o", "outputs/faprotax/functional_table.tsv",
    "-g", "tools/FAPROTAX/FAPROTAX.txt",
    "--collapse_by_metadata", "taxonomy",
    "--group_leftovers_as", "Unassigned",
    "-r", "outputs/faprotax/report.txt",
    "-v"
  )
)
if (status != 0) stop("FAPROTAX failed — check the report file")
```

Use `system2` rather than `system`. Arguments pass as a vector, so paths with
spaces don't need manual quoting, and you get the exit status back.

Keep the report file. It records the assignment rate, which belongs in your
results.

The `microeco` package bundles the same database and does this without leaving
R, via `trans_func$new()` and `cal_spe_func()`. It's more convenient, but the
database version is whatever the installed package ships, so record it.

### Downstream

Most functional groups will be near zero. Filter before testing, so
multiple-testing correction isn't spent on functions that were never detected.

```r
keep_funcs <- func_long %>%
  group_by(function_group) %>%
  summarise(mean_abund = mean(abundance),
            prevalence = mean(abundance > 0), .groups = "drop") %>%
  filter(mean_abund > 0.001, prevalence > 0.25) %>%
  pull(function_group)
```

Then test the survivors between groups and correct with BH, the same as any
other multi-outcome comparison.

---

## 12. Figures

Keep figure generation separate from analysis. The analysis scripts build plot
objects and save them unthemed; one export script applies the theme and writes
every file.

```r
saveRDS(p_alpha, "outputs/plot_alpha.rds")   # in the analysis script
```

The reason is maintenance. If each script saves its own figures, the plotting
code drifts apart, the same panel ends up rendered twice under two different
numbers, and changing a font means editing seven files. With one export script,
figure numbering exists in exactly one place, which is also the only place it
can stay in sync with the manuscript.

### One theme, one palette

```r
theme_pub <- function(base_size = 10, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey92"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(face = "bold", size = rel(0.9)),
      axis.text.x      = element_text(angle = 45, hjust = 1)
    )
}
theme_set(theme_pub())

pal_group <- c(group_a = "#4E79A7", group_b = "#E15759")
```

Name the palette entries rather than relying on position. A named vector keeps a
group the same colour across every figure even when one panel is missing a
level, which positional scales silently get wrong.

Reserve grey for residual categories (`Other`, `Unclassified`, `Unknown`) in
every figure that has them.

### Saving

```r
save_fig <- function(plot, name, width = 180, height = 120, dpi = 300) {
  ggsave(file.path("figures", paste0(name, ".png")), plot,
         width = width, height = height, units = "mm", dpi = dpi)
  ggsave(file.path("figures", paste0(name, ".pdf")), plot,
         width = width, height = height, units = "mm", device = cairo_pdf)
  invisible(plot)
}
```

Always give explicit dimensions. Without them `ggsave` falls back to the current
device size, which differs between machines and is why figures come out with
14pt axis labels on a three-inch panel.

Useful sizes: 180 mm is roughly full width in a two-column journal, 85 mm is a
single column. Design at final print size. Rendering large and shrinking in the
manuscript is what produces unreadable 5pt axis text.

Save a vector copy alongside the PNG. Text stays selectable and scales cleanly.

### Multi-panel

```r
p_combined <- (p_alpha | p_nmds) +
  plot_layout(widths = c(1.4, 1), guides = "collect") +
  plot_annotation(tag_levels = "A")
```

`guides = "collect"` merges shared legends instead of drawing the same one
twice. `tag_levels` generates the panel letters rather than having you type them
into the plot.

---

## 13. Things that go wrong

A checklist of the failures that cost the most time, most of which fail quietly
rather than with an error.

**Sample IDs stop matching.** `read.table` prefixes IDs starting with a digit
with `X` unless you pass `check.names = FALSE`. `estimate_richness` converts
hyphens to dots. Both cause joins to silently produce `NA`, and phyloseq to
silently drop samples. Print `nsamples()` after building the object, and
`stopifnot` after every join.

**Filtering removes more than you meant.** Any `subset_taxa` comparison against
a string also removes every `NA` at that rank. Guard with `is.na()`.

**Bray-Curtis on raw counts.** Separates samples by sequencing depth. Transform
to proportions first.

**Alpha diversity on proportions.** Returns numbers, means nothing. Needs
integer counts.

**A significant PERMANOVA that's actually a dispersion difference.** Run
`betadisper` every time.

**Unstable NMDS between runs.** Set a seed. Report stress; above 0.20 the plot
isn't interpretable.

**FEAST changing your working directory.** Capture the root and restore after
every call.

**Figures rendered at the wrong size.** Pass explicit `width`, `height`, `units`
and `dpi` to `ggsave`.

**Missing install instructions.** phyloseq is on Bioconductor; FEAST, microeco
and pairwiseAdonis are on GitHub. `install.packages()` on any of them fails.
