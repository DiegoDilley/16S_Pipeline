# Amplicon Analysis in R

A walkthrough of what happens to 16S/18S amplicon data after denoising: from a
raw feature table through diversity analysis, source tracking, functional
prediction, and publication figures.

**Start here: [GUIDE.md](GUIDE.md)**

The guide explains each step, the decision you're making at that point, and the
code that implements it. The `R/` folder holds the same code as standalone
scripts if you'd rather copy from those.

This is a reference, not a reproducible pipeline. No data is included, and
filenames and group names are placeholders.

## Covered

- Building and filtering a phyloseq object from QIIME2 exports
- Rarefaction, and when it applies
- Alpha diversity: metrics, test selection, multiple-testing correction
- Beta diversity: distances, NMDS/PCoA, PERMANOVA, dispersion checks
- Taxonomic composition: agglomeration, top-N plus Other, stacked bars
- Group-specific taxa and subset characterisation
- Source tracking with FEAST
- Functional prediction with FAPROTAX
- Consistent theming and predictable figure export

## Dependencies

Most packages come from CRAN. `phyloseq` is on Bioconductor, and three
dependencies install from GitHub:

```r
BiocManager::install("phyloseq")
remotes::install_github("cozygene/FEAST")
remotes::install_github("ChiLiubio/microeco")
remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
```

Full list in [`R/00_dependencies.R`](R/00_dependencies.R).
