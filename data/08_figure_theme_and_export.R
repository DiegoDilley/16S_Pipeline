# ---------------------------------------------------------------------------
# 08 — Theming and export
#
# Input : outputs/plot_*.rds
# Output: figures/*.png, figures/*.pdf
#
# This script does no analysis. Earlier scripts build plot objects and save
# them unthemed; this one applies a single theme and writes every file.
#
# The reason to separate them: if each analysis script saves its own figures,
# the plotting code drifts apart, the same panel ends up rendered twice under
# two different numbers, and changing a font means editing seven files. One
# export script keeps figure numbering in one place — which is also the only
# place it can be kept in sync with the manuscript.
# ---------------------------------------------------------------------------

library(ggplot2)
library(patchwork)

dir.create("figures", showWarnings = FALSE)

# --- 1. One theme ----------------------------------------------------------
theme_pub <- function(base_size = 10, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey92"),
      panel.border     = element_rect(linewidth = 0.4, colour = "grey30"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(face = "bold", size = rel(0.9)),
      axis.text        = element_text(colour = "grey20"),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.key.size  = unit(0.4, "cm"),
      legend.title     = element_text(face = "bold", size = rel(0.9)),
      plot.title       = element_text(face = "bold", size = rel(1.1))
    )
}

theme_set(theme_pub())

# --- 2. One palette --------------------------------------------------------
# Fixed, named vectors rather than positional ones. A named palette keeps a
# group the same colour across every figure even when a panel is missing a
# level — which positional scales silently get wrong.
pal_group <- c(group_a = "#4E79A7", group_b = "#E15759")

# Categorical palette for taxa. Beyond ~12 categories, colour stops carrying
# information; lump to "Other" instead of extending the palette.
pal_taxa <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948",
  "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#86BCB6", "#D37295"
)

scale_fill_group  <- function(...) scale_fill_manual(values = pal_group, ...)
scale_col_group   <- function(...) scale_colour_manual(values = pal_group, ...)

# Grey for the residual categories, in every figure that has them.
fill_taxa <- function(levels_in_plot) {
  n_named <- length(setdiff(levels_in_plot, c("Other", "Unclassified", "Unknown")))
  cols <- setNames(pal_taxa[seq_len(n_named)],
                   setdiff(levels_in_plot, c("Other", "Unclassified", "Unknown")))
  c(cols, Unclassified = "#BFBFBF", Other = "#8C8C8C", Unknown = "#666666")
}

# --- 3. One save function --------------------------------------------------
# Explicit dimensions in mm, and a dpi that survives print. Defaults are the
# reason figures come out with 14pt axis labels on a 3-inch panel: ggsave
# falls back to the current device size, which differs between machines.
save_fig <- function(plot, name, width = 180, height = 120, dpi = 300) {
  # 180 mm ~ full width in a two-column journal; 85 mm ~ single column.
  ggsave(file.path("figures", paste0(name, ".png")), plot,
         width = width, height = height, units = "mm", dpi = dpi)
  # Vector copy for submission. Text stays selectable and scales cleanly.
  ggsave(file.path("figures", paste0(name, ".pdf")), plot,
         width = width, height = height, units = "mm", device = cairo_pdf)
  invisible(plot)
}

# --- 4. Render -------------------------------------------------------------
# Figure numbering lives here and nowhere else.

p_alpha <- readRDS("outputs/plot_alpha.rds") + scale_fill_group() +
  theme(legend.position = "none")     # x-axis already labels the groups
save_fig(p_alpha, "Figure_1_alpha_diversity", width = 180, height = 80)

p_nmds <- readRDS("outputs/plot_nmds.rds") + scale_col_group()
save_fig(p_nmds, "Figure_2_nmds", width = 120, height = 110)

p_comp <- readRDS("outputs/plot_composition.rds")
p_comp <- p_comp +
  scale_fill_manual(values = fill_taxa(levels(p_comp$data$taxon_lumped)))
save_fig(p_comp, "Figure_3_composition", width = 200, height = 120)

p_sources <- readRDS("outputs/plot_sources.rds")
save_fig(p_sources, "Figure_4_source_tracking", width = 180, height = 110)

p_func <- readRDS("outputs/plot_functions.rds") + scale_fill_group()
save_fig(p_func, "Figure_5_functional_groups", width = 180, height = 140)

# --- 5. Multi-panel --------------------------------------------------------
# patchwork composes saved plot objects. plot_layout(guides = "collect")
# merges shared legends instead of drawing the same one twice.
p_combined <- (p_alpha | p_nmds) +
  plot_layout(widths = c(1.4, 1), guides = "collect") +
  plot_annotation(tag_levels = "A")   # panel letters, generated not typed

save_fig(p_combined, "Figure_6_diversity_combined", width = 200, height = 90)

# --- 6. A note on sizing ---------------------------------------------------
# Design at final print size. Making a figure large and shrinking it in the
# manuscript is what produces unreadable 5pt axis text. If the panel will be
# 85 mm wide in print, render it at 85 mm and check legibility there.
