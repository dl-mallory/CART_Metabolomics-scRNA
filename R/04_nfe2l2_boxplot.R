# Stage 4 -- Figure S12, Nfe2l2 expression by cell subset.
# Inputs : results/objects/nk_annotated.rds
# Outputs: results/figures/published/figS12_nfe2l2_by_celltype.svg
#          results/tables/figS12_nfe2l2_expression.csv.gz  (every plotted point)
#          results/tables/figS12_nfe2l2_pvalues.csv         (the 7 bracket p-values)
#
# Nfe2l2 expression is compared between the adipocyte- and fibroblast-cocultured
# arms within each cell subset. Zero-valued cells are dropped before testing and
# plotting; the bracket p-value is a Wilcoxon rank-sum test on the non-zero cells.

suppressPackageStartupMessages({
    library(Seurat); library(tidyverse); library(ggpubr); library(rstatix)
})
source("R/00_setup.R")
source("R/lib/export.R")

nk <- readRDS(file.path(DIR$objects, "nk_annotated.rds"))

# The figure shows 7 subsets; Mast Cells (34 cells) is not among them.
FIG_S12_CELLTYPES <- c("Dendritic Cells", "Macrophages", "Monocytes",
                       "Myeloid Progenitors", "Neutrophils", "NK Cells", "T Cells")

# Adipo left/salmon, Fibro right/turquoise -- ggplot's default two-colour fill.
FIG_S12_FILL <- c(Adipo_Cocultured = "#F8766D", Fibro_Cocultured = "#00BFC4")

# ---- data: Nfe2l2 per cell, with arm and subset -----------------------------
expr <- FetchData(nk, vars = c("SampleID", "CellType1", "Nfe2l2")) |>
    rownames_to_column("Cell") |>
    filter(CellType1 %in% FIG_S12_CELLTYPES) |>
    mutate(
        CellType1 = factor(CellType1, levels = FIG_S12_CELLTYPES),
        SampleID  = factor(SampleID, levels = c("Adipo_Cocultured", "Fibro_Cocultured"))
    )

# Zero-valued cells are dropped before testing and plotting.
expr_nz <- expr |> filter(Nfe2l2 != 0)

# ---- statistics: Wilcoxon rank-sum p per subset -----------------------------
# Computed once, here, so the number exported and the number drawn on the panel
# are guaranteed identical (stat_pvalue_manual, below).
stats <- expr_nz |>
    group_by(CellType1) |>
    wilcox_test(Nfe2l2 ~ SampleID) |>
    ungroup() |>
    adjust_pvalue(method = "BH") |>
    arrange(CellType1)

# ---- CSV exports ------------------------------------------------------------
write_table(
    expr |> mutate(Nfe2l2 = round(Nfe2l2, 6),
                   nonzero = Nfe2l2 != 0,
                   fill = unname(FIG_S12_FILL[as.character(SampleID)])),
    "figS12_nfe2l2_expression.csv.gz", DIR$tables)

pval_out <- stats |>
    transmute(CellType1, group1, group2,
              n_adipo = n1, n_fibro = n2,
              statistic, p = p, p_adj_BH = p.adj)
write_table(pval_out, "figS12_nfe2l2_pvalues.csv", DIR$tables)

# ---- figure -----------------------------------------------------------------
# Jitter is seeded so the panel is pixel-stable across runs.
set.seed(42)

brackets <- stats |>
    mutate(p.label = ifelse(p < 1e-3, formatC(p, format = "e", digits = 2),
                            formatC(p, format = "f", digits = 3)),
           y.position = max(expr_nz$Nfe2l2) * 1.05)

fig <- ggplot(expr_nz, aes(x = SampleID, y = Nfe2l2, fill = SampleID)) +
    geom_boxplot(show.legend = FALSE, outlier.shape = NA) +
    geom_violin(show.legend = FALSE, alpha = 0.3) +
    geom_jitter(show.legend = FALSE, width = 0.2, size = 0.6, alpha = 0.5) +
    stat_pvalue_manual(brackets, label = "p.label", tip.length = 0.01) +
    scale_fill_manual(values = FIG_S12_FILL) +
    facet_wrap(~CellType1, nrow = 1, strip.position = "bottom") +
    labs(x = NULL, y = "Nfe2l2 Expression") +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(file.path(DIR$figures, "figS12_nfe2l2_by_celltype.svg"),
       fig, width = 16, height = 8)

message(sprintf("stage 4 complete: Figure S12, %d subsets, %d non-zero cells plotted",
                nlevels(expr_nz$CellType1), nrow(expr_nz)))
