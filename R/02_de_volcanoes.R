# Stage 2 -- differential expression, panels E and G.
# Inputs : results/objects/nk_annotated.rds
# Outputs: results/figures/published/panelE_tcells_volcano.svg
#          results/figures/published/panelG_macrophages_volcano.svg
#          results/objects/de_{macrophages,tcells}.csv
#
# Verified: with no pre-plot filtering, panelG regenerates byte-identically to the
# published macrophage_de.svg (8,806 points; 7,801 NS / 839 FC-only / 57 p-only /
# 109 both). Panel E matches on every count and label; only ggrepel's unseeded
# label placement differs.

suppressPackageStartupMessages({
    library(Seurat); library(tidyverse); library(EnhancedVolcano)
})
source("R/00_setup.R")

nk <- readRDS(file.path(DIR$objects, "nk_annotated.rds"))

# A volcano must show every tested gene. Do NOT filter on p_val_adj before
# plotting -- doing so removes the entire null cloud and leaves a plot where every
# point sits above the significance line by construction.
de_for <- function(cell_type) {
    FindMarkers(subset(nk, subset = CellType1 == cell_type),
                group.by = "SampleID", ident.1 = IDENT_1, ident.2 = IDENT_2)
}

macro <- de_for("Macrophages")
tcell <- de_for("T Cells")

stopifnot(
    nrow(macro) == EXPECTED$de_genes[["Macrophages"]],
    nrow(tcell) == EXPECTED$de_genes[["T Cells"]]
)

write.csv(macro, file.path(DIR$objects, "de_macrophages.csv"))
write.csv(tcell, file.path(DIR$objects, "de_tcells.csv"))

# y = p_val_adj is Seurat's Bonferroni correction over all 26,182 genes in the
# object, so the axis is -log10(adjusted p), not -log10(P). Labelled accordingly.
volcano <- function(de, p_cutoff, fc_cutoff, select_lab = NULL) {
    EnhancedVolcano(
        toptable = de, lab = rownames(de),
        x = "avg_log2FC", y = "p_val_adj",
        pCutoff = p_cutoff, FCcutoff = fc_cutoff,
        ylab = bquote(~-Log[10] ~ italic(P)[adj]),
        title = NULL, subtitle = NULL, caption = NULL,
        gridlines.major = FALSE, gridlines.minor = FALSE,
        legendPosition = "none", pointSize = 3, labSize = 6,
        drawConnectors = TRUE, arrowheads = FALSE, widthConnectors = 0.5,
        lengthConnectors = unit(1, "npc"),
        selectLab = select_lab
    )
}

# Only the most extreme macrophage genes are labelled; everything is still plotted.
macro_labels <- rownames(macro)[macro$p_val_adj < 1e-100]

ggsave(file.path(DIR$figures, "panelG_macrophages_volcano.svg"),
       volcano(macro, p_cutoff = 1e-50, fc_cutoff = 1, select_lab = macro_labels),
       width = 9, height = 9)

ggsave(file.path(DIR$figures, "panelE_tcells_volcano.svg"),
       volcano(tcell, p_cutoff = 1e-5, fc_cutoff = 0.75),
       width = 7, height = 7)

message(sprintf("stage 2 complete: %d macrophage genes, %d T cell genes plotted",
                nrow(macro), nrow(tcell)))
