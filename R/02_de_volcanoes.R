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
source("R/lib/export.R")

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

# ---- CSV exports: the DE result, and the points each volcano actually draws --
# ident.1 = Adipo_Cocultured, ident.2 = Fibro_Cocultured, so a positive
# avg_log2FC means higher in the adipocyte arm.
for (nm in c("macrophages", "tcells")) {
    de <- if (nm == "macrophages") macro else tcell
    write_table(data.frame(Gene = rownames(de), de, row.names = NULL),
                sprintf("de_%s.csv", nm), DIR$tables)
}

# The volcano source tables: every plotted point, its position, and the class
# that determines its colour. Labelled == the genes EnhancedVolcano annotates.
panel_source <- function(de, p_cutoff, fc_cutoff, labelled) {
    data.frame(
        Gene = rownames(de),
        avg_log2FC = de$avg_log2FC,
        p_val = de$p_val, p_val_adj = de$p_val_adj,
        neg_log10_p_adj = -log10(de$p_val_adj),
        pct.1 = de$pct.1, pct.2 = de$pct.2,
        Class = volcano_class(de, p_cutoff, fc_cutoff),
        Labelled = rownames(de) %in% labelled,
        row.names = NULL
    )
}

macro_labels <- rownames(macro)[macro$p_val_adj < 1e-100]

panelG <- panel_source(macro, 1e-50, 1,    macro_labels)
panelE <- panel_source(tcell, 1e-5,  0.75, rownames(tcell)[volcano_class(tcell, 1e-5, 0.75) == "p_adj_and_log2FC"])
write_table(panelG, "panelG_macrophages_volcano.csv", DIR$tables)
write_table(panelE, "panelE_tcells_volcano.csv",      DIR$tables)

# The class census is what was checked against the published SVG's circle colours.
stopifnot(nrow(panelG) == EXPECTED$de_genes[["Macrophages"]])
message("  panel G class census: ",
        paste(sprintf("%s=%d", names(table(panelG$Class)), table(panelG$Class)), collapse = " "))

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

# Only the most extreme macrophage genes are labelled (macro_labels, above);
# everything is still plotted.
ggsave(file.path(DIR$figures, "panelG_macrophages_volcano.svg"),
       volcano(macro, p_cutoff = 1e-50, fc_cutoff = 1, select_lab = macro_labels),
       width = 9, height = 9)

ggsave(file.path(DIR$figures, "panelE_tcells_volcano.svg"),
       volcano(tcell, p_cutoff = 1e-5, fc_cutoff = 0.75),
       width = 7, height = 7)

message(sprintf("stage 2 complete: %d macrophage genes, %d T cell genes plotted",
                nrow(macro), nrow(tcell)))
