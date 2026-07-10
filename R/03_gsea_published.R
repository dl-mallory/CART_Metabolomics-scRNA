# Stage 3 -- GSEA, panels F and H, reproduced AS PUBLISHED.
# Inputs : results/objects/nk_annotated.rds
# Outputs: results/figures/published/panel{F,H}_*.svg
#          results/objects/gsea_{tcells,macrophages}.csv
#
# READ THIS BEFORE USING THE OUTPUT.
#
# Panel F is filed as "upregulated in adipocyte co-culture" but was never filtered
# on NES. Its top two rows, RRNA METABOLIC PROCESS (NES = -2.05) and RRNA
# PROCESSING (NES = -2.04), are enriched in the FIBROBLAST arm. Panel H escaped
# only because GSEA_KEYWORDS happens to exclude every negatively-enriched set.
#
# This script exists to document and regenerate what was published. It is not a
# recommendation. To make panel F say what its caption says, filter on NES sign
# before selecting rows, and map NES to the fill so a wrong-direction row is
# visibly wrong rather than invisible. The warning() at the foot of this script
# names the offending rows on every run.

suppressPackageStartupMessages({
    library(Seurat); library(tidyverse); library(msigdbr)
    library(clusterProfiler); library(BiocParallel); library(scales)
})
source("R/00_setup.R")
source("R/lib/gsea_utils.R")

stopifnot(packageVersion("msigdbr") == "24.1.0")   # MSigDB v2024.1; see README

nk <- readRDS(file.path(DIR$objects, "nk_annotated.rds"))

gene_sets <- msigdbr(species = "Mus musculus",
                     collection = "C5", subcollection = "GO:BP") |>
    dplyr::rename(gene_id = gene_symbol) |>
    dplyr::select(gs_name, gene_id)

gsea_for <- function(cell_type) {
    ranked <- rank_genes(subset(nk, subset = CellType1 == cell_type))
    run_gsea(ranked, gene_sets)
}

tcell_gsea <- gsea_for("T Cells")
macro_gsea <- gsea_for("Macrophages")

write.csv(tcell_gsea, file.path(DIR$objects, "gsea_tcells.csv"), row.names = FALSE)
write.csv(macro_gsea, file.path(DIR$objects, "gsea_macrophages.csv"), row.names = FALSE)

# Panel F: no keyword filter. This is the omission that lets negatively-enriched
# gene sets into a panel captioned "upregulated".
panel_f <- format_gsea_panel(tcell_gsea)

# Panel H: keyword filter applied.
panel_h <- format_gsea_panel(macro_gsea[grepl(GSEA_KEYWORDS, macro_gsea$ID), ])

message("reproduction against the archived panels:")
check_panel(panel_f, "tcell_gsea_upreg_terms.txt", "panel F (T cells)")
check_panel(panel_h, "macro_gsea_upreg_terms.txt", "panel H (macrophages)")

for (nm in c("F", "H")) {
    df <- if (nm == "F") panel_f else panel_h
    cell <- if (nm == "F") "tcells" else "macrophages"
    for (leg in c(TRUE, FALSE)) {
        suffix <- if (leg) "with_legend" else "without_legend"
        ggsave(file.path(DIR$figures, sprintf("panel%s_gsea_%s_%s.svg", nm, cell, suffix)),
               gsea_dotplot(df, legend = leg), width = 8, height = 5)
    }
}

# Make the direction visible in the numbers even though the plot hides it.
report <- function(df, gsea, label) {
    ids <- paste0("GOBP_", gsub(" ", "_", as.character(df$Description)))
    nes <- gsea$NES[match(ids, gsea$ID)]
    neg <- which(!is.na(nes) & nes < 0)
    if (length(neg))
        warning(sprintf("%s: %d of %d rows have NEGATIVE NES (enriched in %s): %s",
                        label, length(neg), nrow(df), IDENT_2,
                        paste(as.character(df$Description)[neg], collapse = "; ")),
                call. = FALSE)
}
report(panel_f, tcell_gsea, "panel F")
report(panel_h, macro_gsea, "panel H")

message("stage 3 complete")
