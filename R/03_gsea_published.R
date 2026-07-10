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
source("R/lib/export.R")

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

# ---- CSV exports: full GSEA result, and the 18 rows each panel draws ---------
for (nm in c("tcells", "macrophages")) {
    g <- if (nm == "tcells") tcell_gsea else macro_gsea
    write_table(g, sprintf("gsea_%s.csv.gz", nm), DIR$tables)
}

# Rows are written in plot order (bottom of the y-axis first, as the factor
# levels were set). PercentSet is the x position, Count the point size, qvalue
# the fill. NES is carried even though the panel never plots it -- it is the
# column that reveals panel F's mixed directions.
# df carries the untouched ID column; do not reconstruct it from Description,
# which format_gsea_panel() rewrites (it strips MEDIATED BY ANTIMICROBIAL PEPTIDE).
panel_rows <- function(df, gsea) {
    i <- match(df$ID, gsea$ID)
    stopifnot(!anyNA(i))
    data.frame(
        PlotOrder   = seq_len(nrow(df)),
        Description = as.character(df$Description),
        ID          = df$ID,
        PercentSet  = round(df$PercentSet, 6),
        Count       = df$Count,
        setSize     = df$setSize,
        qvalue      = df$qvalue,
        pvalue      = gsea$pvalue[i],
        p.adjust    = gsea$p.adjust[i],
        NES         = round(gsea$NES[i], 6),
        enrichmentScore = round(gsea$enrichmentScore[i], 6),
        core_enrichment = gsea$core_enrichment[i],
        row.names = NULL
    )
}
write_table(panel_rows(panel_f, tcell_gsea), "panelF_gsea_tcells.csv",      DIR$tables)
write_table(panel_rows(panel_h, macro_gsea), "panelH_gsea_macrophages.csv", DIR$tables)

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
    # match on the carried ID, not on a Description that format_gsea_panel rewrote
    nes <- gsea$NES[match(df$ID, gsea$ID)]
    stopifnot(!anyNA(nes))
    neg <- which(nes < 0)
    if (length(neg))
        warning(sprintf("%s: %d of %d rows have NEGATIVE NES (enriched in %s): %s",
                        label, length(neg), nrow(df), IDENT_2,
                        paste(as.character(df$Description)[neg], collapse = "; ")),
                call. = FALSE)
}
report(panel_f, tcell_gsea, "panel F")
report(panel_h, macro_gsea, "panel H")

message("stage 3 complete")
