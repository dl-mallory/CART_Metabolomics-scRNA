# GSEA helpers.
#
# The published panels were produced by applying format_gsea_melt() straight to
# the GSEA result. That function selects the top 20 gene sets BY Q-VALUE -- not by
# NES -- and then keeps the 18 with the smallest PercentSet. Direction is never
# used, and NES is never plotted, so a set enriched in the opposite condition is
# indistinguishable on the page. See README, "Known defects".

# Keyword filter over gene set IDs. Note what this actually does: no
# negatively-enriched set in either comparison matches it (0/11 in T cells, 0/20
# in macrophages), so it silently doubles as a direction filter.
GSEA_KEYWORDS <- paste0(
    "(LEUKOCYTE|LYMPHOCYTE|MYELOID|GRANULOCYTE|NEUTROPHIL|MONOCYTE|MACROPHAGE|",
    "DENDRITIC|CYTOKINE|INTERLEUKIN|CHEMOKINE|INFLAMMATORY|INFLAMMATION|",
    "HEMATOPOIETIC|ERYTHROID|STROMAL|PHAGOCYTE|IMMUNE|T_CELL|B_CELL|NK_CELL|",
    "ADIPOCYTE|PROGENITOR|MIGRATION|APOPTOSIS|RECRUITMENT|FATTY|LIPID|INTERFERON|",
    "PHOSPHOLIPID|TRANSFORMING_GROWTH_FACTOR|ICOSANOID|P38MAPK|PROSTANOID|",
    "EPITHELIAL|MESENCHYMAL|FIBROBLAST|ADIPOCYTE)"
)

#' Rank genes for GSEA by log2 fold change over the full tested gene set.
rank_genes <- function(seurat_obj, ident.1 = IDENT_1, ident.2 = IDENT_2) {
    de <- Seurat::FindMarkers(seurat_obj, group.by = "SampleID",
                              ident.1 = ident.1, ident.2 = ident.2)
    sort(setNames(de$avg_log2FC, rownames(de)), decreasing = TRUE)
}

#' Run GSEA deterministically, and precisely enough that the ranking is stable.
#'
#' Two independent requirements, often confused:
#'
#' 1. DETERMINISM. seed = TRUE makes DOSE call set.seed(123) internally;
#'    SerialParam fixes the BiocParallel RNG stream. Together these make a rerun
#'    bit-identical to the last run. They do NOT make it correct.
#'
#' 2. PRECISION. format_gsea_panel() cuts the top 20 by q-value. Around that cut
#'    the macrophage q-values are ~1e-4 apart. At fgsea's default tail precision
#'    their Monte Carlo sd is ~8e-5, so the cut lands on a different gene set
#'    from run to run -- a seed only freezes WHICH wrong answer you get.
#'    eps = 0 and sampleSize = 1001 raise the precision of fgseaMultilevel's tail
#'    estimate until the ranking converges: the q-value CV on the deciding gene
#'    set (GOBP_AMEBOIDAL_TYPE_CELL_MIGRATION) falls 38% -> 9%, it lands at rank
#'    18 in every replicate, and the top-20 selection stops moving (Jaccard 1.000
#'    across independent unseeded runs). With this, panel H reproduces 18/18.
#'
#' nPermSimple governs the simple-permutation pass that normalises ES into NES,
#' not the tail p-values. Raising it to 1e5 stabilises NES but leaves the q-value
#' CV at 39% and panel H at 1/4. It is left at its default deliberately.
#'
#' Cost: ~3 min for both cell types, against ~15 s at default precision.
run_gsea <- function(gene_list, term2gene, n_perm_simple = 1000) {
    BiocParallel::register(BiocParallel::SerialParam(RNGseed = 123))
    res <- clusterProfiler::GSEA(
        geneList      = gene_list,
        TERM2GENE     = term2gene,
        pvalueCutoff  = 0.05,
        pAdjustMethod = "BH",
        seed          = TRUE,
        verbose       = FALSE,
        nPermSimple   = n_perm_simple,
        eps           = 0,
        sampleSize    = 1001
    )@result
    res$Count <- vapply(strsplit(res$core_enrichment, "/"), length, integer(1))
    res
}

#' Select and label rows exactly as the published panels did.
#'
#' Faithful reproduction of the original format_gsea_melt(). Two behaviours are
#' preserved on purpose and are both defects: selection is by q-value regardless
#' of NES sign, and the final slice_head() discards the two gene sets with the
#' LARGEST leading-edge fraction.
format_gsea_panel <- function(df, n_by_qvalue = 20, n_final = 18) {
    df |>
        dplyr::arrange(qvalue) |>
        dplyr::slice_head(n = n_by_qvalue) |>
        dplyr::mutate(PercentSet = Count / setSize, LogQ = -log10(qvalue)) |>
        dplyr::arrange(PercentSet) |>
        dplyr::mutate(
            Description = gsub("GOBP_", "", Description),
            Description = gsub("_", " ", Description),
            Description = trimws(gsub("MEDIATED BY ANTIMICROBIAL PEPTIDE", "", Description))
        ) |>
        dplyr::mutate(Description = factor(Description, levels = Description)) |>
        dplyr::slice_head(n = n_final) |>
        droplevels()
}

#' Format a number as an expression like 6.94 x 10^-8, for colour bar labels.
#'
#' The original hard-coded these as typed strings while computing `breaks` from
#' the data, so the published colour bars print q-values that are not the
#' q-values of the plotted points. Derived here instead.
sci_labels <- function(x) {
    txt <- vapply(x, function(v) {
        e <- floor(log10(v)); m <- signif(v / 10^e, 3)
        sprintf("%s %%*%% 10^{%d}", m, e)
    }, character(1))
    parse(text = txt)   # expression vector; a list of calls breaks grid's plotmath
}

#' The published GSEA dot plot: x = leading-edge fraction, size = leading-edge
#' size, fill = q-value. NES is deliberately absent, as in the original.
gsea_dotplot <- function(df, legend = TRUE) {
    q_breaks <- c(min(df$qvalue), stats::median(df$qvalue), max(df$qvalue))
    c_breaks <- unique(round(c(min(df$Count), stats::median(df$Count), max(df$Count))))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = PercentSet, y = Description, fill = qvalue)) +
        ggplot2::geom_point(ggplot2::aes(size = Count), colour = "black", pch = 21) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(colour = "black", size = 12),
            axis.text.x = ggplot2::element_text(colour = "black", size = 10),
            plot.margin = ggplot2::margin(l = 10, r = 10)
        ) +
        ggplot2::scale_fill_gradient(
            low = "#edb081", high = "#aa3a6f",
            breaks = q_breaks, limits = range(df$qvalue),
            labels = sci_labels(q_breaks),
            guide = ggplot2::guide_colorbar(ticks = TRUE)
        ) +
        ggplot2::scale_x_continuous(labels = scales::percent, expand = c(0.1, 0)) +
        ggplot2::scale_size_continuous(range = c(2, 8), breaks = c_breaks) +
        ggplot2::labs(x = "% of Gene Set", y = NULL) +
        ggplot2::guides(size = ggplot2::guide_legend(order = 1))

    if (!legend) p <- p + ggplot2::theme(legend.position = "none")
    p
}

#' Compare a reconstructed panel against the archived term list.
check_panel <- function(df, fixture, label) {
    expected <- trimws(readLines(file.path(DIR$checks, fixture)))
    got      <- trimws(as.character(df$Description))
    overlap  <- length(intersect(got, expected))
    message(sprintf("  %-28s %2d/%d terms match published; order exact: %s",
                    label, overlap, length(expected),
                    identical(got, rev(expected))))
    invisible(overlap)
}
