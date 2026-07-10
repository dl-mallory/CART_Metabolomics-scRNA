# Stage 1 -- QC, normalisation, integration, clustering, annotation.
# Inputs : data/bd_pipeline_out/  (stage 0 output)
# Outputs: results/objects/nk_annotated.rds
#          results/figures/published/  panels B, C, D
#
# Verified: reproduces the published object bit-identically -- 13,325 cells,
# ARI 1.0 against its clustering, and PCA/Harmony/UMAP embeddings with
# max|diff| = 0.

suppressPackageStartupMessages({
    library(Seurat); library(scran); library(tidyverse)
    library(harmony); library(ggrastr); library(scattermore)
})
source("R/00_setup.R")
source("R/lib/export.R")

# ---- load counts and sample tags -------------------------------------------
counts <- Read10X(data.dir = BD$mex)[["Gene Expression"]]
nk <- CreateSeuratObject(counts = counts, assay = "RNA")

nk[["percent.mt"]] <- PercentageFeatureSet(nk, pattern = "^mt-")
nk[["percent.rb"]] <- PercentageFeatureSet(nk, pattern = "^Mrp|Rp[sl]")

sample_tags <- read.csv(BD$tags, skip = BD_TAG_HEADER_LINES) |>
    column_to_rownames("Cell_Index") |>
    mutate(SampleID = dplyr::coalesce(SAMPLE_MAP[Sample_Name], Sample_Name)) |>
    dplyr::select(SampleID)

nk <- AddMetaData(nk, sample_tags)

# ---- quality filtering ------------------------------------------------------
nk <- subset(nk, subset =
    nCount_RNA   < QC$max_counts &
    nFeature_RNA > QC$min_features &
    percent.mt   < QC$max_percent_mt &
    percent.rb   < QC$max_percent_rb &
    !(SampleID %in% c("Multiplet", "Undetermined")) & !is.na(SampleID))

# ---- scran pooled-size-factor normalisation ---------------------------------
sce <- as.SingleCellExperiment(nk)
sce <- computeSumFactors(sce, clusters = quickCluster(sce))
sce <- logNormCounts(sce)
nk  <- as.Seurat(sce)

# ---- remove NIH-3T3 fibroblast carry-over -----------------------------------
# FIBROBLAST_GENES is a character vector, so this produces FibroblastScore1..11,
# one module per gene. FibroblastScore1 is Col1a1 alone. This is what the
# published object was gated on and is preserved for exact reproduction.
nk <- AddModuleScore(nk, features = FIBROBLAST_GENES, name = "FibroblastScore")
nk <- subset(nk, subset = FibroblastScore1 < QC$fibroblast_cutoff)

stopifnot(ncol(nk) == EXPECTED$cells_total)

# Diagnostic only -- the correct 11-gene signature, computed AFTER the gate so it
# cannot influence which cells survive. Compare against FibroblastScore1.
nk$FibroblastScoreFull <- AddModuleScore(
    nk, features = list(FIBROBLAST_GENES), name = "FullSig")$FullSig1

# ---- dimensionality reduction, integration, clustering ----------------------
nk <- FindVariableFeatures(nk, selection.method = "vst",
                           nfeatures = QC$n_variable_features, verbose = FALSE)
nk <- ScaleData(nk, verbose = FALSE)
nk <- RunPCA(nk, npcs = QC$n_pcs, verbose = FALSE)
nk <- RunHarmony(nk, group.by.vars = "SampleID")

nk <- RunUMAP(nk, reduction = "harmony", dims = 1:QC$n_pcs, seed.use = QC$umap_seed)
nk <- FindNeighbors(nk, reduction = "harmony", dims = 1:QC$n_pcs)
nk <- FindClusters(nk, resolution = QC$cluster_resolution)

# ---- annotate ---------------------------------------------------------------
# unname(): the lookup carries cluster ids as names, which Seurat would otherwise
# match against cell barcodes.
stopifnot(all(levels(nk$seurat_clusters) %in% names(CLUSTER_ANNOTATION)))
nk$CellType1 <- as.factor(unname(CLUSTER_ANNOTATION[as.character(nk$seurat_clusters)]))
Idents(nk) <- "CellType1"

stopifnot(
    nlevels(nk$seurat_clusters) == EXPECTED$clusters,
    identical(c(table(nk$CellType1)[names(EXPECTED$celltype_counts)]), EXPECTED$celltype_counts),
    identical(c(table(nk$SampleID)[names(EXPECTED$sample_counts)]),   EXPECTED$sample_counts)
)
stopifnot(setequal(names(CELLTYPE_COLORS), levels(nk$CellType1)))

# Save a slimmed object: `data` + reductions + metadata. Dropped are scale.data
# (a DENSE 2000 x 13325 block), the RNA_nn/RNA_snn graphs, and the `counts` layer
# -- 140 MB -> 58 MB. All three are recomputable scratch, not data: counts are in
# data/bd_pipeline_out/, and stages 2-3 use only the `data` layer (FindMarkers
# returns identical results either way, verified). `nk` in memory keeps
# everything, so the panels below are unaffected.
saveRDS(DietSeurat(nk, layers = "data", dimreducs = c("pca", "harmony", "umap"),
                   graphs = NULL),
        file.path(DIR$objects, "nk_annotated.rds"))

# ---- CSV exports: the object, as tables ------------------------------------
# nk_annotated.rds is the authoritative object. These reconstitute it: the
# per-cell-type matrices carry the processed expression, cells_metadata.csv the
# per-cell annotation, and the embeddings the reductions. R/combine_matrices.R
# rebuilds a Seurat object from them and checks it against the original.
message("  exporting per-cell metadata and embeddings")
cells_metadata <- data.frame(
    Cell            = colnames(nk),
    SampleID        = as.character(nk$SampleID),
    seurat_clusters = as.character(nk$seurat_clusters),
    CellType1       = as.character(nk$CellType1),
    nCount_RNA      = nk$nCount_RNA,
    nFeature_RNA    = nk$nFeature_RNA,
    percent.mt      = round(nk$percent.mt, 4),
    percent.rb      = round(nk$percent.rb, 4),
    FibroblastScore1    = round(nk$FibroblastScore1, 6),
    FibroblastScoreFull = round(nk$FibroblastScoreFull, 6),
    round(Embeddings(nk, reduction = "umap"), 6),
    row.names = NULL
)
stopifnot(nrow(cells_metadata) == EXPECTED$cells_total)
write_table(cells_metadata, "cells_metadata.csv", DIR$tables)

for (red in c("pca", "harmony")) {
    e <- data.frame(Cell = colnames(nk), round(Embeddings(nk, reduction = red), 6),
                    row.names = NULL)
    write_table(e, sprintf("embeddings_%s.csv.gz", red), DIR$tables)
}

message("  exporting processed expression, one dense CSV per cell type")
manifest <- export_celltype_matrices(nk, DIR$matrices)

# ---- panels -----------------------------------------------------------------
# geom_jitter_rast() jitters without a seed; fix it so panels are pixel-stable.
set.seed(42)

umap <- data.frame(Embeddings(nk, reduction = "umap")) |>
    mutate(SampleID = nk$SampleID, CellType1 = nk$CellType1)

# Panel B and panel D use different point sizes and raster resolutions; keep them
# as parameters rather than baking in one panel's values.
base_umap <- function(mapping, colours, alpha_col = NULL,
                      pointsize = c(6, 5), raster.dpi = 450, scale = 0.75) {
    aes_full <- if (is.null(alpha_col)) ggplot2::aes(x = umap_1, y = umap_2)
                else ggplot2::aes(x = umap_1, y = umap_2, alpha = .data[[alpha_col]])
    ggplot(umap, aes_full) +
        geom_scattermore(pointsize = pointsize[1], colour = "darkgray") +
        geom_scattermore(pointsize = pointsize[2], colour = "white") +
        ggrastr::geom_jitter_rast(mapping, raster.dpi = raster.dpi, scale = scale) +
        scale_color_manual(values = colours) +
        theme_classic() +
        theme(panel.background = element_blank(),
              legend.key.height = unit(0.4, "cm"),
              legend.spacing.y  = unit(0.01, "cm"),
              legend.text = element_text(margin = margin(t = 0, b = -1))) +
        guides(color = guide_legend(byrow = TRUE))
}

# Panel B -- both variants derive from ONE plot object, so they cannot drift apart
panel_b <- base_umap(aes(color = CellType1), CELLTYPE_COLORS)

# Plotted data for panel B: one row per point, with the colour it is drawn in.
write_table(
    data.frame(Cell = rownames(umap), umap_1 = round(umap$umap_1, 6),
               umap_2 = round(umap$umap_2, 6),
               CellType1 = as.character(umap$CellType1),
               colour = unname(CELLTYPE_COLORS[as.character(umap$CellType1)])),
    "panelB_umap_celltype.csv", DIR$tables)
ggsave(file.path(DIR$figures, "panelB_umap_with_legend.svg"),
       panel_b + theme(legend.position = "right"), width = 5, height = 5)
ggsave(file.path(DIR$figures, "panelB_umap_without_legend.svg"),
       panel_b + theme(legend.position = "none"), width = 5, height = 5)

# Panel D -- one arm highlighted at a time
umap$AdipoAlpha <- if_else(umap$SampleID == "Adipo_Cocultured", 1, 0.3)
umap$FibroAlpha <- if_else(umap$SampleID == "Fibro_Cocultured", 1, 0.3)

# Plotted data for panel D: same points as B, keyed on arm, with each variant's
# alpha. The two SVGs differ only in which arm is opaque.
write_table(
    data.frame(Cell = rownames(umap), umap_1 = round(umap$umap_1, 6),
               umap_2 = round(umap$umap_2, 6),
               SampleID = as.character(umap$SampleID),
               AdipoAlpha = umap$AdipoAlpha, FibroAlpha = umap$FibroAlpha),
    "panelD_umap_sample.csv", DIR$tables)

for (arm in c("Adipo", "Fibro")) {
    p <- base_umap(aes(color = SampleID),
                   if (arm == "Adipo") ADIPO_HIGHLIGHT else FIBRO_HIGHLIGHT,
                   alpha_col = paste0(arm, "Alpha"),
                   pointsize = c(4, 3.5), raster.dpi = 600, scale = 0.4) +
        theme(legend.position = "none")
    ggsave(file.path(DIR$figures, sprintf("panelD_umap_sample_%s.svg", tolower(arm))),
           p, width = 4, height = 4)
}

# Panel C -- per-marker expression, scaled to each gene's own maximum
expr <- data.frame(t(GetAssayData(nk, layer = "data")[PANEL_C_GENES, ])) |>
    rownames_to_column("Cell")

panel_c_data <- data.frame(Embeddings(nk, reduction = "umap")) |>
    rownames_to_column("Cell") |>
    left_join(expr, by = "Cell") |>
    pivot_longer(all_of(PANEL_C_GENES), names_to = "Gene", values_to = "Expression") |>
    mutate(Gene = factor(Gene, levels = PANEL_C_GENES)) |>
    group_by(Gene) |>
    mutate(ExpressionScaled = Expression / max(Expression, na.rm = TRUE)) |>
    ungroup()

# Plotted data for panel C: one row per cell per marker gene. ExpressionScaled is
# what the colour encodes -- expression divided by that gene's own maximum, so
# colours are not comparable between facets.
write_table(panel_c_data |> mutate(across(where(is.numeric), \(x) round(x, 6))),
            "panelC_marker_expression.csv.gz", DIR$tables)

panel_c <- ggplot(panel_c_data, aes(x = umap_1, y = umap_2)) +
    geom_scattermore(pointsize = 4, colour = "black") +
    geom_scattermore(pointsize = 3.5, colour = "white") +
    ggrastr::geom_jitter_rast(aes(color = ExpressionScaled), raster.dpi = 600, scale = 0.4) +
    scale_color_gradientn(colors = c("lightyellow", "seagreen", "darkblue")) +
    theme_classic() +
    theme(panel.background = element_blank(), panel.border = element_blank(),
          legend.position = "none", strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 20),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.line = element_blank()) +
    facet_wrap(~Gene, ncol = 6)

ggsave(file.path(DIR$figures, "panelC_per_marker_umap.svg"), panel_c, width = 10, height = 6)

message("stage 1 complete: ", ncol(nk), " cells, ", nlevels(nk$seurat_clusters), " clusters")
