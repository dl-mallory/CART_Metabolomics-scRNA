# Rebuild the full processed matrix, and a complete Seurat object, from the
# per-cell-type CSV exports. This is the inverse of the export in 01_qc_annotate.R.
#
#   Rscript R/combine_matrices.R              # rebuild, verify, save
#   Rscript R/combine_matrices.R --no-save    # rebuild and verify only
#
# Outputs: results/objects/nk_from_csv.rds  (data-only, ~54 MB)
#
# Every expr_*.csv.gz carries all genes in identical row order, so the files
# concatenate on columns. That invariant is asserted, not assumed -- if a future
# export filters genes per cell type, this stops rather than silently misaligning.
#
# The CSVs hold log-normalised expression only. Raw counts are not derivable from
# them and are not attached -- they already sit, losslessly, in
# data/bd_pipeline_out/. See the note at the assay construction below.

suppressPackageStartupMessages({
    library(Seurat); library(data.table); library(Matrix)
})
source("R/00_setup.R")

save_obj <- !("--no-save" %in% commandArgs(trailingOnly = TRUE))

manifest <- fread(file.path(DIR$matrices, "manifest.csv"), header = TRUE)
meta     <- fread(file.path(DIR$tables, "cells_metadata.csv"), header = TRUE)
# Cell barcodes are BD cell indices -- all digits, so fread types them integer.
# Matrix dimnames are always character; compare like with like.
meta[, Cell := as.character(Cell)]
stopifnot(nrow(manifest) > 0, nrow(meta) == EXPECTED$cells_total)

# ---- read and concatenate ---------------------------------------------------
genes <- NULL; blocks <- list()
for (i in seq_len(nrow(manifest))) {
    f  <- file.path(DIR$matrices, manifest$file[i])
    # header = TRUE is required, not cosmetic: every column after `gene` is
    # numeric in both the header and the first data row, so fread's type-based
    # detection reads the header as data and invents V1..Vn names.
    dt <- fread(f, showProgress = FALSE, header = TRUE)
    g  <- dt[[1]]; dt[, 1 := NULL]

    if (is.null(genes)) genes <- g
    # the invariant that makes column-binding sound
    stopifnot(identical(g, genes))
    stopifnot(ncol(dt) == manifest$cells[i])

    blocks[[i]] <- Matrix(as.matrix(dt), sparse = TRUE)
    message(sprintf("  %-32s %5d cells", manifest$file[i], ncol(dt)))
    rm(dt); invisible(gc(verbose = FALSE))
}

expr <- do.call(cbind, blocks)
rownames(expr) <- genes
rm(blocks); invisible(gc(verbose = FALSE))

stopifnot(ncol(expr) == EXPECTED$cells_total,
          setequal(colnames(expr), meta$Cell))

# ---- restore cell order, metadata, reductions -------------------------------
expr <- expr[, meta$Cell, drop = FALSE]

# The export holds log-normalised expression, not counts. Build a data-only assay:
# CreateSeuratObject(counts = expr) would file log values under `counts`, storing
# the matrix twice and mislabelling one copy as raw.
#
# No `counts` layer is attached. Raw counts are already in the repo, losslessly,
# in data/bd_pipeline_out/ -- duplicating them here would add ~46 MB of the same
# numbers. To recover them exactly (verified bit-identical to the counts layer of
# the original object):
#
#   cts <- Read10X(data.dir = BD$mex)[["Gene Expression"]][genes, meta$Cell]
#
# `genes` is column 1 of any expr_*.csv.gz; `meta$Cell` is cells_metadata.csv.
assay <- SeuratObject::CreateAssay5Object(data = expr)
nk <- CreateSeuratObject(assay, assay = "RNA")

md <- as.data.frame(meta); rownames(md) <- md$Cell; md$Cell <- NULL
md$CellType1       <- factor(md$CellType1)
md$seurat_clusters <- factor(md$seurat_clusters,
                             levels = as.character(sort(as.integer(unique(md$seurat_clusters)))))
nk <- AddMetaData(nk, md)
Idents(nk) <- "CellType1"

for (red in c("pca", "harmony")) {
    f <- file.path(DIR$tables, sprintf("embeddings_%s.csv.gz", red))
    if (!file.exists(f)) next
    e <- as.data.frame(fread(f, header = TRUE))
    rownames(e) <- as.character(e$Cell); e$Cell <- NULL
    nk[[red]] <- CreateDimReducObject(embeddings = as.matrix(e[colnames(nk), ]),
                                      key = paste0(red, "_"), assay = "RNA")
}
umap <- as.matrix(md[colnames(nk), c("umap_1", "umap_2")])
nk[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "umap_", assay = "RNA")

# ---- verify against the expectations the pipeline asserts elsewhere ----------
stopifnot(
    ncol(nk) == EXPECTED$cells_total,
    nlevels(nk$seurat_clusters) == EXPECTED$clusters,
    identical(c(table(nk$CellType1)[names(EXPECTED$celltype_counts)]), EXPECTED$celltype_counts),
    identical(c(table(nk$SampleID)[names(EXPECTED$sample_counts)]),    EXPECTED$sample_counts)
)

# If the authoritative object is present, check the reconstruction against it.
orig_path <- file.path(DIR$objects, "nk_annotated.rds")
if (file.exists(orig_path)) {
    orig <- readRDS(orig_path)
    o <- GetAssayData(orig, layer = "data")[genes, colnames(nk)]
    d <- max(abs(o - GetAssayData(nk, layer = "data")))
    message(sprintf("\n  vs nk_annotated.rds: max |diff| in expression = %.2e  (export rounds to 4 dp)", d))
    message(sprintf("  cell types identical: %s   clusters identical: %s",
                    identical(as.character(nk$CellType1), as.character(orig$CellType1)),
                    identical(as.character(nk$seurat_clusters), as.character(orig$seurat_clusters))))
    stopifnot(d <= 1e-4)
    rm(orig, o); invisible(gc(verbose = FALSE))
} else {
    message("\n  nk_annotated.rds absent; skipped the against-original check")
}

message(sprintf("\n  rebuilt: %d genes x %d cells, %d cell types, reductions: %s",
                nrow(nk), ncol(nk), nlevels(nk$CellType1),
                paste(Reductions(nk), collapse = ", ")))

if (save_obj) {
    saveRDS(nk, file.path(DIR$objects, "nk_from_csv.rds"))
    message("  wrote ", file.path(DIR$objects, "nk_from_csv.rds"))
}
message("combine_matrices.R complete")
