# Export helpers: every number behind every panel, as CSV.
#
# Two kinds of output:
#   results/tables/    small, analysis-level tables (one row per cell / gene / set)
#   results/matrices/  processed expression, one dense CSV per cell type
#
# The matrices carry ALL genes in identical row order in every file, so they can
# be cbind()-ed without re-aligning. Filtering each file to its own detected
# genes would be smaller and would silently misalign anyone who concatenated
# them. See R/combine_matrices.R.

suppressPackageStartupMessages({library(data.table); library(Matrix)})

slugify <- function(x) gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(x)))

#' Write a data frame to results/tables/<name>.csv (gzipped if name ends .gz).
write_table <- function(df, name, dir) {
    path <- file.path(dir, name)
    data.table::fwrite(df, path)
    message(sprintf("    %-42s %6.2f MB", name, file.size(path) / 1048576))
    invisible(path)
}

#' Dense genes x cells CSV, gzipped, written in gene-chunks so the full dense
#' matrix is never materialised. Values are log-normalised expression, rounded.
write_expr_csv <- function(mat, path, digits = 4, chunk = 2000L) {
    csv <- sub("\\.gz$", "", path)
    if (file.exists(csv)) unlink(csv)
    starts <- seq(1L, nrow(mat), by = chunk)
    for (s in starts) {
        idx <- s:min(s + chunk - 1L, nrow(mat))
        blk <- as.matrix(mat[idx, , drop = FALSE])
        dt  <- data.table::data.table(gene = rownames(mat)[idx])
        dt  <- cbind(dt, data.table::as.data.table(round(blk, digits)))
        data.table::fwrite(dt, csv, append = (s != starts[1]), col.names = (s == starts[1]))
        rm(blk, dt); gc(verbose = FALSE)
    }
    if (file.exists(path)) unlink(path)
    ok <- system2("gzip", c("-9", "-f", shQuote(csv)), stdout = FALSE, stderr = FALSE)
    if (ok != 0L || !file.exists(path))
        stop("gzip failed for ", csv, " -- left uncompressed")
    invisible(path)
}

# Cells per file, and the size we refuse to exceed.
#
# This repo is pushed with the git CLI, whose hard per-file limit is 100 MB (it
# warns above 50 MB). GitHub's *web* uploader caps files at 25 MB, but that path
# cannot take data/bd_pipeline_out/matrix.mtx.gz (72 MB) either, so it is not a
# constraint we design around. The largest cell type, macrophages, is 6,585 cells
# -> 39 MB in one file: comfortably under the limit, so nothing splits today.
#
# The split mechanism is kept for a future, larger dataset. It splits on CELLS,
# never on genes, so every part keeps the full gene set in the same row order and
# parts concatenate exactly like whole cell types do.
CELLS_PER_FILE <- 100000L
MAX_FILE_MB    <- 90

#' Dense CSVs per cell type, split into parts, + a manifest. Returns the manifest.
export_celltype_matrices <- function(nk, dir, digits = 4) {
    dat <- SeuratObject::GetAssayData(nk, layer = "data")
    rows <- list()
    for (ct in levels(nk$CellType1)) {
        m <- dat[, nk$CellType1 == ct, drop = FALSE]
        starts <- seq(1L, ncol(m), by = CELLS_PER_FILE)
        for (p in seq_along(starts)) {
            idx <- starts[p]:min(starts[p] + CELLS_PER_FILE - 1L, ncol(m))
            mp  <- m[, idx, drop = FALSE]
            f <- file.path(dir, if (length(starts) == 1L)
                    sprintf("expr_%s.csv.gz", slugify(ct))
                else sprintf("expr_%s_p%02d.csv.gz", slugify(ct), p))
            write_expr_csv(mp, f, digits)
            rows[[length(rows) + 1L]] <- data.frame(
                CellType = ct, part = p, file = basename(f),
                genes = nrow(mp), cells = ncol(mp),
                nonzero = length(mp@x),
                pct_nonzero = round(100 * length(mp@x) / (as.numeric(nrow(mp)) * ncol(mp)), 2),
                file_mb = round(file.size(f) / 1048576, 2)
            )
            message(sprintf("    %-42s %6.2f MB  (%d x %d)",
                            basename(f), file.size(f) / 1048576, nrow(mp), ncol(mp)))
            rm(mp); gc(verbose = FALSE)
        }
    }
    manifest <- do.call(rbind, rows)
    stopifnot(sum(manifest$cells) == ncol(nk))
    # every file must carry the same genes, or combine_matrices.R is unsound
    stopifnot(length(unique(manifest$genes)) == 1L)
    # git hard-rejects files over 100 MB; stay well clear
    if (any(manifest$file_mb > MAX_FILE_MB))
        warning(sprintf("a matrix file exceeds %d MB (git rejects >100 MB); lower CELLS_PER_FILE",
                        MAX_FILE_MB), call. = FALSE)
    data.table::fwrite(manifest, file.path(dir, "manifest.csv"))
    manifest
}

#' EnhancedVolcano's four colour classes, recovered as a column so the volcano
#' source table says why each point is the colour it is.
volcano_class <- function(de, p_cutoff, fc_cutoff) {
    sig <- de$p_val_adj < p_cutoff
    fc  <- abs(de$avg_log2FC) > fc_cutoff
    ifelse(sig & fc, "p_adj_and_log2FC",
    ifelse(sig, "p_adj_only",
    ifelse(fc, "log2FC_only", "NS")))
}
