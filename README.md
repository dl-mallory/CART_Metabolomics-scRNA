# CART Metabolomics SingleCell RNAseq

Everything needed to rebuild panels **B–H** of the CAR-T metabolomics figure from
raw reads, and nothing else. 

NK/myeloid cells co-cultured with adipocytes vs fibroblasts, one BD Rhapsody WTA
run, two multiplexed mice.

```
FASTQ ──stage 0──▶ BD Rhapsody v2.2.1 ──▶ counts + sample tags
                                             │
                   stage 1 ─────────────────▶│ QC, scran, Harmony, cluster, annotate
                                             │   → nk_annotated.rds  (13,325 cells)
                                             │   → panels B, C, D
                   stage 2 ─────────────────▶│ FindMarkers → panels E, G
                   stage 3 ─────────────────▶│ GSEA        → panels F, H
```

## Run it

```sh
./run_all.sh              # stages 1-3, ~16 min (stage 3 runs GSEA at eps = 0)
```

Nothing to install first: the pinned `msigdbr 24.1.0` ships in `env/rlib/`, and
`.Rprofile` puts it ahead of the system library. `Rscript env/install.R` rebuilds
it from `env/msigdbr_24.1.0.tar.gz` if you need to.

Stage 0 is separate and not included here (it needs the mouse WTA reference and
~64 GB RAM). Its outputs are pre-staged in `data/bd_pipeline_out/`. For the
record, they were produced with the **BD Rhapsody Sequence Analysis Pipeline
v2.2.1** (CWL) against `RhapRef_Mouse_WTA_2023-02`, Sample Tag version `mm`.

**Push with the `git` CLI, not the browser.** The raw counts
(`data/bd_pipeline_out/.../matrix.mtx.gz`, 72 MB) and the largest processed matrix
(`results/matrices/expr_macrophages.csv.gz`, 37 MB) both exceed the 25 MB cap on
GitHub's web uploader. Both are well under git's 100 MB per-file hard limit, so a
normal `git push` takes them (it prints a size advisory above 50 MB).

## Layout

| path | contents |
|---|---|
| `data/fastq/` | raw reads (**filenames still to be supplied**) |
| `data/bd_pipeline_out/` | the two Rhapsody outputs the analysis actually consumes |
| `R/00_setup.R` | all paths, thresholds, colours, cluster→cell-type map, expected counts |
| `R/01_qc_annotate.R` | panels B, C, D + `nk_annotated.rds` |
| `R/02_de_volcanoes.R` | panels E, G |
| `R/03_gsea_published.R` | panels F, H — **reproduces the published panels, defects included** |
| `R/combine_matrices.R` | rebuilds the full Seurat object from the CSV exports |
| `results/tables/` | the data plotted in every panel, as CSV |
| `results/matrices/` | processed expression, one dense CSV per cell type |
| `checks/` | the archived term lists, used as regression fixtures |
| `env/rlib/` | pinned `msigdbr 24.1.0` |

## Data exports

Every number behind every panel is written as CSV, so the figures can be checked
without running R.

`results/tables/` — one row per plotted element:

| file | panel | one row per |
|---|---|---|
| `panelB_umap_celltype.csv` | B | cell: UMAP xy, cell type, the hex colour drawn |
| `panelC_marker_expression.csv.gz` | C | cell × marker: expression, and expression scaled to that gene's own max (which is what the colour encodes, so facets are not comparable) |
| `panelD_umap_sample.csv` | D | cell: UMAP xy, arm, the alpha used in each variant |
| `panelE_tcells_volcano.csv` | E | gene: position, colour class, whether labelled |
| `panelG_macrophages_volcano.csv` | G | gene: as above |
| `panelF_gsea_tcells.csv` | F | plotted gene set, in plot order |
| `panelH_gsea_macrophages.csv` | H | plotted gene set, in plot order |

`Class` in the volcano tables reconstructs EnhancedVolcano's four colour
categories. For panel G it reproduces the published SVG's circle census exactly:
7,801 `NS`, 839 `log2FC_only`, 57 `p_adj_only`, 109 `p_adj_and_log2FC`.

The GSEA panel tables carry `NES` even though the dot plot never encodes it —
that is the column that exposes panel F's three fibroblast-enriched rows.

Also written: `cells_metadata.csv` (per-cell annotation + UMAP coordinates),
`embeddings_{pca,harmony}.csv.gz`, the full `de_*.csv`, and the full
`gsea_*.csv.gz`.

### Per-cell-type matrices, and putting them back together

`results/matrices/expr_<celltype>.csv.gz` holds processed (scran log-normalised)
expression, genes × cells, one file per cell type, with `manifest.csv` recording
the dimensions of each. Values are rounded to 4 decimal places — which, on this
dataset, loses no non-zeros at all. 74 MB across 8 files, the largest 37 MB
(macrophages), all far below git's 100 MB per-file limit.

Every file carries **all 26,182 genes in identical row order**, so the files
concatenate on columns with no re-alignment. Filtering each file to its own
detected genes would save ~1% and would silently misalign anyone who `cbind`-ed
them, so it is not done. `combine_matrices.R` asserts the invariant on every file.

A cell type above `CELLS_PER_FILE` would be split into `_p01`, `_p02` … parts,
splitting on **cells, not genes**, so parts concatenate exactly like whole cell
types. Nothing splits at the current threshold; the export warns if any file
approaches the git limit.

```sh
Rscript R/combine_matrices.R            # -> results/objects/nk_from_csv.rds
Rscript R/combine_matrices.R --no-save  # rebuild and verify only
```

### Raw counts

The CSVs hold log-normalised expression only, and no `counts` layer is attached to
either object — the raw counts already sit, losslessly, in
`data/bd_pipeline_out/`. Duplicating 22.1 M integers into an `.rds` would add
~46 MB of numbers the repo already has.

They are recovered exactly, whenever needed, with one line:

```r
cts <- Read10X(data.dir = BD$mex)[["Gene Expression"]][genes, meta$Cell]
```

where `genes` is column 1 of any `expr_*.csv.gz` and `meta$Cell` is the `Cell`
column of `cells_metadata.csv`. Verified against the counts layer of the original
pre-slimming object: same dimensions, same dimnames, same sparsity pattern,
`max |diff| = 0`, **bit-identical**.

It rebuilds expression, metadata, cluster and cell-type factors, `Idents()`, and
the PCA / Harmony / UMAP reductions; re-checks the pipeline's own guards (13,325
cells, 16 clusters, exact per-cell-type and per-arm counts); and, when
`nk_annotated.rds` is present, diffs the expression against it. Observed
`max |diff| = 5.0e-05`, which is exactly the 4-dp rounding and nothing else.

### What the saved objects contain

Both are `data`-only: log-normalised expression, metadata, and the three
reductions. Neither carries a `counts` layer — **raw counts live in
`data/bd_pipeline_out/`** and are not reconstructible from these exports.

| | size | dropped |
|---|---|---|
| `nk_annotated.rds` | 58 MB | `scale.data` (dense 2000 × 13325), `RNA_nn`/`RNA_snn` graphs, `counts`, loadings, command log |
| `nk_from_csv.rds` | 54 MB | — (built data-only from the CSVs) |

Everything dropped is recomputable scratch from `ScaleData`, `FindNeighbors` and
`RunPCA`, and none of it is read by stages 2–3, which use only the `data` layer.
`FindMarkers` returns identical results either way (8,806 macrophage genes, same
ranking). Keeping it cost 82 MB.

Both objects are under the gitignored `results/objects/`; the committed record is
the CSVs.


## Environment

The one thing that must be pinned is **`msigdbr 24.1.0` (MSigDB v2024.1)**. The
figure was built 2025-06-27; `msigdbr 25.1.0` (MSigDB v2025.1) shipped 2025-07-03.
Between releases, gene sets grew:

```
                            v2024.1   v2025.1
INFLAMMATORY_RESPONSE          467       507   ← crosses maxGSSize = 500
CYTOKINE_PRODUCTION            489       504   ← crosses maxGSSize = 500
NUCLEIC_ACID_CATABOLIC_PROCESS 499       234
```

Under v2025.1 two of panel F's rows are excluded from the analysis outright.
`.Rprofile` puts `env/rlib` ahead of the system library; `03_gsea_published.R`
asserts the version before doing anything.

Everything else in the GSEA stack (clusterProfiler 4.12.6, DOSE 3.30.5,
fgsea 1.30.0, qvalue 2.36.0) predates the figure and needs no pin.


