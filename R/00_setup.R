# Shared configuration. Sourced by every stage.

# msigdbr 24.1.0 (MSigDB v2024.1) is pinned in env/rlib and MUST take precedence.
# The published figure was built on v2024.1; v2025.1 changes gene set membership
# enough to alter which terms are selected. See README, "Environment".
.libPaths(c(normalizePath("env/rlib", mustWork = FALSE), .libPaths()))

DIR <- list(
    bd       = "data/bd_pipeline_out",
    objects  = "results/objects",
    figures  = "results/figures/published",
    checks   = "checks"
)

BD <- list(
    mex  = file.path(DIR$bd, "_1_read3-read4-test_RSEC_MolsPerCell_MEX"),
    tags = file.path(DIR$bd, "_1_read3-read4-test_Sample_Tag_Calls.csv")
)

# BD writes 7 lines of pipeline metadata above the Sample_Tag_Calls header
BD_TAG_HEADER_LINES <- 7

# Sample tag -> condition. Tag names come from the BD run: "1-mouse1; 2-mouse2".
SAMPLE_MAP <- c(mouse1 = "Fibro_Cocultured", mouse2 = "Adipo_Cocultured")

# Differential expression is always Adipocyte vs Fibroblast co-culture, so a
# positive avg_log2FC means enriched in the ADIPOCYTE arm.
IDENT_1 <- "Adipo_Cocultured"
IDENT_2 <- "Fibro_Cocultured"

QC <- list(
    max_counts   = 18000,
    min_features = 500,
    max_percent_mt = 15,
    max_percent_rb = 15,
    fibroblast_cutoff = 0.25,
    n_variable_features = 2000,
    n_pcs = 30,
    cluster_resolution = 0.7,
    umap_seed = 500
)

# Passed to AddModuleScore() as a bare character vector, which makes Seurat build
# ONE MODULE PER GENE (FibroblastScore1..11), not one module for the signature.
# FibroblastScore1 is therefore the module score of Col1a1 alone. The published
# object was gated on that. Preserved deliberately -- see README, "Known defects".
FIBROBLAST_GENES <- c("Col1a1", "Col1a2", "Fn1", "Dcn", "Thy1", "Pdgfra",
                      "Acta2", "Vim", "Tagln", "Serpine1", "Cxcl12")

# Cluster -> cell type at resolution 0.7. Recovered from the published object;
# replaying it reproduces CellType1 for all 13,325 cells with zero disagreements.
CLUSTER_ANNOTATION <- c(
    "0"  = "Macrophages",         "1"  = "Monocytes",
    "2"  = "Macrophages",         "3"  = "Macrophages",
    "4"  = "Macrophages",         "5"  = "Monocytes",
    "6"  = "Myeloid Progenitors", "7"  = "T Cells",
    "8"  = "NK Cells",            "9"  = "T Cells",
    "10" = "Myeloid Progenitors", "11" = "T Cells",
    "12" = "T Cells",             "13" = "Neutrophils",
    "14" = "Dendritic Cells",     "15" = "Mast Cells"
)

CELLTYPE_COLORS <- c(
    "Dendritic Cells"     = "#d3bcd4",
    "Macrophages"         = "#f08592",
    "Mast Cells"          = "#afa200",
    "Monocytes"           = "#9897df",
    "Myeloid Progenitors" = "#004583",
    "Neutrophils"         = "#c38200",
    "NK Cells"            = "#8a001b",
    "T Cells"             = "#f6d100"
)

SAMPLE_COLORS <- c(Adipo_Cocultured = "#055118", Fibro_Cocultured = "#3888a4")
ADIPO_HIGHLIGHT <- c(Adipo_Cocultured = "#055118", Fibro_Cocultured = "darkgray")
FIBRO_HIGHLIGHT <- c(Adipo_Cocultured = "darkgray", Fibro_Cocultured = "#3888a4")

# Panel C marker genes, in plotted order
PANEL_C_GENES <- c("Fscn1", "Cd3e", "Ncr1", "Gzmd", "Mpo", "Elane",
                   "S100a8", "S100a9", "Fn1", "Cd14", "Ccl8", "C1qc")

# Cell counts the pipeline must reproduce. Guards, not decoration: an upstream
# change that silently alters the cell set will stop the run here.
EXPECTED <- list(
    cells_total = 13325L,
    clusters    = 16L,
    celltype_counts = c(
        "Dendritic Cells" = 114L, "Macrophages" = 6585L, "Mast Cells" = 34L,
        "Monocytes" = 3167L, "Myeloid Progenitors" = 1029L, "Neutrophils" = 140L,
        "NK Cells" = 668L, "T Cells" = 1588L
    ),
    sample_counts = c(Adipo_Cocultured = 5854L, Fibro_Cocultured = 7471L),
    de_genes = c(Macrophages = 8806L, `T Cells` = 9877L)
)

dir.create(DIR$objects, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR$figures, recursive = TRUE, showWarnings = FALSE)
