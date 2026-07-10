# Package installation. Run once, from the repo root:  Rscript env/install.R
#
# The ONLY version that must be pinned is msigdbr. The published figure was built
# against MSigDB v2024.1 (msigdbr 24.1.0, CRAN 2025-05-13). msigdbr 25.1.0 ships
# MSigDB v2025.1, under which gene set membership changed enough that two sets
# (INFLAMMATORY_RESPONSE 467 -> 507, CYTOKINE_PRODUCTION 489 -> 504) cross the
# maxGSSize = 500 default and drop out of the analysis entirely.

lib <- "env/rlib"
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran <- c("Seurat", "SeuratObject", "tidyverse", "harmony", "ggrastr",
          "scattermore", "svglite", "reshape2", "ggpubr", "scales", "BiocManager")
bioc <- c("scran", "scuttle", "SingleCellExperiment", "EnhancedVolcano",
          "clusterProfiler", "DOSE", "fgsea", "qvalue", "BiocParallel")

for (p in cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
for (p in bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)

# pinned
if (!requireNamespace("msigdbr", quietly = TRUE) ||
    packageVersion("msigdbr") != "24.1.0") {
    tarball <- "env/msigdbr_24.1.0.tar.gz"
    if (!file.exists(tarball))
        download.file("https://cran.r-project.org/src/contrib/Archive/msigdbr/msigdbr_24.1.0.tar.gz",
                      tarball, mode = "wb")
    install.packages(tarball, repos = NULL, type = "source", lib = lib)
}

stopifnot(packageVersion("msigdbr", lib.loc = lib) == "24.1.0")
writeLines(capture.output(sessionInfo()), "env/session_info.txt")
message("environment ready; msigdbr ", packageVersion("msigdbr", lib.loc = lib))
