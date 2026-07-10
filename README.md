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
./run_all.sh              # stages 1-4, ~16 min (stage 3 runs GSEA at eps = 0)
```

Nothing to install first: the pinned `msigdbr 24.1.0` ships in `env/rlib/`, and
`.Rprofile` puts it ahead of the system library. `Rscript env/install.R` rebuilds
it from `env/msigdbr_24.1.0.tar.gz` if you need to.

Stage 0 is separate (needs the mouse WTA reference and ~64 GB RAM); its outputs
are pre-staged in `data/bd_pipeline_out/`. See
`pipeline/00_bd_rhapsody/README.md`.

`data/bd_pipeline_out/.../matrix.mtx.gz` is 72 MB. That is under git's 100 MB
limit but over the 25 MB cap on GitHub's web uploader, so push this repo with
the `git` CLI rather than by dragging it into the browser.

## Layout

| path | contents |
|---|---|
| `pipeline/00_bd_rhapsody/` | CWL inputs + provenance of the published Rhapsody run |
| `data/fastq/` | raw reads (**filenames still to be supplied**) |
| `data/bd_pipeline_out/` | the two Rhapsody outputs the analysis actually consumes |
| `R/00_setup.R` | all paths, thresholds, colours, cluster→cell-type map, expected counts |
| `R/01_qc_annotate.R` | panels B, C, D + `nk_annotated.rds` |
| `R/02_de_volcanoes.R` | panels E, G |
| `R/03_gsea_published.R` | panels F, H — **reproduces the published panels, defects included** |
| `checks/` | the archived term lists, used as regression fixtures |
| `env/rlib/` | pinned `msigdbr 24.1.0` |


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


