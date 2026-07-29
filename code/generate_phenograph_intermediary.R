#!/usr/bin/env Rscript
# =============================================================================
# Generate the frozen phenograph-clustering intermediaries (the "reference"
# clusterings) in the DOCKER environment.
#
# WHY: Rphenoannoy performs Annoy (approximate nearest-neighbour, C++) kNN and
# then igraph::cluster_louvain community detection. Neither is bit-reproducible
# across igraph versions, compilers or CPUs, so the SAME seed + SAME input can
# yield DIFFERENT clusters on a different machine/container (exactly the caveat
# the uwot/UMAP note describes for `calculateUMAP`). To make results identical
# everywhere we compute the clustering ONCE in the reference (docker) image and
# ship the labels; every other run loads them instead of recomputing.
#
# HOW TO USE (run inside the docker image that produced the published results):
#   1. Set `output_dir` below to the folder that holds your spe_*.rds objects
#      (i.e. the pipeline's params$output), or export MAPMET_OUTPUT.
#   2. Rscript code/generate_phenograph_intermediary.R
#   3. Ship the four files written to <output_dir>/R_intermediary/ : commit them
#      to data/R_intermediary/ (they are small) or upload to Zenodo, then place
#      them next to the input data (params$input/R_intermediary/) on any machine.
#      02_QC_1.Rmd and 04_phenotyping.Rmd then load them automatically.
#
# The matrix construction and seeds below are copied verbatim from
# analysis(_AD)/02_QC_1.Rmd (lost cells) and analysis(_AD)/04_phenotyping.Rmd
# (phenotyping), so the output is identical to a fresh pipeline run.
# =============================================================================

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(Rphenoannoy)
})

# Folder holding spe_*.rds (== the pipeline's params$output). Edit if needed.
output_dir <- Sys.getenv("MAPMET_OUTPUT", "/mnt/data/output")
inter_dir  <- file.path(output_dir, "R_intermediary")
dir.create(inter_dir, showWarnings = FALSE, recursive = TRUE)

# Marker exclusion pattern used for the phenotyping matrix (04_phenotyping.Rmd).
pheno_exclude <- paste0(
  "PARP|IDO|DAPI|IF2_GD2|IF3_CD45_CD56|DNA|H4K12Ac|H3K9Ac|",
  "CXCR2|PNMT|Fibro|FOXP3"
)

## ---- lost-cells clustering (as in 02_QC_1.Rmd) -----------------------------
gen_lostcells <- function(spe_path, out_name, seed = 230619, k = 45) {
  if (!file.exists(spe_path)) {
    message("skip lost-cells: ", spe_path, " not found"); return(invisible())
  }
  spe <- readRDS(spe_path)

  intensity_mat <- assay(spe, "counts")
  colnames(intensity_mat) <- spe$sample_id
  m80 <- intensity_mat[grepl("mean-80", rownames(intensity_mat)), ]
  imc <- m80[!grepl("IF2_GD2|IF3_CD45|DAPI", rownames(m80)), ]
  non_zeros <- as.logical(colSums(imc) != 0)
  spe <- spe[, non_zeros]
  mat <- m80[, non_zeros]
  mat <- mat[!grepl("IF2_GD2|IF3_CD45", rownames(mat)), ]
  mat <- t(mat)

  set.seed(seed)
  out <- Rphenoannoy(mat, k = k)
  labels <- setNames(as.character(membership(out[[2]])), colnames(spe))

  saveRDS(labels, file.path(inter_dir, out_name))
  message(sprintf("wrote %s  (%d cells, %d clusters)",
                  out_name, length(labels), length(unique(labels))))
}

## ---- phenotyping clustering (as in 04_phenotyping.Rmd) ----------------------
gen_phenotyping <- function(spe_path, out_name, seed = 221228,
                            k_values = c(15, 30, 45, 60)) {
  if (!file.exists(spe_path)) {
    message("skip phenotyping: ", spe_path, " not found"); return(invisible())
  }
  spe <- readRDS(spe_path)

  intensity_mat <- assay(spe, "counts")
  colnames(intensity_mat) <- colnames(spe)
  intensity_mat <- intensity_mat[!grepl("mean-80", rownames(intensity_mat)), ]
  hq <- intensity_mat[!grepl(pheno_exclude, rownames(intensity_mat)), ]
  non_zeros <- as.logical(colSums(hq) != 0)
  hq  <- hq[, non_zeros]
  spe <- spe[, non_zeros]
  mat <- hq

  df <- data.frame(row.names = colnames(spe))
  set.seed(seed)
  for (k in k_values) {
    out <- Rphenoannoy(t(mat), k = k)
    df[[paste0("pg_clusters_k", k)]] <- as.character(membership(out[[2]]))
  }

  saveRDS(df, file.path(inter_dir, out_name))
  message(sprintf("wrote %s  (%d cells, k = %s)",
                  out_name, nrow(df), paste(k_values, collapse = ",")))
}

# --- AD pipeline ------------------------------------------------------------
gen_lostcells(  file.path(output_dir, "spe_01_read_data_AD.rds"),
                "pg_clusters_lostcells_AD.rds")
gen_phenotyping(file.path(output_dir, "spe_03_QC_2_AD.rds"),
                "pg_clusters_phenotyping_AD.rds")

# --- PT/BM pipeline ---------------------------------------------------------
gen_lostcells(  file.path(output_dir, "spe_01_read_data.rds"),
                "pg_clusters_lostcells_PTBM.rds")
gen_phenotyping(file.path(output_dir, "spe_03_QC_2.rds"),
                "pg_clusters_phenotyping_PTBM.rds")

message("\nDone. Intermediaries written to: ", inter_dir)
