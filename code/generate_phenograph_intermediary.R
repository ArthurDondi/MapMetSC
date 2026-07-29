#!/usr/bin/env Rscript
# =============================================================================
# Extract the frozen phenograph clusterings from the DOCKER (reference) SPE
# objects. NOTHING is re-clustered here: the labels already computed in docker
# (pg_clusters_lostcells and pg_clusters_k30) are copied out and keyed by the
# unique cell id (sample_id_ObjectNumber), so they are immune to the Annoy /
# igraph::cluster_louvain non-reproducibility across machines (the same caveat
# uwot documents for calculateUMAP).
#
# WHY freeze: Rphenoannoy (Annoy approximate-NN kNN + igraph Louvain) is not
# bit-reproducible across igraph versions / compilers / CPUs. 02_QC_1.Rmd and
# 04_phenotyping.Rmd load these files when present instead of recomputing.
#
# HOW TO USE (run once in the docker image that produced the published objects):
#   1. Point `output_dir` at the folder holding your spe_*.rds (params$output),
#      or export MAPMET_OUTPUT.
#   2. Rscript code/generate_phenograph_intermediary.R
#   3. Ship the files written to <output_dir>/R_intermediary/ : commit them to
#      data/R_intermediary/ (they are small) or upload to Zenodo, then place
#      them in params$input/R_intermediary/ on any machine.
#
# Files produced (one named character vector each, names = unique cell id):
#   pg_clusters_lostcells_AD.rds / _PTBM.rds   (from 02_QC_1: k=45 lost cells)
#   pg_clusters_AD.rds           / _PTBM.rds   (from 04: k=30 phenograph)
#
# NOTE the clustering is frozen as-is from these objects. If you regenerated the
# objects after changing an upstream step (e.g. the QC area threshold), extract
# from the regenerated objects so the frozen labels match the new cell set.
# =============================================================================

suppressPackageStartupMessages(library(SpatialExperiment))

output_dir <- Sys.getenv("MAPMET_OUTPUT", "/mnt/data/output")  # holds spe_*.rds
inter_dir  <- file.path(output_dir, "R_intermediary")
dir.create(inter_dir, showWarnings = FALSE, recursive = TRUE)

## ---- phenotyping: copy the k=30 labels out of the 04 object ----------------
extract_pheno <- function(spe_path, out_name, col = "pg_clusters_k30") {
  if (!file.exists(spe_path)) {
    message("skip phenotyping: ", spe_path, " not found"); return(invisible())
  }
  spe <- readRDS(spe_path)
  if (!col %in% names(colData(spe)))
    stop(col, " not found in ", spe_path, " (pg_* columns present: ",
         paste(grep("pg_clusters", names(colData(spe)), value = TRUE),
               collapse = ", "), "). ",
         "For PT/BM use the base k=30 column, not the post-subclustering ",
         "`pg_clusters`.")
  v <- setNames(as.character(colData(spe)[[col]]), colnames(spe))
  saveRDS(v, file.path(inter_dir, out_name))
  message(sprintf("wrote %s  (%d cells, %d clusters, from column %s)",
                  out_name, length(v), length(unique(v)), col))
}

## ---- lost cells: reconstruct full post-QC1 labels (retained + ghosts) ------
# spe_02 stores pg_clusters_lostcells only for the RETAINED cells (the ghost
# cluster was excluded before saving). The 02 load path needs a label for every
# post-QC1 cell, so we recompute the deterministic QC-step-1 cell set from
# spe_01 (no clustering) and label the missing cells - the ghosts - with
# `ghost_label`. That is exactly what the original clustering assigned them, so
# the reconstructed vector reproduces the docker `pg_clusters_lostcells`.
extract_lostcells <- function(spe01_path, spe02_path, out_name, ghost_label) {
  if (!file.exists(spe01_path) || !file.exists(spe02_path)) {
    message("skip lost-cells: ", spe01_path, " / ", spe02_path, " not found")
    return(invisible())
  }
  spe01 <- readRDS(spe01_path)
  im  <- assay(spe01, "counts")
  m80 <- im[grepl("mean-80", rownames(im)), ]
  imc <- m80[!grepl("IF2_GD2|IF3_CD45|DAPI", rownames(m80)), ]
  nz  <- as.logical(colSums(imc) != 0)
  post_qc1 <- colnames(spe01)[nz]

  spe02 <- readRDS(spe02_path)
  if (!"pg_clusters_lostcells" %in% names(colData(spe02)))
    stop("pg_clusters_lostcells not found in ", spe02_path)
  retained <- colnames(spe02)
  if (!all(retained %in% post_qc1))
    stop("Retained cells in ", spe02_path, " are not a subset of the ",
         "recomputed QC-step-1 set from ", spe01_path,
         "; are these the matching AD/PTBM objects?")

  lab <- setNames(rep(as.character(ghost_label), length(post_qc1)), post_qc1)
  lab[retained] <- as.character(colData(spe02)$pg_clusters_lostcells)
  saveRDS(lab, file.path(inter_dir, out_name))
  message(sprintf("wrote %s  (%d post-QC1 cells = %d retained + %d ghost [= %s])",
                  out_name, length(lab), length(retained),
                  length(post_qc1) - length(retained), ghost_label))
}

# --- AD pipeline ------------------------------------------------------------
extract_lostcells(file.path(output_dir, "spe_01_read_data_AD.rds"),
                  file.path(output_dir, "spe_02_QC_1_AD.rds"),
                  "pg_clusters_lostcells_AD.rds", ghost_label = "17")
extract_pheno(file.path(output_dir, "spe_04_phenotyping_AD.rds"),
              "pg_clusters_AD.rds")

# --- PT/BM pipeline ---------------------------------------------------------
extract_lostcells(file.path(output_dir, "spe_01_read_data.rds"),
                  file.path(output_dir, "spe_02_QC_1.rds"),
                  "pg_clusters_lostcells_PTBM.rds", ghost_label = "26")
extract_pheno(file.path(output_dir, "spe_04_phenotyping.rds"),
              "pg_clusters_PTBM.rds")

message("\nDone. Intermediaries written to: ", inter_dir)
