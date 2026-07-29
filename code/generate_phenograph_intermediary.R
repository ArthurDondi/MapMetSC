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

## ---- lost cells: copy the retained cells' labels out of spe_02 -------------
# spe_02 stores pg_clusters_lostcells for the RETAINED cells only (the ghost
# cluster was excluded before saving). That is all we need: 02_QC_1.Rmd rebuilds
# the ghost cells at load time (colnames(spe) there is still the full post-QC1
# set, so any cell not in this file is a ghost, labelled `lostcells_ghost`).
# 02_QC_1.Rmd also verifies that ghost number is not itself a retained cluster,
# so nothing about the ghost number is needed here.
extract_lostcells <- function(spe02_path, out_name) {
  if (!file.exists(spe02_path)) {
    message("skip lost-cells: ", spe02_path, " not found"); return(invisible())
  }
  spe02 <- readRDS(spe02_path)
  if (!"pg_clusters_lostcells" %in% names(colData(spe02)))
    stop("pg_clusters_lostcells not found in ", spe02_path)
  v <- setNames(as.character(colData(spe02)$pg_clusters_lostcells), colnames(spe02))
  saveRDS(v, file.path(inter_dir, out_name))
  message(sprintf("wrote %s  (%d retained cells, %d clusters; ghosts rebuilt at load)",
                  out_name, length(v), length(unique(v))))
}

# --- AD pipeline ------------------------------------------------------------
extract_lostcells(file.path(output_dir, "spe_02_QC_1_AD.rds"),
                  "pg_clusters_lostcells_AD.rds")
extract_pheno(file.path(output_dir, "spe_04_phenotyping_AD.rds"),
              "pg_clusters_AD.rds")

# --- PT/BM pipeline ---------------------------------------------------------
extract_lostcells(file.path(output_dir, "spe_02_QC_1.rds"),
                  "pg_clusters_lostcells_PTBM.rds")
extract_pheno(file.path(output_dir, "spe_04_phenotyping.rds"),
              "pg_clusters_PTBM.rds")

message("\nDone. Intermediaries written to: ", inter_dir)
