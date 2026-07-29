#!/usr/bin/env Rscript
# =============================================================================
# Extract the frozen phenograph outputs from the DOCKER (reference) SPE objects.
# NOTHING is re-clustered here: we copy out the kept-cell list (02) and the k=30
# phenograph labels (04) that were already produced in docker, keyed by the
# unique cell id (sample_id_ObjectNumber), so they are immune to the Annoy /
# igraph::cluster_louvain non-reproducibility across machines (the same caveat
# uwot documents for calculateUMAP).
#
# WHY freeze: Rphenoannoy (Annoy approximate-NN kNN + igraph Louvain) is not
# bit-reproducible across igraph versions / compilers / CPUs. 02_QC_1.Rmd loads
# the kept-cell list to reproduce the lost-cell exclusion; 04_phenotyping.Rmd
# loads the k=30 labels instead of re-clustering.
#
# HOW TO USE (run once in the docker image that produced the published objects):
#   1. Point `output_dir` at the folder holding your spe_*.rds (params$output),
#      or export MAPMET_OUTPUT.
#   2. Rscript code/generate_phenograph_intermediary.R
#   3. Ship the files written to <output_dir>/R_intermediary/ : commit them to
#      data/R_intermediary/ (they are small) or upload to Zenodo, then place
#      them in params$input/R_intermediary/ on any machine.
#
# Files produced:
#   kept_cells_AD.rds  / _PTBM.rds   (character vector of kept cell ids, from
#                                     the columns of the reference spe_02)
#   pg_clusters_AD.rds / _PTBM.rds   (named vector, cell id -> k=30 phenograph
#                                     cluster, from spe_04's pg_clusters_k30)
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

## ---- lost cells: the kept-cell list is just spe_02's cells -----------------
# 02_QC_1.Rmd now excludes lost cells by loading this list and keeping those
# cells, instead of re-identifying the (non-reproducible) lost-cells cluster.
# The kept-cell list is simply the columns of the reference spe_02.
extract_kept <- function(spe02_path, out_name) {
  if (!file.exists(spe02_path)) {
    message("skip kept-cells: ", spe02_path, " not found"); return(invisible())
  }
  spe02 <- readRDS(spe02_path)
  saveRDS(colnames(spe02), file.path(inter_dir, out_name))
  message(sprintf("wrote %s  (%d kept cells)", out_name, ncol(spe02)))
}

# --- AD pipeline ------------------------------------------------------------
extract_kept(file.path(output_dir, "spe_02_QC_1_AD.rds"), "kept_cells_AD.rds")
extract_pheno(file.path(output_dir, "spe_04_phenotyping_AD.rds"), "pg_clusters_AD.rds")

# --- PT/BM pipeline ---------------------------------------------------------
extract_kept(file.path(output_dir, "spe_02_QC_1.rds"), "kept_cells_PTBM.rds")
extract_pheno(file.path(output_dir, "spe_04_phenotyping.rds"), "pg_clusters_PTBM.rds")

message("\nDone. Intermediaries written to: ", inter_dir)
