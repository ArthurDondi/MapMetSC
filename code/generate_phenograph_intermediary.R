#!/usr/bin/env Rscript
# =============================================================================
# Extract frozen phenograph outputs from the reference SPE objects.
#
# Phenograph (Annoy kNN + igraph Louvain) is not bit-reproducible across
# machines, so we copy the kept-cell list (from spe_02) and the k=30 cluster
# labels (from spe_04), keyed by cell id, for the notebooks to load instead of
# re-clustering. Run once in the environment that produced the published
# objects, then ship <output_dir>/R_intermediary/ to params$input/R_intermediary/.
#
# Produces: kept_cells_{AD,PTBM}.rds, pg_clusters_{AD,PTBM}.rds
# =============================================================================

suppressPackageStartupMessages(library(SpatialExperiment))

output_dir <- Sys.getenv("MAPMET_OUTPUT", "/mnt/data/output")  # holds spe_*.rds
inter_dir  <- file.path(output_dir, "R_intermediary")
dir.create(inter_dir, showWarnings = FALSE, recursive = TRUE)

## ---- phenotyping: copy the k=30 labels out of the 04 object ----------------
# Use the first spe path that exists; prefer the pg_clusters_k30 column, else
# pg_clusters. For PT/BM this is the clustering before GMM subclustering.
extract_pheno <- function(spe_paths, out_name,
                          cols = c("pg_clusters_k30", "pg_clusters")) {
  spe_path <- spe_paths[file.exists(spe_paths)][1]
  if (is.na(spe_path)) {
    message("skip phenotyping: none found of ", paste(spe_paths, collapse = ", "))
    return(invisible())
  }
  spe <- readRDS(spe_path)
  col <- cols[cols %in% names(colData(spe))][1]
  if (is.na(col))
    stop("none of (", paste(cols, collapse = ", "), ") found in ", spe_path,
         " (pg_* columns present: ",
         paste(grep("pg_clusters", names(colData(spe)), value = TRUE),
               collapse = ", "), ")")
  v <- setNames(as.character(colData(spe)[[col]]), colnames(spe))
  saveRDS(v, file.path(inter_dir, out_name))
  message(sprintf("wrote %s  (%d cells, %d clusters; from %s col %s)",
                  out_name, length(v), length(unique(v)), basename(spe_path), col))
}

## ---- lost cells: the kept-cell list is just spe_02's cells -----------------
# 02_QC_1.Rmd loads this list to reproduce the lost-cell exclusion.
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
# AD's final clustering is k=15, stored in pg_clusters.
extract_pheno(file.path(output_dir, c("spe_final_clustering_AD.rds",
                                      "spe_04_phenotyping_AD.rds")),
              "pg_clusters_AD.rds",
              cols = c("pg_clusters", "pg_clusters_k15"))

# --- PT/BM pipeline ---------------------------------------------------------
extract_kept(file.path(output_dir, "spe_02_QC_1.rds"), "kept_cells_PTBM.rds")
extract_pheno(file.path(output_dir, c("spe_final_clustering.rds",
                                      "spe_04_phenotyping.rds")),
              "pg_clusters_PTBM.rds")

message("\nDone. Intermediaries written to: ", inter_dir)
