#!/usr/bin/env Rscript
# ============================================================================
# MapMetSC — master script
#
# Renders the R-based single-cell analysis pipeline (Lazic et al.) end to end,
# in order. Each step reads the object saved by the previous step from `output`
# and writes its own object and figures back to `output`.
#
# Reviewer setup
# --------------
#   1. Download and extract the MapMetIP-processed single-cell data from Zenodo
#      (see the README). The small annotation/metadata files (Supplementary
#      Tables, celltype_order.csv, protein2gene.csv, NB_Panel.csv) and the
#      public reference datasets (Fetahu/Jansky/Lee) are expected inside the
#      extracted Zenodo folders as well.
#   2. Point the four paths in the CONFIG block below at that data. The defaults
#      assume the Docker mount `-v <singlecelldata>:/mnt/data` (see README).
#   3. From the repository root run:   Rscript run_analysis.R
#      (or source("run_analysis.R") inside the RStudio session of the provided
#      Docker image).
#
# Only the four paths in the CONFIG block need to be set. Rendered HTML is
# written to docs/; intermediate objects (spe_*.rds) and figures to `output`.
#
# NOTE: to rebuild the versioned workflowr website instead of a one-off render,
#       use workflowr::wflow_publish() with the same file order.
# ============================================================================

# ---------------------------- CONFIG (edit me) ------------------------------
# Point these at your data. Either export the MAPMET_* environment variables
# (e.g. `export MAPMET_INPUT=/abs/path/to/MapMetIP_ProcessedDataset`), or just
# replace the fallback strings after the commas below. The same variables drive
# the per-Rmd `params:` defaults, so `Rscript run_analysis.R` and knitting a
# single .Rmd interactively both pick them up.
input    <- Sys.getenv("MAPMET_INPUT",  "/mnt/data/input_all")        # main cohort: the MapMetIP_ProcessedDataset folder (steinbock data + NB_Panel.csv; celltype_order.csv lives in R_intermediary)
input_AD <- input                                                    # adrenal-gland reference is the SAME dataset; analysis_AD/ filters for the AD samples
public   <- Sys.getenv("MAPMET_PUBLIC", "/mnt/data/public_datasets")  # public scRNA-seq references (Fetahu / Jansky / Lee) + protein2gene.csv
output   <- Sys.getenv("MAPMET_OUTPUT", "/mnt/data/output")           # writable directory for intermediate objects & figures
# ----------------------------------------------------------------------------

dir.create(output, showWarnings = FALSE, recursive = TRUE)

render_step <- function(rmd, params) {
  message("\n=====================================================================")
  message("Rendering: ", rmd)
  message("=====================================================================")
  rmarkdown::render(
    input      = rmd,
    output_dir = "docs",
    params     = params,
    envir      = new.env(),
    quiet      = FALSE
  )
}

main_params <- list(input = input,    output = output)
corr_params <- list(input = input,    output = output, public = public)
ad_params   <- list(input = input_AD, output = output)

# --- Adrenal-gland (AD) reference pipeline -----------------------------------
# Runs before the correlation step: produces spe_phenotyping_AD.rds, which is
# consumed by analysis/06.1_correlation_AD.Rmd.
render_step("analysis_AD/01_read_data.Rmd",   ad_params)
render_step("analysis_AD/02_QC_1.Rmd",        ad_params)
render_step("analysis_AD/03_QC_2.Rmd",        ad_params)
render_step("analysis_AD/04_phenotyping.Rmd", ad_params)

# --- Main single-cell pipeline ----------------------------------------------
# Shared colors live in analysis/00_MapMet_color_code.R, which each script
# sources automatically - it is a helper, not a rendered step.
render_step("analysis/01_read_data.Rmd",              main_params)   # incl. patient cohort overview
render_step("analysis/02_QC_1.Rmd",                   main_params)
render_step("analysis/03_QC_2.Rmd",                   main_params)
render_step("analysis/04_phenotyping.Rmd",            main_params)   # writes spe_04_phenotyping.rds
render_step("analysis/05.1_PT_DE_cellularcomm.Rmd",   main_params)   # writes spe_05.1_PT_DE_cellularcomm.rds
render_step("analysis/05.2_PT_DE_DA.Rmd",             main_params)   # reads  spe_05.1_PT_DE_cellularcomm.rds

# --- Validation by correlation with reference datasets (run last) ------------
render_step("analysis/06.1_correlation_AD.Rmd",       corr_params)   # needs AD + main phenotyping objects
render_step("analysis/06.2_correlation_Fetahu.Rmd",   corr_params)
render_step("analysis/06.3_correlation_Jansky.Rmd",   corr_params)
render_step("analysis/06.4_correlation_Lee.Rmd",      corr_params)

# --- Spatial analysis --------------------------------------------------------
# Uses imcRtools::testInteractions (the fix that once required a separate spatial
# image is in the imcRtools 1.10.0 / Bioconductor 3.19 this image pins) and the
# allcelltype_community column produced by 05.1, so it now runs in the same pass.
render_step("analysis/10_spatial_analysis.Rmd",     main_params)

message("\nAll steps completed. Rendered HTML is in docs/; objects & figures in ", output)
