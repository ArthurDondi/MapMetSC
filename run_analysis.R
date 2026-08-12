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

# --- Supplementary: AD region / marker images (customised params) ------------
# Standalone AD image generator (mirrors analysis/08_generate_images.Rmd for the
# main cohort). Reads spe_04_phenotyping_AD.rds, so it runs after AD phenotyping.
# Customise the image, zoom crop and marker list here (not ad_params).
ad_image_params <- list(
  input    = input_AD,
  output   = output,
  image_id = "20250623_19700101_25-XDD397_AD_001",
  markers  = c("DNA1_Ir193_mean", "ELAVL4_Yb174_mean", "Vimentin_Pt196_mean",  # DNA2 / ELAVL4 / VIM
               "Ki-67_Tm169_mean", "CXCR4_Tb159_mean", "CHGA_Dy164_mean",       # KI67 / CXCR4 / CHGA
               "GD2_Gd155_mean"),                                               # GD2
  marker_colours = c("DNA1" = "blue", "ELAVL4" = "green", "Vimentin" = "red",
                     "Ki-67" = "darkorange", "CXCR4" = "cyan", "CHGA" = "yellow",
                     "GD2" = "purple"),
  # Only Ki-67 is left out of the merge (per request); that keeps the merge at
  # cytomapper's 6-marker max. GD2 is purple (not orange) so where it overlaps
  # ELAVL4 the blend is pale, not the green -> yellow an orange channel gave.
  merge_exclude = c("Ki-67"),
  crop     = c(2000, 3000, 2000, 3000)   # zoom c(xmin, xmax, ymin, ymax); NULL = whole image only
)
render_step("analysis_AD/05_generate_images.Rmd", ad_image_params)

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
render_step("analysis/07_spatial_analysis.Rmd",     main_params)

# --- Supplementary: region / marker images (customised params) ---------------
# Standalone image-generation step. Unlike the steps above it does NOT reuse
# main_params: set the image, the zoom crop and the marker list here. Reads
# spe_07_spatial_analysis.rds (compartment + cellular-neighborhood annotations),
# so it runs after the spatial analysis. Renders both the whole image and the
# zoomed region.
image_params <- list(
  input    = input,
  output   = output,
  image_id = "20220926_20220809_16-006_TU_003",
  markers  = c("Vimentin_Pt196_mean", "CD3_Sm152_mean", "DNA1_Ir193_mean",  # Ir193 = DNA2
               "CD14_Lu175_mean", "GD2_Gd155_mean", "CD15_Bi209_mean"),
  # colour per marker (individual image + merge), keyed by short marker name
  marker_colours = c(CD3 = "red", Vimentin = "purple", CD14 = "cyan",
                     DNA1 = "blue", GD2 = "green", CD15 = "yellow"),
  # all six markers go into the merge (GD2 takes CD56's green slot; none excluded)
  merge_exclude = character(0),
  crop     = c(2000, 3000, 2000, 3000)   # zoom c(xmin, xmax, ymin, ymax); NULL = whole image only
)
render_step("analysis/08_generate_images.Rmd", image_params)

message("\nAll steps completed. Rendered HTML is in docs/; objects & figures in ", output)
