#!/usr/bin/env Rscript

# ============================================================================
# MapMetSC - central color definitions
#
# Shared color palettes for the single-cell (IMC/IF) analysis, collected here so
# that every script uses the same colors. Source this file near the top of a
# script (path-robust across working directories) and then assign the vector(s)
# you need, e.g.
#
# the loop variable name below is deliberately obscure so it never clobbers a
# plot object named `p` in the caller's environment:
#   for (mapmet_cc_path in c("analysis/00_MapMet_color_code.R", "00_MapMet_color_code.R",
#                            "../analysis/00_MapMet_color_code.R"))
#     if (file.exists(mapmet_cc_path)) { source(mapmet_cc_path); break }
#   col_celltype <- COL_CELLTYPE
#
# All colors are base R rgb()/hex/named colors, so no extra packages are needed.
# ============================================================================

# ---- Main cohort: cell types (04_phenotyping.Rmd) --------------------------
COL_CELLTYPE <- c(
  ## Tumor
  "Ki67+ Tu"                  = rgb(128, 128,   0, maxColorValue = 255),  # olive
  "CD44+ Tu"                  = rgb( 50, 205,  50, maxColorValue = 255),  # bright lime green
  "GD2neg Tu"                 = rgb(189, 183, 107, maxColorValue = 255),  # khaki
  "early neuroblast-like Tu"  = rgb(152, 251, 152, maxColorValue = 255),  # pale green
  "adrenergic marker-low Tu"  = rgb( 34, 100,  60, maxColorValue = 255),  # dark green
  "neuroblast-like Tu"        = rgb(  0, 145, 150, maxColorValue = 255),  # teal / blue-green

  ## Proliferation
  "Ki67+ cells"               = rgb(255, 200,   0, maxColorValue = 255),
  "CXCR4+ cells"              = rgb(255, 245, 160, maxColorValue = 255),

  ## Progenitor
  "CD34+ HSPCs"               = rgb(255,  95,   0, maxColorValue = 255),

  ## T cells
  "CD8+ T cells"              = rgb(220,  20,  60, maxColorValue = 255),  # crimson
  "CD4+ T cells"              = rgb(255, 127,  80, maxColorValue = 255),  # coral
  "CD3+ T cells"              = rgb(139,   0,   0, maxColorValue = 255),  # dark red
  "CD8+ S100B+ T cells"       = rgb(230,  20, 210, maxColorValue = 255),  # fuchsia
  "CD8+ GZMB+ T cells"        = rgb(255, 182, 193, maxColorValue = 255),  # light pink
  "dense T cell region"       = rgb(199,  21, 133, maxColorValue = 255),  # magenta

  ## Myeloid
  "granulocytes"              = rgb(  8,  48, 107, maxColorValue = 255),
  "CD14+ MO"                  = rgb( 33, 102, 172, maxColorValue = 255),
  "myeloid cells"             = rgb( 64, 180, 170, maxColorValue = 255),
  "GZMB+ DC/NK cells"         = rgb(107, 174, 214, maxColorValue = 255),

  ## B cells
  "B cells"                   = rgb(136,  86, 167, maxColorValue = 255),

  ## Stromal
  "Schwann cells"             = rgb(166,  97,  26, maxColorValue = 255),
  "fibroblast/endothel"       = rgb(101,  67,  33, maxColorValue = 255),

  ## Other
  "other"                     = rgb(120, 120, 120, maxColorValue = 255)
)

# ---- Main cohort: metaclusters (04_phenotyping.Rmd) ------------------------
COL_METACLUSTER <- c(
  "tumor"               = "#238B45",
  "T cells"             = "#E34A33",
  "granulocytes"        = "#08306B",
  "myeloid/NK"          = "#41B6C4",
  "B cells"             = "#88419D",
  "mesenchymal"         = "#8C510A",
  "proliferating cells" = "#FFD700",
  "progenitors"         = "#F59E0B",
  "other"               = "#969696"
)

# ---- Adrenal-gland (AD) cohort: cell types (analysis_AD/04_phenotyping.Rmd) -
COL_CELLTYPE_AD <- c(
  ## Neural crest lineage (green tones)
  "Schwann cell progenitors" = rgb(166,  97,  26, maxColorValue = 255),  # brown (matches "Schwann cells" in the main cohort)
  "early neuroblasts"        = rgb(152, 251, 152, maxColorValue = 255),  # pale green
  "GD2+ neuroblasts"         = rgb( 50, 205,  50, maxColorValue = 255),  # bright lime green
  "CD15+ chromaffin cells"   = "#ffec17",  # bright yellow
  "CD15- chromaffin cells"   = "#e6c200",  # darker gold-yellow (distinct from CD15+)

  ## Mesenchymal (brown tones)
  "fibroblasts"             = rgb(210, 140,  80, maxColorValue = 255),  # tan (distinct from the Schwann brown)
  "LUM+ mesenchymal cells"  = rgb(140,  81,  10, maxColorValue = 255),  # dark brown
  "CD44+ mesenchymal cells" = rgb(101,  67,  33, maxColorValue = 255),  # deep brown

  ## Immune
  "immune cells"            = "#6aa9d2",  # blue

  ## Endothelial
  "endothelial cells"       = rgb(255,  95,   0, maxColorValue = 255),  # orange

  ## Other (greyish / blue-grey tones)
  "other"          = rgb(120, 120, 120, maxColorValue = 255),  # medium grey
  "other (Ki67+)"  = rgb( 80,  80,  80, maxColorValue = 255),  # dark grey
  "other (CXCR4+)" = rgb(195, 195, 195, maxColorValue = 255),  # mid grey
  "other (CD15+)"  = rgb( 96, 125, 139, maxColorValue = 255)   # blue-grey (was CHGA- CD15+ chromaffin progenitors)
)

# ---- AD cohort: metaclusters (analysis_AD/04_phenotyping.Rmd) ---------------
COL_METACLUSTER_AD <- c(
  "neural crest lineage" = "#238B45",  # green
  "mesenchymal"          = "#8C510A",  # brown
  "immune cells"         = "#6aa9d2",  # blue
  "endothelial"          = "#F59E0B",  # orange
  "other"                = "#969696"   # grey
)

# ---- Sample-level annotations: MYCN status / aberrations / progression /
#      community numbers (05.1_PT_DE_cellularcomm.Rmd, 10_spatial_analysis.Rmd) -
# Union of the keys used across scripts; color assignment is by name, so extra
# keys are harmless.
COL_CODE <- c(
  "MNA"     = rgb(202,   0,  32, maxColorValue = 255),
  "het"     = rgb(146, 197, 222, maxColorValue = 255),
  "nMNA"    = rgb(  5, 113, 176, maxColorValue = 255),
  "unknown" = "white",
  "yes"     = rgb(202,   0,  32, maxColorValue = 255),
  "no"      = rgb(  5, 113, 176, maxColorValue = 255),
  "Prog"    = rgb(202,   0,  32, maxColorValue = 255),
  "NoProg"  = rgb(  5, 113, 176, maxColorValue = 255),
  "1"       = "blue",
  "2"       = "red",
  "3"       = "green"
)

# ---- Fetahu et al. reference: RNA cluster colors (06.2_correlation_Fetahu.Rmd)
COL_FETAHU_CLUSTERS <- c(
  "B (3)"       = rgb(145,   2, 144, maxColorValue = 255),
  "B (16)"      = rgb(253, 218, 236, maxColorValue = 255),
  "B (19)"      = rgb(240,   0, 240, maxColorValue = 255),
  "B (11)"      = rgb(188, 128, 189, maxColorValue = 255),
  "B (7)"       = rgb(190, 174, 212, maxColorValue = 255),
  "M (1)"       = rgb( 28, 146, 255, maxColorValue = 255),
  "M (10)"      = rgb(  0, 255, 255, maxColorValue = 255),
  "M (15)"      = rgb( 39, 130, 187, maxColorValue = 255),
  "M (2)"       = rgb(166, 206, 227, maxColorValue = 255),
  "T (18)"      = rgb(227,  26,  28, maxColorValue = 255),
  "T (5)"       = rgb(254, 142, 175, maxColorValue = 255),
  "T (6)"       = rgb(242,   3, 137, maxColorValue = 255),
  "T (9)"       = rgb(251, 128, 114, maxColorValue = 255),
  "SC (14)"     = rgb(149, 144,  89, maxColorValue = 255),
  "SC (17)"     = rgb(210, 205, 126, maxColorValue = 255),
  "SC (20)"     = rgb(230, 245, 201, maxColorValue = 255),
  "NB (8)"      = rgb( 60, 119,  97, maxColorValue = 255),
  "E (13)"      = rgb(254, 217, 166, maxColorValue = 255),
  "NK (4)"      = rgb(212, 143,  75, maxColorValue = 255),
  "pDC (12)"    = rgb(253, 180,  98, maxColorValue = 255),
  "other (21)"  = rgb(174, 129,  37, maxColorValue = 255)
)

# ---- Fetahu et al. reference: predicted metacluster colors ------------------
COL_FETAHU_METACLUSTER <- c(
  "NB"    = "palegreen4",
  "NK"    = "#FB8072",
  "T"     = "indianred",
  "E"     = "#80B1D3",
  "M"     = "#8DD3C7",
  "B"     = "#BC80BD",
  "other" = "goldenrod1",
  "pDC"   = "tan4",
  "SC"    = "tan2"
)

# ---- Patient cohort overview tracks (01_read_data.Rmd cohort figure) -------
COL_YNH <- c(
  "yes" = rgb(202,   0,  32, maxColorValue = 255),  # red
  "no"  = rgb(  5, 113, 176, maxColorValue = 255),  # blue
  "het" = rgb(146, 197, 222, maxColorValue = 255)   # light blue
)
COL_SEX <- c(
  "f"   = rgb(202,   0,  32, maxColorValue = 255),  # red
  "m"   = rgb(  5, 113, 176, maxColorValue = 255)   # blue
)

# ---- Lee et al. reference: extra label colors (06.4_correlation_Lee.Rmd) ----
# Added on top of the per-celltype colors read from Lee/celltype_metacluster.csv.
COL_LEE_EXTRA <- c("tumor" = "snow2", "other" = "goldenrod1")

# ---- Generic qualitative palette for random sample/cluster colors (02_QC_1) -
# All named R colors except greys, used to assign arbitrary colors to samples /
# clusters where the specific hue does not carry meaning.
COL_VECTOR_433 <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = TRUE)]
