# =============================================================================
# Shared cluster bubble-plot helper for the phenotyping notebooks (PT/BM and AD).
#
# plot_cluster_bubble() draws a per-cluster bubble plot; it auto-detects tissue
# (per-tissue count barplots when colData(spe) has a `tissue` column, else one).
# Expects `spe`, `mat`, `plot_dir`, `date_prefix` in the calling environment.
# =============================================================================

library(ComplexHeatmap)
library(circlize)
library(grid)
library(viridis)

# Base seed for the (cosmetic) random cluster palettes.
color_seed_base <- 221228

# Canonical ordering of cluster levels (numeric if all-digit, else alphabetical).
get_cluster_levels <- function(values) {
  cluster_levels <- unique(as.character(values))
  cluster_levels <- cluster_levels[!is.na(cluster_levels)]
  cluster_levels <- cluster_levels[nzchar(cluster_levels)]

  if (length(cluster_levels) && all(grepl("^\\d+$", cluster_levels))) {
    cluster_levels <- as.character(sort(as.integer(cluster_levels)))
  } else {
    cluster_levels <- sort(cluster_levels)
  }

  cluster_levels
}

# Deterministic random palette for a set of cluster levels.
generate_cluster_colors <- function(cluster_levels, seed = color_seed_base) {
  non_gray_colors <- colors()[
    grep("\\bgr[ae]y", grDevices::colors(), invert = TRUE)
  ]

  set.seed(seed)
  setNames(
    sample(
      non_gray_colors,
      length(cluster_levels),
      replace = length(cluster_levels) > length(non_gray_colors)
    ),
    cluster_levels
  )
}

# S3 print method: no-op so the caller's print() does not re-draw.
print.heatmap_drawn <- function(x, ...) invisible(NULL)

plot_cluster_bubble <- function(cluster_col, save_name,
                                annotation_level = c("cluster", "metacluster"),
                                col_annotation  = NULL,
                                metacluster_col = "metacluster",
                                ncells_break_at = 50000,
                                column_title    = "",
                                plot_width      = 16,
                                plot_height     = 14,
                                barplot_height  = NULL) {
  annotation_level <- match.arg(annotation_level)

  # Plot size: bar/dendrogram heights scale with plot_height (defaults reproduce
  # 4 cm bar / 1.5 cm dendrogram); barplot_height overrides.
  if (is.null(barplot_height)) {
    barplot_height <- unit(plot_height * 4 / 14, "cm")
  } else if (!inherits(barplot_height, "unit")) {
    barplot_height <- unit(barplot_height, "cm")
  }
  dend_height <- unit(plot_height * 1.5 / 14, "cm")

  # Cluster vector (+ optional tissue split)
  cluster_vec <- as.character(colData(spe)[[cluster_col]])
  tissue_vec  <- if ("tissue" %in% names(colData(spe))) {
    as.character(colData(spe)$tissue)
  } else {
    NULL
  }

  # Expression matrix
  mat_df <- as.data.frame(t(mat))

  exp_mat_ <- aggregate(
    mat_df,
    by = list(cluster = cluster_vec),
    FUN = function(x) {
      pos_vals <- x[x > 0]
      if (length(pos_vals) == 0) return(0)
      median(log1p(pos_vals))
    }
  )

  cluster_names <- exp_mat_$cluster

  exp_mat <- exp_mat_[, -1, drop = FALSE]
  rownames(exp_mat) <- cluster_names
  exp_mat <- t(as.matrix(exp_mat))

  # Row-scale expression to 0-1
  rmin <- apply(exp_mat, 1, min, na.rm = TRUE)
  rmax <- apply(exp_mat, 1, max, na.rm = TRUE)
  den  <- rmax - rmin

  exp_mat <- (exp_mat - rmin) / den

  exp_mat[is.nan(exp_mat)] <- 0
  exp_mat[is.infinite(exp_mat)] <- 0

  clusters <- colnames(exp_mat)

  # Fraction of cells expressing each marker, per cluster
  percent_mat <- sapply(seq_len(nrow(mat)), function(i) {

    row_vals <- mat[i, ]
    pos <- row_vals[row_vals > 0]

    if (length(pos) == 0) {

      out <- rep(0, length(clusters))
      names(out) <- clusters
      return(out)
    }

    thresh <- quantile(pos, 0.2, na.rm = TRUE)

    row_bin <- ifelse(row_vals > thresh, 1, 0)

    tapply(
      row_bin,
      factor(cluster_vec, levels = clusters),
      mean
    )
  })

  percent_mat <- t(percent_mat)

  percent_mat <- percent_mat[, clusters, drop = FALSE]

  percent_mat[is.na(percent_mat)] <- 0
  percent_mat[percent_mat < 0.1] <- 0.1

  # Cell counts (per tissue if a `tissue` column is present, else one group)
  if (is.null(tissue_vec)) {
    count_list <- list(AD = table(factor(cluster_vec, levels = clusters)))
    bar_fill   <- c(AD = "goldenrod3")
  } else {
    count_list <- list(
      PT = table(factor(cluster_vec[tissue_vec == "PT"], levels = clusters)),
      BM = table(factor(cluster_vec[tissue_vec == "BM"], levels = clusters))
    )
    bar_fill   <- c(PT = "brown", BM = "aquamarine4")
  }

  # Barplot info helper (optional axis break)
  make_barplot_info <- function(v) {
    v_num   <- as.numeric(v)
    max_val <- max(v_num, na.rm = TRUE)

    if (max_val <= ncells_break_at) {
      y_axis_max <- max(10000, ceiling(max_val / 10000) * 10000)
      at <- seq(0, y_axis_max, by = 10000)
      return(list(
        values      = v_num,
        upper       = y_axis_max,
        axis_param  = list(
          at     = at,
          labels = ifelse(at == 0, "0", paste0(at / 1000, "k")),
          gp     = gpar(fontsize = 12)
        ),
        needs_break = FALSE
      ))
    }

    # Break case: compress above-break bars into a small top band
    top_frac   <- 0.18
    break_frac <- 1 - top_frac
    top_ext    <- ncells_break_at * (top_frac / break_frac)
    ylim_max   <- ncells_break_at + top_ext

    # Above-break bars are drawn at 90% of the top band
    above_display <- ncells_break_at + top_ext * 0.9
    values_plot   <- ifelse(v_num > ncells_break_at, above_display, v_num)

    at_below   <- seq(0, ncells_break_at, by = 10000)
    at_all     <- c(at_below, above_display)
    labels_all <- c(
      ifelse(at_below == 0, "0", paste0(at_below / 1000, "k")),
      paste0(round(max_val / 1000), "k")
    )

    list(
      values           = values_plot,
      upper            = ylim_max,
      axis_param       = list(
        at     = at_all,
        labels = labels_all,
        gp     = gpar(fontsize = 12)
      ),
      needs_break      = TRUE,
      break_y_npc      = ncells_break_at / ylim_max,
      above_break_orig = which(v_num > ncells_break_at),
      n_total          = length(v_num)
    )
  }

  # Cluster / metacluster colours
  if (annotation_level == "metacluster") {
    # one metacluster per cluster (aggregate with max, as in the original)
    metaclust_by_cluster <- aggregate(
      colData(spe)[[metacluster_col]],
      list(cluster = cluster_vec),
      max
    )
    meta_by_cluster <- setNames(
      as.character(metaclust_by_cluster$x),
      as.character(metaclust_by_cluster$cluster)
    )
    annot_values <- meta_by_cluster[clusters]
    annot_name <- "metacluster"
    annot_colors <- metadata(spe)$color_vectors$col_metacluster
  } else {
    color_key <- paste0("col_", cluster_col)
    # deterministic per-name seed offset
    seed_offset <- sum(utf8ToInt(cluster_col))

    existing_colors <- metadata(spe)$color_vectors[[color_key]]
    if (is.null(existing_colors) ||
        !all(clusters %in% names(existing_colors)) ||
        any(is.na(existing_colors[clusters]))) {
      existing_colors <- generate_cluster_colors(
        cluster_levels = clusters,
        seed = color_seed_base + seed_offset
      )
      metadata(spe)$color_vectors[[color_key]] <- existing_colors
    }

    annot_values <- clusters
    annot_name <- cluster_col
    annot_colors <- existing_colors[clusters]
  }

  if (!is.null(col_annotation)) {
    annot_colors <- col_annotation
  }

  levels_from_colors <- names(annot_colors)
  annot_levels <- if (!is.null(levels_from_colors)) {
    intersect(levels_from_colors, unique(annot_values))
  } else {
    unique(annot_values)
  }
  if (length(annot_levels) == 0) {
    annot_levels <- unique(annot_values)
  }

  col_for_annot <- setNames(list(annot_colors), annot_name)

  info_list <- lapply(count_list, make_barplot_info)
  cluster_annotation_list <- setNames(
    list(factor(annot_values, levels = annot_levels)),
    annot_name
  )

  col_fun <- circlize::colorRamp2(
    c(0, 0.5, 1),
    viridis(100)[c(1, 50, 100)]
  )

  # Column annotation (one count barplot per tissue group)
  bar_annos <- setNames(
    lapply(names(info_list), function(nm) {
      info <- info_list[[nm]]
      anno_barplot(
        info$values,
        ylim       = c(0, info$upper),
        height     = barplot_height,
        gp         = gpar(fill = bar_fill[[nm]]),
        axis_param = info$axis_param
      )
    }),
    paste0("ncells_", names(info_list))
  )

  column_ha <- do.call(
    HeatmapAnnotation,
    c(
      cluster_annotation_list,
      bar_annos,
      list(
        col = col_for_annot,
        show_legend = TRUE
      )
    )
  )

  # Bubble cell: circle sized by percent expressed, coloured by expression
  cell_fun <- function(j, i, x, y, w, h, fill) {

    grid.rect(
      x, y, w, h,
      gp = gpar(col = NA, fill = NA)
    )

    grid.circle(
      x = x,
      y = y,
      r = percent_mat[i, j] / 2 * min(unit.c(w, h)),
      gp = gpar(
        fill = col_fun(exp_mat[i, j]),
        col = NA
      )
    )
  }

  lgd_list <- list(

    Legend(
      # labels shown as percentages; circle radii stay on the 0-1 scale
      labels = c(0, 25, 50, 75, 100),
      title = "percentage expressed",

      graphics = lapply(
        c(0, 0.25, 0.5, 0.75, 1),
        function(r) {

          function(x, y, w, h)
            grid.circle(
              x = x,
              y = y,
              r = r * unit(2, "mm"),
              gp = gpar(fill = "black")
            )
        }
      )
    )
  )

  hp <- Heatmap(
    exp_mat,

    heatmap_legend_param = list(title = "expression"),

    column_title = column_title,

    top_annotation = column_ha,

    col = col_fun,

    rect_gp = gpar(type = "none"),

    clustering_distance_columns = "pearson",

    clustering_method_columns = "average",

    row_labels = sapply(
      strsplit(rownames(exp_mat), "_"),
      "[[",
      1
    ),

    cell_fun = cell_fun,

    show_row_dend = FALSE,

    column_dend_height = dend_height,

    row_names_gp = gpar(fontsize = 15),

    column_names_gp = gpar(fontsize = 15),

    column_names_rot = 45,

    # raise the column-label height cap so long rotated labels aren't truncated
    column_names_max_height = max_text_width(
      colnames(exp_mat),
      gp = gpar(fontsize = 15)
    ),

    border = "black"
  )

  # Break-mark decoration helper
  add_break_marks <- function(anno_name, info, drawn_hp) {
    if (!info$needs_break) return(invisible(NULL))

    col_ord <- column_order(drawn_hp)
    if (is.list(col_ord)) col_ord <- unlist(col_ord)
    n_tot <- info$n_total
    b_y   <- info$break_y_npc
    gap   <- 0.04   # white erasure gap height (npc)

    decorate_annotation(anno_name, {
      pushViewport(viewport(clip = "off"))

      for (orig_idx in info$above_break_orig) {
        vis <- match(orig_idx, col_ord)
        if (is.na(vis)) next
        xc <- (vis - 0.5) / n_tot
        bw <- 0.8  / n_tot

        # white gap to break the bar
        grid.rect(
          x = xc, y = b_y - gap / 2,
          width = bw, height = gap,
          just  = c("center", "bottom"),
          default.units = "npc",
          gp = gpar(fill = "white", col = NA)
        )

        # zigzag across the bar
        xs <- c(xc - bw/2, xc - bw/4, xc, xc + bw/4, xc + bw/2)
        ys <- b_y + c(-gap * 0.6, gap * 0.6, -gap * 0.6, gap * 0.6, -gap * 0.6)
        grid.polyline(
          x = xs, y = ys,
          default.units = "npc",
          gp = gpar(col = "black", lwd = 1.5)
        )
      }

      # zigzag break mark on the y-axis
      grid.polyline(
        x = unit(c(-2.5, -1.5, -2.5), "mm"),
        y = unit(c(b_y - gap * 0.6, b_y, b_y + gap * 0.6), "npc"),
        gp = gpar(col = "black", lwd = 1.5)
      )

      popViewport()
    })
  }

  # Draw, decorate, save
  full_path <- file.path(plot_dir, paste0(date_prefix, "_", save_name))
  pdf(full_path, width = plot_width, height = plot_height)
  hp_drawn <- draw(hp, annotation_legend_list = lgd_list,
                   padding = unit(c(2, 55, 2, 2), "mm"))
  for (nm in names(info_list)) {
    add_break_marks(paste0("ncells_", nm), info_list[[nm]], hp_drawn)
  }
  dev.off()
  message("Saved: ", full_path)

  # Draw again for interactive display; return a silent marker so the caller's
  # print() does not trigger another re-draw.
  hp_drawn <- draw(hp, annotation_legend_list = lgd_list,
                   padding = unit(c(2, 55, 2, 2), "mm"))
  for (nm in names(info_list)) {
    add_break_marks(paste0("ncells_", nm), info_list[[nm]], hp_drawn)
  }

  invisible(structure(list(), class = "heatmap_drawn"))
}
