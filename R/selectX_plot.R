###
# Author  : Arvind Iyer, Miljan Petrovic
# Project : SelectSim
# Desc    : The file contains plot related functions
# Version : 0.1.6
###

# Suppress R CMD check notes for ggplot2/ggridges column name variables
utils::globalVariables(c(
  "density", "fill", "freq", "group", "iscale",
  "name", "pairs", "scale", "x", "x_actual", "y", "ymin"
))


#' A clean ggplot2 theme for publication-quality plots
#'
#' @param base_size size of fonts
#' @param base_family type of font
#' @return A ggplot2 theme object.
#'
#' @examples
#' library(ggplot2)
#' ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) +
#'   geom_point() + theme_Publication()
#'
#' @export
theme_Publication <- function(base_size = 14, base_family = "sans") {
  (theme(
    text = element_text(family = base_family, size = base_size),
    plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.title = element_text(face = "bold", size = rel(1)),
    axis.title.y = element_text(angle = 90, vjust = 2),
    axis.title.x = element_text(vjust = -0.2),
    axis.text = element_text(size = rel(0.8)),
    axis.ticks = element_line(),
    panel.grid.major = element_line(colour = "#f0f0f0"),
    panel.grid.minor = element_blank(),
    legend.key = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = rel(0.8)),
    legend.key.size = unit(4, "mm"),
    legend.direction = "vertical",
    legend.spacing = unit(0, "cm"),
    legend.title = element_text(face = "italic"),
    plot.margin = margin(5, 5, 5, 5, "mm"),
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_text(face = "bold")
  ))
}


#' Scatter plot of observed vs expected weighted co-mutation
#'
#' Plots each gene pair as a point with observed weighted co-mutation on the y-axis
#' and expected (null model mean) on the x-axis. Significant co-mutations (CO) and
#' mutual exclusivities (ME) are coloured; non-significant pairs are grey.
#'
#' @import ggpubr
#' @import ggplot2
#' @importFrom  dplyr %>%
#' @importFrom  dplyr mutate
#' @importFrom  dplyr case_when
#' @param result Result table from \code{selectX()$result}.
#' @param title Plot title.
#' @return A ggplot2 object.
#'
#' @examples
#' \donttest{
#' data(luad_result, package = "SelectSim")
#' obs_exp_scatter(result = luad_result, title = "TCGA LUAD")
#' }
#'
#' @export

obs_exp_scatter <- function(result, title) {
  result$log_overlap <- log10(result$w_overlap + 1)
  result$log_r_overlap <- log10(result$w_r_overlap + 1)
  result <- result %>% mutate(cat = case_when(FDR == TRUE ~ type, FDR == FALSE ~ "NS"))

  plot <- ggscatter(result,
    x = "log_r_overlap",
    y = "log_overlap",
    color = "cat",
    palette = c("forestgreen", "purple", "grey93"),
    repel = TRUE,
    alpha = 1
  ) + xlab("Random weighted co-mutation (log10)") + ylab("Actual weighted co-mutation(log10)") + geom_abline(slope = 1, intercept = 0) +
    scale_y_continuous(breaks = seq(0, max(result$log_overlap), 1)) +
    scale_x_continuous(breaks = seq(0, max(result$log_overlap), 1))

  plot <- plot + ggtitle(title) + theme_Publication(base_size = 16) +
    theme(
      legend.position = "top",
      legend.title = element_blank()
    )

  return(plot)
}


#' Extract null-model weighted overlap distribution for a gene pair
#'
#' Returns the vector of weighted co-mutation values from all null-model
#' permutations for a specific pair of genes.
#'
#' @param gene1 Name of the first gene/alteration.
#' @param gene2 Name of the second gene/alteration.
#' @param obj SelectX object (list returned by \code{selectX()$obj}).
#' @return Numeric vector of length \code{obj$nSim} with the null-model
#'   weighted overlap values for the pair.
#'
#' @examples
#' \donttest{
#' data(luad_run_data, package = "SelectSim")
#' result <- selectX(M = luad_run_data$M,
#'                   sample.class = luad_run_data$sample.class,
#'                   alteration.class = luad_run_data$alteration.class,
#'                   n.cores = 1, min.freq = 10, n.permut = 10,
#'                   verbose = FALSE)
#' genes <- rownames(result$obj$al$am$full)[1:2]
#' overlap_pair_extract(genes[1], genes[2], result$obj)
#' }
#'
#' @export

overlap_pair_extract <- function(gene1, gene2, obj) {
  wes_dist <- c()
  panel_dist <- c()
  for (i in c(1:obj$nSim)) {
    rownames(obj$wrobs.co[[i]]) <- rownames(obj$al$am$full)
    colnames(obj$wrobs.co[[i]]) <- rownames(obj$al$am$full)
    wes_dist <- c(wes_dist, obj$wrobs.co[[i]][gene1, gene2])
  }
  return(wes_dist)
}


#' Ridge plot of null-model background distribution for significant gene pairs
#'
#' For each significant evolutionary dependency in \code{result_df}, plots the
#' null-model weighted-overlap distribution as a ridge, with vertical lines
#' marking the observed overlap (red) and mean background (blue).
#'
#' @import ggpubr
#' @import ggplot2
#' @import ggridges
#' @importFrom  dplyr %>%
#' @importFrom  dplyr case_when
#' @importFrom  dplyr group_by
#' @importFrom  dplyr ungroup
#' @importFrom  dplyr mutate
#' @importFrom  dplyr row_number
#' @importFrom  dplyr filter
#' @importFrom  dplyr left_join
#' @importFrom  dplyr arrange
#' @importFrom  dplyr distinct
#' @importFrom  stats as.dist
#' @importFrom  stats density
#' @importFrom  stats ecdf
#' @importFrom  stats runif
#' @importFrom  stats setNames
#' @param result_df Subset of the result table from \code{selectX()$result} containing
#'   only the pairs you want to display (e.g., significant hits).
#' @param obj SelectX object (\code{selectX()$obj}).
#' @return A ggplot2 ridge-plot object.
#'
#' @examples
#' \donttest{
#' data(luad_run_data, package = "SelectSim")
#' result <- selectX(M = luad_run_data$M,
#'                   sample.class = luad_run_data$sample.class,
#'                   alteration.class = luad_run_data$alteration.class,
#'                   n.cores = 1, min.freq = 10, n.permut = 10,
#'                   verbose = FALSE)
#' sig_pairs <- head(result$result[result$result$FDR, ], 3)
#' if (nrow(sig_pairs) > 0) ridge_plot_ed(sig_pairs, result$obj)
#' }
#'
#' @export

ridge_plot_ed <- function(result_df, obj) {
  pair_mat <- matrix(0, nrow = nrow(result_df), ncol = obj$nSim)
  for (row in c(1:nrow(result_df))) {
    gene1 <- (result_df[row, "SFE_1"])
    gene2 <- (result_df[row, "SFE_2"])
    pair <- (result_df[row, "name"])
    pair_mat[row, ] <- overlap_pair_extract(gene1, gene2, obj)
  }
  pair_mat_plot <- t(pair_mat)
  colnames(pair_mat_plot) <- result_df$name
  df <- setNames(as.data.frame.table(pair_mat_plot, stringsAsFactors = FALSE), c("rows", "pairs", "freq"))
  p <- ggplot(df, aes(x = freq, y = pairs)) +
    xlab("Co-mutated samples") +
    geom_density_ridges_gradient(quantile_lines = FALSE, quantiles = 2, color = "black", fill = "lightgrey") +
    theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
    scale_x_continuous(breaks = seq(0, max(df$freq), by = 2))

  q <- ggplot_build(p)$data[[1]]
  density_lines <- q %>%
    group_by(group) %>%
    dplyr::filter(density == max(density)) %>%
    ungroup()
  temp <- df %>%
    distinct(pairs) %>%
    mutate(number = row_number())
  density_lines_complete <- left_join(density_lines, temp, by = c("group" = "number"))
  density_lines_complete$x_actual <- result_df$w_overlap

  a <- ifelse(result_df$type == "ME", "purple", "forestgreen")
  plot <- ggplot(
    df,
    aes(
      x = freq,
      y = pairs
    )
  ) +
    xlab("Weighted overlap") +
    ylab("") +
    geom_density_ridges_gradient(quantile_lines = FALSE, quantiles = 2, color = "black", fill = "lightgrey") +
    theme_ridges(grid = TRUE, center_axis_labels = TRUE, font_size = 18) +
    scale_x_continuous(breaks = seq(0, max(df$freq), by = 5)) +
    geom_segment(data = density_lines_complete, aes(x = x_actual, xend = x_actual, y = ymin, yend = ymin + density * scale * iscale, color = "Actual overlap")) +
    geom_segment(data = density_lines_complete, aes(x = x, xend = x, y = ymin, yend = ymin + density * scale * iscale, color = "Mean Background overlap")) +
    scale_color_manual("Legend", values = c("red", "blue")) +
    theme(axis.text.y = element_text(colour = a))

  return(plot)
}

#' Ridge plot comparing null-model distributions for two datasets
#'
#' For each gene pair in \code{result_df}, overlays the null-model weighted-overlap
#' distributions from two separate \code{selectX} runs, allowing visual comparison
#' of evolutionary dependencies across cohorts or conditions.
#'
#' @import ggpubr
#' @import ggplot2
#' @import ggridges
#' @importFrom  dplyr %>%
#' @importFrom  dplyr case_when
#' @importFrom  dplyr group_by
#' @importFrom  dplyr ungroup
#' @importFrom  dplyr mutate
#' @importFrom  dplyr row_number
#' @importFrom  dplyr filter
#' @importFrom  dplyr left_join
#' @importFrom  dplyr arrange
#' @importFrom  stats as.dist
#' @importFrom  stats density
#' @importFrom  stats ecdf
#' @importFrom  stats runif
#' @importFrom  stats setNames
#' @param result_df A data frame of gene pairs to display, with columns
#'   \code{SFE_1}, \code{SFE_2}, \code{name}, \code{type}, \code{dataset1_w_overlap},
#'   and \code{dataset2_w_overlap}.
#' @param obj1 SelectX object for dataset 1 (\code{selectX()$obj}).
#' @param obj2 SelectX object for dataset 2 (\code{selectX()$obj}).
#' @param name1 Display label for dataset 1.
#' @param name2 Display label for dataset 2.
#' @return A ggplot2 ridge-plot object.
#'
#' @examples
#' \donttest{
#' data(luad_run_data, package = "SelectSim")
#' r1 <- selectX(M = luad_run_data$M,
#'               sample.class = luad_run_data$sample.class,
#'               alteration.class = luad_run_data$alteration.class,
#'               n.cores = 1, min.freq = 10, n.permut = 10, verbose = FALSE)
#' r2 <- selectX(M = luad_run_data$M,
#'               sample.class = luad_run_data$sample.class,
#'               alteration.class = luad_run_data$alteration.class,
#'               n.cores = 1, min.freq = 10, n.permut = 10, verbose = FALSE)
#' common <- head(r1$result, 3)
#' common$dataset1_w_overlap <- common$w_overlap
#' common$dataset2_w_overlap <- common$w_overlap
#' ridge_plot_ed_compare(common, r1$obj, r2$obj, "Run1", "Run2")
#' }
#'
#' @export

ridge_plot_ed_compare <- function(result_df, obj1, obj2, name1, name2) {
  pair_mat <- matrix(0, nrow = nrow(result_df), ncol = obj1$nSim)
  for (row in c(1:nrow(result_df))) {
    gene1 <- (result_df[row, "SFE_1"])
    gene2 <- (result_df[row, "SFE_2"])
    pair <- (result_df[row, "name"])
    pair_mat[row, ] <- overlap_pair_extract(gene1, gene2, obj1)
  }
  pair_mat_plot <- t(pair_mat)
  colnames(pair_mat_plot) <- result_df$name
  df1 <- setNames(as.data.frame.table(pair_mat_plot, stringsAsFactors = FALSE), c("rows", "pairs", "freq"))

  pair_mat <- matrix(0, nrow = nrow(result_df), ncol = obj2$nSim)
  for (row in c(1:nrow(result_df))) {
    gene1 <- (result_df[row, "SFE_1"])
    gene2 <- (result_df[row, "SFE_2"])
    pair <- (result_df[row, "name"])
    pair_mat[row, ] <- overlap_pair_extract(gene1, gene2, obj2)
  }
  pair_mat_plot <- t(pair_mat)
  colnames(pair_mat_plot) <- result_df$name
  df2 <- setNames(as.data.frame.table(pair_mat_plot, stringsAsFactors = FALSE), c("rows", "pairs", "freq"))

  df1$name <- rep(name1, nrow(df1))
  df2$name <- rep(name2, nrow(df2))

  df <- rbind(df1, df2)

  p <- ggplot(df, aes(x = freq, y = pairs, fill = name)) +
    xlab("Co-mutated samples") +
    geom_density_ridges_gradient(quantile_lines = FALSE, quantiles = 2, color = "black") +
    scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
    theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
    scale_x_continuous(breaks = seq(0, max(df$freq), by = 5))

  q <- ggplot_build(p)$data[[1]]
  density_lines <- q %>%
    group_by(fill, y) %>%
    dplyr::filter(density == max(density)) %>%
    arrange(fill) %>%
    ungroup()

  density_lines$pairs <- c(result_df$name, result_df$name)
  density_lines$x_actual <- c(result_df$dataset1_w_overlap, result_df$dataset2_w_overlap)
  density_lines$name <- c(rep(name2, nrow(result_df)), rep(name1, nrow(result_df)))
  a <- ifelse(result_df$type == "ME", "purple", "forestgreen")
  plot <- ggplot(df,
    aes(
      x = freq,
      y = pairs,
      fill = name
    ),
    alpha = 0.2
  ) +
    xlab("Weighted overlap") +
    ylab("") +
    geom_density_ridges_gradient(quantile_lines = FALSE, quantiles = 2) +
    theme_ridges(grid = TRUE, center_axis_labels = TRUE, font_size = 18) +
    scale_x_continuous(breaks = seq(0, max(df$freq), by = 5)) +
    geom_segment(data = density_lines, aes(x = x_actual, xend = x_actual, y = ymin, yend = ymin + density * scale * iscale, color = name, linetype = "Actual overlap")) +
    geom_segment(data = density_lines, aes(x = x, xend = x, y = ymin, yend = ymin + density * scale * iscale, color = name, linetype = "Mean Background overlap")) +
    scale_color_manual("Line Color", values = c("red", "blue")) +
    scale_linetype_manual("Linetype", values = c("Actual overlap" = 1, "Mean Background overlap" = 2)) +
    scale_fill_manual("Dataset", values = c("#E69F00", "#56B4E9")) +
    theme(axis.text.y = element_text(colour = a))

  return(plot)
}
