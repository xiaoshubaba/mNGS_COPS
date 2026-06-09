library(ggplot2)
library(dplyr)
library(ggsci)
library(patchwork)

pdf(NULL)

config <- read.delim("mul.cfg.txt", header = TRUE)
config_project_order <- unique(config$project)
n_projects <- length(config_project_order)

all_data <- list()
for (i in 1:nrow(config)) {
  proj <- config$project[i]
  file_path <- config$status[i]
  source_name <- config$source[i]
  sample_count <- config$number[i]
  cat(sprintf("读取: %s -> %s (%s)\n", proj, source_name, file_path))
  if (!file.exists(file_path)) {
    warning(sprintf("文件不存在: %s", file_path))
    next
  }
  file_data <- read.delim(file_path, comment.char = "#", header = TRUE)
  if ("adjusted_signal_ratio" %in% names(file_data)) {
    signal_col <- "adjusted_signal_ratio"
  } else if ("adjust_signal_ratio" %in% names(file_data)) {
    signal_col <- "adjust_signal_ratio"
  } else {
    warning("No signal column")
    next
  }
  if ("#sample" %in% names(file_data)) {
    sample_col <- "#sample"
  } else if ("sample" %in% names(file_data)) {
    sample_col <- "sample"
  } else {
    warning("No sample column")
    next
  }
  extracted_data <- data.frame(
    project = proj, source = source_name, number = sample_count,
    sample = file_data[[sample_col]], signal_ratio = file_data[[signal_col]],
    stringsAsFactors = FALSE
  )
  all_data[[length(all_data) + 1]] <- extracted_data
}
df_final <- do.call(rbind, all_data)
if (is.null(df_final) || nrow(df_final) == 0) stop("No data")
df_final$project <- factor(df_final$project, levels = config_project_order)

plot_list <- list()

# 保持您认为宽度合适的原始参数（高度仍为原始值，后续通过调整p值偏移解决重叠）
if (n_projects <= 3) {
  base_width_per_plot <- 3.8
  base_height <- 6.0
  text_size_factor <- 1.0
  jitter_width <- 0.12
  point_size <- 0.9
} else if (n_projects <= 6) {
  base_width_per_plot <- 3.5
  base_height <- 5.5
  text_size_factor <- 0.9
  jitter_width <- 0.15
  point_size <- 0.8
} else if (n_projects <= 9) {
  base_width_per_plot <- 3.2
  base_height <- 5.0
  text_size_factor <- 0.85
  jitter_width <- 0.18
  point_size <- 0.7
} else {
  base_width_per_plot <- 3.0
  base_height <- 4.5
  text_size_factor <- 0.8
  jitter_width <- 0.2
  point_size <- 0.6
}

if (n_projects <= 3) {
  ncol_layout <- n_projects
  nrow_layout <- 1
} else if (n_projects <= 6) {
  ncol_layout <- 3
  nrow_layout <- ceiling(n_projects / 3)
} else if (n_projects <= 9) {
  ncol_layout <- 3
  nrow_layout <- ceiling(n_projects / 3)
} else {
  ncol_layout <- 4
  nrow_layout <- ceiling(n_projects / 4)
}

cat(sprintf("项目数量: %d, 布局: %d列 x %d行\n", n_projects, ncol_layout, nrow_layout))

for (i in seq_along(config_project_order)) {
  proj <- config_project_order[i]
  proj_config <- config %>% filter(project == proj)
  proj_data <- df_final %>% filter(project == proj)
  if (nrow(proj_data) == 0) next

  source_order <- proj_config$source
  proj_data$x_lab <- ifelse(proj_data$source == "Base", "Base",
                            paste0(proj_data$source, "\n(n=", proj_data$number, ")"))
  x_lab_levels <- c()
  for (src in source_order) {
    if (src == "Base") {
      x_lab_levels <- c(x_lab_levels, "Base")
    } else {
      n_val <- proj_config$number[proj_config$source == src]
      if (!is.na(n_val) && n_val != "NA") {
        x_lab_levels <- c(x_lab_levels, paste0(src, "\n(n=", n_val, ")"))
      } else {
        x_lab_levels <- c(x_lab_levels, src)
      }
    }
  }
  proj_data$x_lab <- factor(proj_data$x_lab, levels = x_lab_levels)
  proj_data$source <- factor(proj_data$source, levels = source_order)
  show_y_label <- (i %% ncol_layout == 1) || (ncol_layout == 1)

  p <- ggplot(proj_data, aes(x = x_lab, y = signal_ratio, fill = source)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.35, linewidth = 0.4, color = "black") +
    geom_jitter(width = jitter_width, height = 0, color = "grey60", alpha = 0.35, size = point_size) +
    scale_fill_nejm() +
    labs(x = NULL, y = ifelse(show_y_label, "Signal Ratio (%)", ""), caption = proj) +
    theme_minimal(base_size = 10 * text_size_factor, base_family = "sans") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9 * text_size_factor),
      axis.text.y = element_text(size = 9 * text_size_factor),
      axis.title.y = element_text(size = 10 * text_size_factor, margin = margin(r = 5)),
      plot.caption = element_text(hjust = 0.5, size = 11 * text_size_factor, margin = margin(t = 5), face = "bold"),
      legend.position = "none",
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.8),
      plot.margin = margin(8, 8, 10, 8)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15)))

  x_levels <- levels(proj_data$x_lab)
  base_pos <- which(x_levels == "Base")
  if (length(base_pos) > 0) {
    p <- p + geom_vline(xintercept = base_pos + 0.5, linetype = "dashed", color = "gray50", linewidth = 0.6, alpha = 0.7)
  }

  # ==================== 修改后的p值添加部分 ====================
  if ("Merge" %in% source_order && length(source_order) > 2) {
    other_sources <- setdiff(source_order, c("Base", "Merge"))
    comparison_max_vals <- c()
    for (src in other_sources) {
      merge_data <- proj_data$signal_ratio[proj_data$source == "Merge"]
      src_data <- proj_data$signal_ratio[proj_data$source == src]
      if (length(merge_data) > 0 && length(src_data) > 0) {
        comparison_max <- max(c(merge_data, src_data), na.rm = TRUE)
        comparison_max_vals <- c(comparison_max_vals, comparison_max)
      }
    }
    if (length(comparison_max_vals) > 0) {
      global_y_max <- max(comparison_max_vals, na.rm = TRUE)
      y_offset <- global_y_max * 0.04   # 减小偏移
      for (j in seq_along(other_sources)) {
        src <- other_sources[j]
        merge_data <- proj_data$signal_ratio[proj_data$source == "Merge"]
        src_data <- proj_data$signal_ratio[proj_data$source == src]
        if (length(merge_data) > 0 && length(src_data) > 0) {
          p_val <- wilcox.test(merge_data, src_data, paired = FALSE)$p.value
          merge_label <- as.character(unique(proj_data$x_lab[proj_data$source == "Merge"]))
          src_label <- as.character(unique(proj_data$x_lab[proj_data$source == src]))
          merge_pos <- which(x_levels == merge_label)
          src_pos <- which(x_levels == src_label)
          if (length(merge_pos) > 0 && length(src_pos) > 0) {
            x_pos <- mean(c(merge_pos, src_pos))
            y_pos <- global_y_max + y_offset * (0.5 + (j - 1) * 0.4)
            if (p_val < 0.0001) {
              p_label <- "p < 0.0001"
            } else if (p_val < 0.001) {
              p_label <- "p < 0.001"
            } else if (p_val < 0.01) {
              p_label <- sprintf("p = %.3f", p_val)
            } else {
              p_label <- sprintf("p = %.3f", p_val)
            }
            p <- p +
              annotate("text", x = x_pos, y = y_pos, label = p_label,
                       size = 2.0 * text_size_factor, vjust = 0, family = "sans")
            bracket_height <- y_offset * 0.2
            p <- p +
              annotate("segment", x = merge_pos, xend = src_pos,
                       y = y_pos - bracket_height, yend = y_pos - bracket_height,
                       color = "black", linewidth = 0.3) +
              annotate("segment", x = merge_pos, xend = merge_pos,
                       y = y_pos - bracket_height, yend = y_pos - bracket_height * 1.5,
                       color = "black", linewidth = 0.3) +
              annotate("segment", x = src_pos, xend = src_pos,
                       y = y_pos - bracket_height, yend = y_pos - bracket_height * 1.5,
                       color = "black", linewidth = 0.3)
          }
        }
      }
      if (length(other_sources) > 0) {
        max_y_pos <- global_y_max + y_offset * (0.5 + (length(other_sources) - 1) * 0.4) + y_offset * 0.3
        p <- p + coord_cartesian(ylim = c(0, max_y_pos), clip = "off")
      }
    }
  }
  # ==================== 修改结束 ====================

  plot_list[[proj]] <- p
}

if (length(plot_list) > 0) {
  plot_list <- plot_list[config_project_order]
  final_plot <- wrap_plots(plot_list, ncol = ncol_layout, nrow = nrow_layout, guides = "collect") +
    plot_layout(guides = 'collect') & theme(legend.position = 'none')
  total_width <- base_width_per_plot * min(n_projects, ncol_layout)
  total_height <- base_height * nrow_layout
  cat(sprintf("\n最终图形尺寸: %.1f英寸 x %.1f英寸\n", total_width, total_height))
  ggsave(sprintf("signal_ratio_analysis_%dprojects.png", n_projects), final_plot,
         width = total_width, height = total_height, dpi = 300, bg = "white")
  ggsave(sprintf("signal_ratio_analysis_%dprojects.pdf", n_projects), final_plot,
         width = total_width, height = total_height, dpi = 300, bg = "white", device = "pdf")
  cat("\n✓ 图形已保存\n")
}
