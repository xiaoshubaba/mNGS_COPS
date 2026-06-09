library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(ggsci)
library(scales)

# ==================== 1. Read project config file ====================
config_file <- "anno.cfg.txt"
config <- read.table(config_file, header = TRUE, sep = "", stringsAsFactors = FALSE)
cat("Found", nrow(config), "project configurations:\n")
print(config)

# ==================== 2. Define function with EARLY NA REMOVAL ====================
process_project_data <- function(project_folder, project_name) {
  cat("  Processing project:", project_name, "| Folder:", project_folder, "\n")
  if (!dir.exists(project_folder)) {
    warning(paste("    Directory does not exist, skipping:", project_folder))
    return(NULL)
  }
  # Read data files
  file_pattern <- "^rand\\..*\\.tabs\\.txt$"
  data_files <- list.files(path = project_folder,
                           pattern = file_pattern,
                           full.names = TRUE)
  if (length(data_files) == 0) {
    warning(paste("    No data files found in", project_folder, "skipping."))
    return(NULL)
  }
  all_data <- map_dfr(data_files, ~{
    df <- read.table(.x, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    df$source_file <- basename(.x)
    return(df)
  })
  cat("    Successfully read", length(data_files), "files, total", nrow(all_data), "rows.\n")
  # ===== EARLY NA REMOVAL =====
  # Remove rows with NA in case.F and control.F FIRST (as requested)
  rows_before_na <- nrow(all_data)
  # 1. Remove rows where case.F is NA
  all_data <- all_data %>%
    filter(!is.na(case.F))
  # 2. Remove rows where control.F is NA (early removal as requested)
  all_data <- all_data %>%
    filter(!is.na(control.F))
  rows_after_na <- nrow(all_data)
  na_removed <- rows_before_na - rows_after_na
  if (na_removed > 0) {
    cat("    Removed", na_removed, "rows with case.F = NA or control.F = NA\n")
  }
  # Now remove rows with NA in other key variables
  rows_before_other_na <- nrow(all_data)
  all_data_clean <- all_data %>%
    filter(!is.na(length) & 
           !is.na(case.mean.depth) & !is.na(control.mean.depth))
  rows_after_other_na <- nrow(all_data_clean)
  other_na_removed <- rows_before_other_na - rows_after_other_na
  if (other_na_removed > 0) {
    cat("    Removed", other_na_removed, "rows with NA in length or depth variables\n")
  }
  cat("    Clean data has", nrow(all_data_clean), "rows after removing all NAs\n")
  cat("    Clean data annotation distribution:\n")
  print(table(all_data_clean$annotation))
  
  # ===== 保存原始数据用于虚线计算 =====
  # 我们需要保存length列用于后续计算500bp和1000bp的归一化位置
  length_raw_data <- all_data_clean$length
  
  # ===== WINSORIZATION to handle extreme values =====
  winsorize <- function(x, lower_percentile = 0.01, upper_percentile = 0.99) {
    if (all(is.na(x))) return(x)
    lower_bound <- quantile(x, probs = lower_percentile, na.rm = TRUE)
    upper_bound <- quantile(x, probs = upper_percentile, na.rm = TRUE)
    x[x < lower_bound] <- lower_bound
    x[x > upper_bound] <- upper_bound
    return(x)
  }
  
  # Apply winsorization to key variables
  all_data_winsorized <- all_data_clean %>%
    mutate(
      length_win = winsorize(length),
      case_depth_win = winsorize(case.mean.depth),
      control_depth_win = winsorize(control.mean.depth)
    )
  
  # ===== 计算每个项目的归一化参数（用于虚线位置计算） =====
  # 计算length的log10转换和归一化参数
  length_log_values <- log10(length_raw_data + 1)
  length_min <- min(length_log_values, na.rm = TRUE)
  length_max <- max(length_log_values, na.rm = TRUE)
  
  # 计算500bp和1000bp对应的归一化值
  calculate_normalized_position <- function(bp_value) {
    # 1. 对数转换
    bp_log <- log10(bp_value + 1)
    # 2. 归一化到[0,1]
    if (length_max == length_min) {
      return(0.5)
    } else {
      return((bp_log - length_min) / (length_max - length_min))
    }
  }
  
  # 计算500bp和1000bp的归一化位置
  position_500bp <- calculate_normalized_position(500)
  position_1000bp <- calculate_normalized_position(1000)
  
  cat("    Normalization parameters for length:\n")
  cat("      Min(log10):", round(length_min, 4), "\n")
  cat("      Max(log10):", round(length_max, 4), "\n")
  cat("      500bp position (normalized):", round(position_500bp, 4), "\n")
  cat("      1000bp position (normalized):", round(position_1000bp, 4), "\n")
  
  # ===== POP FILTERING (control.F is already non-NA due to early removal) =====
  pop_filtered_data <- all_data_clean %>%
    filter(case.F >= 2 & control.F <= 1)
  cat("    POP filtered data (case.F >=2 & control.F <=1):", 
      nrow(pop_filtered_data), "rows.\n")
  
  # ===== LOG TRANSFORMATION AFTER WINSORIZATION =====
  log_transform <- function(data, use_winsorized = TRUE) {
    if (nrow(data) == 0) return(data.frame())
    if (use_winsorized) {
      data %>%
        mutate(
          length_log = log10(length_win + 1),
          case_depth_log = log10(case_depth_win + 1),
          control_depth_log = log10(control_depth_win + 1)
        )
    } else {
      data %>%
        mutate(
          length_log = log10(length + 1),
          case_depth_log = log10(case.mean.depth + 1),
          control_depth_log = log10(control.mean.depth + 1)
        )
    }
  }
  
  all_data_log <- log_transform(all_data_winsorized, TRUE)
  
  if (nrow(pop_filtered_data) > 0) {
    pop_winsorized <- pop_filtered_data %>%
      mutate(
        length_win = winsorize(length),
        case_depth_win = winsorize(case.mean.depth),
        control_depth_win = winsorize(control.mean.depth)
      )
    pop_data_log <- log_transform(pop_winsorized, TRUE)
  } else {
    pop_data_log <- all_data_log[0, ]
  }
  
  # ===== NORMALIZATION =====
  normalize_to_01 <- function(x) {
    x_no_na <- x[!is.na(x)]
    if (length(x_no_na) == 0) {
      return(rep(NA, length(x)))
    } else if (max(x_no_na) == min(x_no_na)) {
      return(rep(0.5, length(x)))
    } else {
      return((x - min(x_no_na)) / (max(x_no_na) - min(x_no_na)))
    }
  }
  
  prepare_plot_data <- function(data, dataset_type, project_id) {
    if (nrow(data) == 0) return(data.frame())
    data %>%
      mutate(
        length_scaled = normalize_to_01(length_log),
        case_depth_scaled = normalize_to_01(case_depth_log),
        control_depth_scaled = normalize_to_01(control_depth_log)
      ) %>%
      select(annotation, length_scaled, case_depth_scaled, control_depth_scaled) %>%
      pivot_longer(
        cols = c(length_scaled, case_depth_scaled, control_depth_scaled),
        names_to = "variable",
        values_to = "value_scaled"
      ) %>%
      filter(!is.na(value_scaled)) %>%
      mutate(
        variable = factor(variable,
          levels = c("length_scaled", "case_depth_scaled", "control_depth_scaled"),
          labels = c("Sequence Length", "Case Mean Depth", "Control Mean Depth")
        ),
        dataset = dataset_type,
        project = project_id
      )
  }
  
  all_data_plot <- prepare_plot_data(all_data_log, "All Data", project_name)
  pop_data_plot <- prepare_plot_data(pop_data_log, "POP Filtered", project_name)
  
  combined_data <- bind_rows(all_data_plot, pop_data_plot)
  
  # ===== 为每个项目存储虚线位置 =====
  dash_positions <- list(
    project = project_name,
    position_500bp = position_500bp,
    position_1000bp = position_1000bp
  )
  
  # ===== WILCOXON RANK SUM TEST =====
  calculate_wilcox_test <- function(data, var_name) {
    if (nrow(data) == 0 || length(unique(data$annotation)) < 2) {
      return(NA)
    }
    # Extract the two groups
    group_y <- data[[var_name]][data$annotation == "Y"]
    group_n <- data[[var_name]][data$annotation == "N"]
    # Check if both groups have data
    if (length(group_y) == 0 || length(group_n) == 0) {
      return(NA)
    }
    # Perform Wilcoxon rank sum test
    test_result <- tryCatch({
      wilcox.test(group_y, group_n, exact = FALSE, na.action = na.omit)
    }, error = function(e) {
      return(NULL)
    })
    if (is.null(test_result)) return(NA)
    return(test_result$p.value)
  }
  
  # Calculate p-values using Wilcoxon rank sum test
  p_values <- list()
  # For All Data
  if (nrow(all_data_log) > 0) {
    p_values[[paste0(project_name, "_All_length")]] <- calculate_wilcox_test(all_data_log, "length_log")
    p_values[[paste0(project_name, "_All_case")]] <- calculate_wilcox_test(all_data_log, "case_depth_log")
    p_values[[paste0(project_name, "_All_control")]] <- calculate_wilcox_test(all_data_log, "control_depth_log")
  }
  # For POP Filtered
  if (nrow(pop_data_log) > 0) {
    p_values[[paste0(project_name, "_POP_length")]] <- calculate_wilcox_test(pop_data_log, "length_log")
    p_values[[paste0(project_name, "_POP_case")]] <- calculate_wilcox_test(pop_data_log, "case_depth_log")
    p_values[[paste0(project_name, "_POP_control")]] <- calculate_wilcox_test(pop_data_log, "control_depth_log")
  }
  
  # Return results
  return(list(
    data = combined_data,
    p_values = p_values,
    dash_positions = dash_positions,
    counts = c(
      all_data = nrow(all_data_clean),
      pop_data = nrow(pop_filtered_data),
      all_y = sum(all_data_clean$annotation == "Y"),
      all_n = sum(all_data_clean$annotation == "N"),
      pop_y = ifelse(nrow(pop_filtered_data) > 0, sum(pop_filtered_data$annotation == "Y"), 0),
      pop_n = ifelse(nrow(pop_filtered_data) > 0, sum(pop_filtered_data$annotation == "N"), 0),
      removed_na = na_removed + other_na_removed
    )
  ))
}

# ==================== 3. Process ALL projects and combine data ====================
all_results <- list()
all_plot_data <- data.frame()
all_p_values <- list()
all_dash_positions <- list()
project_counts <- data.frame()

cat("\n>>> Starting batch processing of all projects...\n")

for (i in 1:nrow(config)) {
  current_folder <- config$path[i]
  current_name <- config$prj[i]
  result <- process_project_data(current_folder, current_name)
  if (!is.null(result)) {
    all_results[[current_name]] <- result
    if (nrow(result$data) > 0) {
      all_plot_data <- bind_rows(all_plot_data, result$data)
    }
    all_p_values <- c(all_p_values, result$p_values)
    all_dash_positions[[current_name]] <- result$dash_positions
    project_counts <- bind_rows(project_counts, data.frame(
      project = current_name,
      removed_na = result$counts["removed_na"],
      all_total = result$counts["all_data"],
      all_y = result$counts["all_y"],
      all_n = result$counts["all_n"],
      pop_total = result$counts["pop_data"],
      pop_y = result$counts["pop_y"],
      pop_n = result$counts["pop_n"],
      stringsAsFactors = FALSE
    ))
  }
}

cat("\n>>> All projects processed successfully!\n")
cat("Total projects with data:", length(all_results), "\n")
cat("Total rows in combined data:", nrow(all_plot_data), "\n")

# ==================== 4. Prepare p-value labels ====================
p_label_data <- expand.grid(
  project = unique(all_plot_data$project),
  variable = factor(c("Sequence Length", "Case Mean Depth", "Control Mean Depth"),
                    levels = c("Sequence Length", "Case Mean Depth", "Control Mean Depth")),
  dataset = factor(c("All Data", "POP Filtered"), levels = c("All Data", "POP Filtered")),
  stringsAsFactors = FALSE
)

get_p_value <- function(project, variable, dataset) {
  var_map <- c(
    "Sequence Length" = "length",
    "Case Mean Depth" = "case",
    "Control Mean Depth" = "control"
  )
  key <- paste0(project, "_", ifelse(dataset == "All Data", "All", "POP"), "_", var_map[variable])
  if (key %in% names(all_p_values)) {
    return(all_p_values[[key]])
  }
  return(NA)
}

p_label_data <- p_label_data %>%
  rowwise() %>%
  mutate(
    p_value = get_p_value(project, variable, dataset),
    p_label = ifelse(is.na(p_value), "N/A",
                     ifelse(p_value < 0.001, "p < 0.001",
                            paste0("p = ", format(p_value, digits = 3))))
  ) %>%
  ungroup()

# Filter out combinations that don't exist
plot_projects <- unique(all_plot_data$project)
p_label_data <- p_label_data %>%
  filter(project %in% plot_projects)

# ==================== 5. 准备虚线数据 ====================
# 创建包含每个项目虚线位置的数据框
dash_lines_data <- data.frame()
for (proj_name in names(all_dash_positions)) {
  dash_info <- all_dash_positions[[proj_name]]
  dash_lines_data <- bind_rows(dash_lines_data, data.frame(
    project = proj_name,
    position = c(dash_info$position_500bp, dash_info$position_1000bp),
    label = c("500 bp", "1000 bp"),
    stringsAsFactors = FALSE
  ))
}

# ==================== 6. Create the IMPROVED COMBINED plot ====================
cat("\n>>> Creating improved combined visualization...\n")

if (nrow(all_plot_data) == 0) {
  stop("No data available for plotting. Please check your data files.")
}

# 创建基础图形（移除副标题和图例）
combined_plot <- ggplot(all_plot_data, aes(x = annotation, y = value_scaled)) +
  # 添加抖动点
  geom_jitter(aes(color = annotation), width = 0.25, height = 0, alpha = 0.15, size = 0.5) +
  # 添加箱线图
  geom_boxplot(aes(fill = annotation), alpha = 0.7, outlier.shape = NA, width = 0.6) +
  # 添加p值标签
  geom_text(data = p_label_data,
            aes(x = 1.5, y = 1.05, label = p_label),
            inherit.aes = FALSE, size = 2.8, vjust = 0) +
  # 分面显示
  facet_grid(project ~ dataset + variable, scales = "free_y") +
  # NEJM颜色方案
  scale_fill_nejm() +
  scale_color_nejm() +
  # 标题和标签（移除副标题，X轴标签说明Y和N的含义）
  labs(
    title = "Cross-Project Analysis: Distribution by Pathogen Annotation",
    x = "Annotation (Y=Responsible Pathogen, N=Not Responsible)",
    y = "Normalized Value (winsorized + log10 transform)"
  ) +
  # 主题设置（移除图例）
  theme_minimal(base_size = 9) +
  theme(
    # 完全移除图例
    legend.position = "none",
    # 分面标签设置
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold", size = 7),
    strip.text.x = element_text(angle = 0, margin = margin(t = 2, b = 2)),
    strip.text.y = element_text(angle = 0, margin = margin(r = 2, l = 2)),
    # 标题设置
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 8)),
    # 坐标轴设置
    axis.text.x = element_text(size = 8, angle = 0),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 9),
    # 面板设置
    panel.spacing = unit(0.5, "lines"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.3),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  # Y轴设置
  scale_y_continuous(
    limits = c(0, 1.1),
    expand = c(0, 0),
    breaks = seq(0, 1, by = 0.2),
    labels = label_number(accuracy = 0.1)
  )

# ==================== 7. 为每个项目的Sequence Length添加虚线 ====================
# 需要为每个项目在每个数据集（All Data和POP Filtered）的Sequence Length分面中添加虚线
# 这需要更复杂的处理，因为ggplot2的facet_grid中每个分面是独立的

# 方法：为每个项目创建自定义的几何对象
# 首先，我们需要过滤出Sequence Length的数据
sequence_length_data <- all_plot_data %>%
  filter(variable == "Sequence Length")

# 为每个项目创建虚线
for (proj_name in unique(sequence_length_data$project)) {
  # 获取当前项目的虚线位置
  if (proj_name %in% names(all_dash_positions)) {
    dash_info <- all_dash_positions[[proj_name]]
    
    # 为All Data和POP Filtered都添加虚线
    for (dataset_type in c("All Data", "POP Filtered")) {
      # 检查该数据集是否有数据
      if (any(sequence_length_data$project == proj_name & sequence_length_data$dataset == dataset_type)) {
        # 添加500bp虚线（红色）
        combined_plot <- combined_plot +
          geom_hline(
            data = data.frame(
              project = proj_name,
              dataset = dataset_type,
              variable = "Sequence Length",
              yintercept = dash_info$position_500bp
            ),
            aes(yintercept = yintercept),
            linetype = "dashed",
            color = "red",
            alpha = 0.7,
            linewidth = 0.5
          )
        
        # 添加1000bp虚线（蓝色）
        combined_plot <- combined_plot +
          geom_hline(
            data = data.frame(
              project = proj_name,
              dataset = dataset_type,
              variable = "Sequence Length",
              yintercept = dash_info$position_1000bp
            ),
            aes(yintercept = yintercept),
            linetype = "dashed",
            color = "blue",
            alpha = 0.7,
            linewidth = 0.5
          )
        
        # 添加标签（只在一个数据集上添加，避免重复）
        if (dataset_type == "All Data") {
          # 添加500bp标签
          combined_plot <- combined_plot +
            geom_text(
              data = data.frame(
                project = proj_name,
                dataset = dataset_type,
                variable = "Sequence Length",
                x = 1.5,
                y = dash_info$position_500bp + 0.03,
                label = "500 bp"
              ),
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              size = 2.5,
              color = "red",
              vjust = 0
            )
          
          # 添加1000bp标签
          combined_plot <- combined_plot +
            geom_text(
              data = data.frame(
                project = proj_name,
                dataset = dataset_type,
                variable = "Sequence Length",
                x = 1.5,
                y = dash_info$position_1000bp + 0.03,
                label = "1000 bp"
              ),
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              size = 2.5,
              color = "blue",
              vjust = 0
            )
        }
      }
    }
  }
}

# ==================== 8. Save the combined plot ====================
cat("\n>>> Saving improved combined visualization...\n")

n_projects <- length(unique(all_plot_data$project))
plot_width <- 18
plot_height <- max(6, n_projects * 1.8)

# Save as PNG
ggsave(
  filename = "F.S4.ALL_Projects_Improved_Analysis.png",
  plot = combined_plot,
  width = plot_width,
  height = plot_height,
  dpi = 300,
  limitsize = FALSE
)
cat("  Improved plot saved: F.S4.ALL_Projects_Improved_Analysis.png\n")
cat("  Dimensions:", plot_width, "x", plot_height, "inches\n")

# Save as PDF
pdf("F.S4.ALL_Projects_Improved_Analysis.pdf", width = plot_width, height = plot_height)
print(combined_plot)
dev.off()
cat("  PDF report saved: F.S4.ALL_Projects_Improved_Analysis.pdf\n")

# ==================== 9. Generate summary statistics ====================
cat("\n>>> Generating summary statistics...\n")

# Show NA removal summary
cat("\nNA Removal Summary:\n")
cat("Total rows removed due to NA values:", sum(project_counts$removed_na, na.rm = TRUE), "\n")

# Calculate statistics
project_counts <- project_counts %>%
  mutate(
    all_y_percent = round(all_y / all_total * 100, 1),
    pop_y_percent = ifelse(pop_total > 0, round(pop_y / pop_total * 100, 1), NA)
  )

cat("\nProject Statistics:\n")
print(project_counts)

# Calculate overall statistics
cat("\nOverall Statistics:\n")
cat("Total projects:", nrow(project_counts), "\n")
cat("Total sequences (All Data):", sum(project_counts$all_total), "\n")
cat("Total sequences (POP Filtered):", sum(project_counts$pop_total, na.rm = TRUE), "\n")
cat("Mean % Y annotations (All Data):", round(mean(project_counts$all_y_percent, na.rm = TRUE), 1), "%\n")
cat("Mean % Y annotations (POP Filtered):", 
    round(mean(project_counts$pop_y_percent[!is.na(project_counts$pop_y_percent)], na.rm = TRUE), 1), "%\n")

# 显示虚线位置信息
cat("\nReference Line Positions (Normalized):\n")
for (proj_name in names(all_dash_positions)) {
  dash_info <- all_dash_positions[[proj_name]]
  cat("  ", proj_name, ": 500bp =", round(dash_info$position_500bp, 4), 
      ", 1000bp =", round(dash_info$position_1000bp, 4), "\n")
}

# Save summary table
write.csv(project_counts, "Project_Summary_Improved.csv", row.names = FALSE)

# Count significant results
sig_results <- p_label_data %>%
  filter(!is.na(p_value) & p_value < 0.05)

cat("\nStatistical Test Results:\n")
cat("Total comparisons:", nrow(p_label_data %>% filter(!is.na(p_value))), "\n")
cat("Significant results (p < 0.05):", nrow(sig_results), "\n")
cat("Significant at p < 0.01:", nrow(sig_results %>% filter(p_value < 0.01)), "\n")
cat("Significant at p < 0.001:", nrow(sig_results %>% filter(p_value < 0.001)), "\n")

cat("\n✅ Analysis completed with all improvements!\n")
cat("• Removed subtitle for cleaner presentation\n")
cat("• Removed legend, added explanation to X-axis label\n")
cat("• Added dashed reference lines at 500bp and 1000bp for Sequence Length\n")
cat("• Generated improved visualization: ALL_Projects_Improved_Analysis.png/pdf\n")
cat("• Saved summary statistics: Project_Summary_Improved.csv\n")
