# 加载必要的包
library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyr)

# 1. 读取数据
plot_df <- read.table("fw.subSampling.summary.extended.txt", 
                     header = TRUE, 
                     sep = "\t", 
                     stringsAsFactors = FALSE)

# 2. 提取分组信息并设置因子水平
plot_df <- plot_df %>%
  mutate(
    group = gsub("\\..*", "", sample),  # 提取小数点前的部分作为分组
    iteration = as.numeric(gsub(".*\\.iter", "", sample))  # 提取迭代次数
  )

# 设置group的因子水平，固定顺序为R50, R10, R1
plot_df$group <- factor(plot_df$group, levels = c("R50", "R10", "R1"))

# 3. 准备绘图数据
# 左侧面板：单个宏基因组指标
left_data <- plot_df %>%
  select(group, iteration, adjusted_cov, adjust_signal_ratio, Adjusts_reads) %>%
  pivot_longer(
    cols = c(adjusted_cov, adjust_signal_ratio, Adjusts_reads),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    panel = "Single Metagenome",
    metric = factor(metric, 
                    levels = c("adjusted_cov", "adjust_signal_ratio", "Adjusts_reads"),
                    labels = c("Adjusted Coverage (%)", "Adjusted Signal Ratio (%)", "Adjusted Reads"))
  )

# 右侧面板：POPS性能指标
right_data <- plot_df %>%
  select(group, iteration, POPS_pathogen_contigs, POPS_all_contigs, POPS_signal_ratio_per) %>%
  pivot_longer(
    cols = c(POPS_pathogen_contigs, POPS_all_contigs, POPS_signal_ratio_per),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    panel = "POPS Performance",
    metric = factor(metric,
                    levels = c("POPS_pathogen_contigs", "POPS_all_contigs", "POPS_signal_ratio_per"),
                    labels = c("Pathogen Contigs", "All Contigs", "Signal Ratio (%)"))
  )

# 合并数据
plot_data <- bind_rows(left_data, right_data)

# 设置因子水平，确保正确的排序
plot_data$panel <- factor(plot_data$panel, levels = c("Single Metagenome", "POPS Performance"))

# 确保group的顺序保持不变
plot_data$group <- factor(plot_data$group, levels = c("R50", "R10", "R1"))

# 4. 创建使用facet_grid的箱线图
p <- ggplot(plot_data, aes(x = group, y = value, fill = group)) +
  geom_boxplot(
    alpha = 0.8,
    outlier.size = 1.5,
    outlier.color = "gray30",
    width = 0.6
  ) +
  geom_jitter(
    width = 0.15,
    height = 0,
    alpha = 0.5,
    size = 1
  ) +
  scale_fill_nejm() +  # 使用NEJM期刊风格的配色
  facet_grid(metric ~ panel, scales = "free_y") +  # 使用facet_grid分开左右
  labs(
    title = "Performance Comparison Across Sampling Groups",
    x = "Sampling Group",
    y = "Value",
    fill = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "gray90", color = "gray80"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.spacing = unit(1.5, "lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1))

# 5. 显示图形
print(p)

# 6. 输出分组统计摘要
cat("\n分组统计摘要：\n\n")

summary_stats <- plot_data %>%
  group_by(panel, metric, group) %>%
  summarise(
    n = n(),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

# 格式化输出
for(panel_name in unique(summary_stats$panel)) {
  cat(paste0(panel_name, ":\n"))
  panel_data <- summary_stats %>% filter(panel == panel_name)
  
  for(metric_name in unique(panel_data$metric)) {
    cat(paste0("  ", metric_name, ":\n"))
    metric_data <- panel_data %>% filter(metric == metric_name)
    
    for(i in 1:nrow(metric_data)) {
      row <- metric_data[i, ]
      
      # 根据指标类型格式化输出
      if(grepl("%", metric_name)) {
        cat(sprintf("    %s (n=%d): mean=%.1f%%, median=%.1f%%, sd=%.1f, range=[%.1f%%, %.1f%%]\n",
                    row$group, row$n, row$mean, row$median, row$sd, row$min, row$max))
      } else if(metric_name == "Adjusted Reads") {
        cat(sprintf("    %s (n=%d): mean=%.0f, median=%.0f, sd=%.0f, range=[%.0f, %.0f]\n",
                    row$group, row$n, row$mean, row$median, row$sd, row$min, row$max))
      } else {
        cat(sprintf("    %s (n=%d): mean=%.1f, median=%.1f, sd=%.1f, range=[%.1f, %.1f]\n",
                    row$group, row$n, row$mean, row$median, row$sd, row$min, row$max))
      }
    }
  }
  cat("\n")
}
