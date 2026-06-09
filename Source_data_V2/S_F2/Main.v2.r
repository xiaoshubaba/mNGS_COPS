library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# ================================
# Jaccard index across projects (sorted by mean)
# Simulated outbreak cluster size fixed at 5 cases
# ================================
jd_df <- read_tsv("PRJ.rand.seeds.status.txt",
                  col_types = cols(
                    prj = col_character(),
                    outbreak1 = col_character(),
                    outbreak2 = col_character(),
                    intersection = col_double(),
                    symmetric_difference = col_double(),
                    jaccard_index = col_double()
                  ))

# 计算每个项目的平均 Jaccard index
prj_means <- jd_df %>%
  group_by(prj) %>%
  summarise(mean_jd = mean(jaccard_index, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_jd))

# 按均值排序因子水平
jd_df$prj <- factor(jd_df$prj, levels = prj_means$prj)

# 颜色设置
col_main <- "#4E79A7"
col_fill <- alpha(col_main, 0.15)
col_points <- alpha(col_main, 0.6)

# 绘图（标题已移除 "D."，添加副标题）
p <- ggplot(jd_df, aes(x = prj, y = jaccard_index)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey85", linewidth = 0.4) +
  geom_violin(fill = col_fill, color = NA, width = 0.8, alpha = 0.7) +
  geom_jitter(width = 0.2, height = 0, size = 1.2,
              color = col_points, alpha = 0.5) +
  geom_point(data = prj_means, aes(x = prj, y = mean_jd),
             color = "grey25", fill = "white",
             shape = 21, size = 2.5, stroke = 1) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.6, size = 0.4, color = "grey40") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Project similarity distribution (Jaccard index)",
    subtitle = "Simulated outbreak cluster size fixed at 5 cases",
    x = "Project",
    y = "Jaccard index",
    caption = paste0("N = ", nrow(jd_df), " pairwise comparisons\n",
                    "Points: individual comparisons | White circles: mean values")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "plain"),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey93", linewidth = 0.3),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks.x = element_line(color = "grey80"),
    plot.margin = margin(5, 15, 5, 15),
    plot.caption = element_text(size = 8, color = "gray40", hjust = 0),
    legend.position = "none"
  )

# 保存为 PDF 和 PNG，文件名为 F_S2.jaccard_distribution
ggsave("F_S2.jaccard_distribution.pdf", plot = p, width = 10, height = 6, dpi = 300)
ggsave("F_S2.jaccard_distribution.png", plot = p, width = 10, height = 6, dpi = 300)

# 显示图形（可选）

# 统计摘要
cat("\n=== Jaccard Index Analysis (Fixed cluster size = 5) ===\n")
cat(sprintf("Total comparisons: %d\n", nrow(jd_df)))
cat(sprintf("Number of projects: %d\n", length(unique(jd_df$prj))))
cat(sprintf("Overall mean Jaccard index: %.3f\n", mean(jd_df$jaccard_index, na.rm = TRUE)))
cat(sprintf("Overall median Jaccard index: %.3f\n", median(jd_df$jaccard_index, na.rm = TRUE)))
cat("\nProjects sorted by mean Jaccard index (descending):\n")
print(prj_means)
