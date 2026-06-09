library(ggsci)
library(ggpubr)
library(ggplot2)

# ====================== 1. 读取数据 ======================
data <- read.table("PRJ.80.congo.f5.table.txt",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

if ("contro.mean.cov" %in% colnames(data)) {
  colnames(data)[colnames(data) == "contro.mean.cov"] <- "control.mean.cov"
}

# ====================== 2. 数据预处理 ======================
data$annotation_clean <- data$annotation
data$annotation_clean[grepl("Plasmodium falciparum 3D7", data$annotation, ignore.case = TRUE)] <- "Plasmodium falciparum"

data$group <- "others"
data$group[grepl("plasmodium", data$annotation_clean, ignore.case = TRUE)] <- "Plasmodium falciparum"
data$group[data$annotation_clean == "Homo sapiens"] <- "Homo sapiens"
data$group[data$annotation_clean == "TTV-like mini virus"] <- "TTV-like mini virus"
data$group[data$annotation_clean == "Torque teno virus"] <- "Torque teno virus"

eps <- 1e-6
data$log10_depth1 <- log10(data$stage1_mean_depth + eps)
data$log10_depth2 <- log10(data$stage2_mean_depth + eps)

# ====================== 3. 颜色映射 ======================
nejm_colors <- pal_nejm("default")(8)
color_mapping <- c(
  "Plasmodium falciparum" = nejm_colors[1],
  "Homo sapiens"          = nejm_colors[2],
  "TTV-like mini virus"   = nejm_colors[3],
  "Torque teno virus"     = nejm_colors[4],
  "others"                = "grey70"
)
color_mapping <- color_mapping[names(color_mapping) %in% unique(data$group)]

# ====================== 4. 过滤步骤 ======================
data_step1 <- subset(data, stage1_fre >= 2)
data_step2 <- subset(data_step1, control.fre <= 1)

# ====================== 5. 辅助函数 ======================
calc_signal_ratio <- function(df, target = "Plasmodium falciparum") {
  total_len <- sum(df$length, na.rm = TRUE)
  target_len <- sum(df$length[df$group == target], na.rm = TRUE)
  round(target_len / total_len * 100, 2)
}

# ====================== 6. 统计信息 ======================
n_all <- nrow(data); sr_all <- calc_signal_ratio(data)
n_step1 <- nrow(data_step1); sr_step1 <- calc_signal_ratio(data_step1)
n_step2 <- nrow(data_step2); sr_step2 <- calc_signal_ratio(data_step2)

pf_data <- subset(data, group == "Plasmodium falciparum")
others_data <- subset(data, group == "others")
r_pf_depth <- if(nrow(pf_data)>1) cor(pf_data$log10_depth1, pf_data$log10_depth2, use="complete.obs") else NA
r_others_depth <- if(nrow(others_data)>1) cor(others_data$log10_depth1, others_data$log10_depth2, use="complete.obs") else NA
r_pf_cov <- if(nrow(pf_data)>1) cor(pf_data$stage1_mean_cov, pf_data$stage2_mean_cov, use="complete.obs") else NA
r_others_cov <- if(nrow(others_data)>1) cor(others_data$stage1_mean_cov, others_data$stage2_mean_cov, use="complete.obs") else NA

# ====================== 7. 基础主题 ======================
base_theme <- theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )

# ====================== 8. 图A、B、C ======================
# 图A（保留图例用于底部共享）
p_all <- ggscatter(data, x = "gc", y = "stage1_mean_depth",
                   color = "group", palette = color_mapping, size = "length",
                   xlab = "GC content", ylab = "Stage1 Mean Depth (log scale)",
                   title = "All contigs", alpha = 0.9) +
  scale_y_log10() +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("n = ", n_all, "\nSR = ", sr_all, "%"),
           size = 3.5, fontface = "bold") +
  base_theme +
  theme(legend.position = "right")

# 图B
p_step1 <- ggscatter(data_step1, x = "gc", y = "stage1_mean_depth",
                     color = "group", palette = color_mapping, size = "length",
                     xlab = "GC content", ylab = "Stage1 Mean Depth (log scale)",
                     title = "Stage1_fre ≥ 2", alpha = 0.9) +
  scale_y_log10() +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("n = ", n_step1, "\nSR = ", sr_step1, "%"),
           size = 3.5, fontface = "bold") +
  base_theme

# 图C：注释放在右侧中间位置，避免挤压
x_pos_c <- max(data_step2$gc, na.rm = TRUE) * 0.98
y_pos_c <- median(data_step2$stage1_mean_depth, na.rm = TRUE)
p_step2 <- ggscatter(data_step2, x = "gc", y = "stage1_mean_depth",
                     color = "group", palette = color_mapping, size = "length",
                     xlab = "GC content", ylab = "Stage1 Mean Depth (log scale)",
                     title = "Control.fre ≤ 1", alpha = 0.9) +
  scale_y_log10() +
  annotate("text", x = x_pos_c, y = y_pos_c, hjust = 1, vjust = 0.5,
           label = paste0("n = ", n_step2, "\nSR = ", sr_step2, "%"),
           size = 3.5, fontface = "bold") +
  base_theme

# ====================== 9. 图D、E、F ======================
current_colors <- color_mapping[names(color_mapping) %in% unique(data_step1$group)]

p_comp_cov <- ggplot(data_step1, aes(x = stage1_mean_cov, y = control.mean.cov, color = group)) +
  geom_point(size = 1.5, alpha = 0.8) + scale_color_manual(values = current_colors) +
  labs(x = "Stage1 Mean Coverage", y = "Controls Mean Coverage") +
  ggtitle("Coverage Comparison") + base_theme

p_comp_freq <- ggplot(data_step1, aes(x = stage1_fre, y = control.fre, color = group)) +
  geom_point(size = 1.5, alpha = 0.8) + scale_color_manual(values = current_colors) +
  labs(x = "Stage1 Frequency", y = "Controls Frequency") +
  ggtitle("Frequency Comparison") + base_theme

p_comp_depth <- ggplot(data_step1, aes(x = stage1_mean_depth, y = control.mean.depth, color = group)) +
  geom_point(size = 1.5, alpha = 0.8) + scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = current_colors) +
  labs(x = "Stage1 Mean Depth (log scale)", y = "Controls Mean Depth (log scale)") +
  ggtitle("Depth Comparison") + base_theme

# ====================== 10. 图G、H、I ======================
all_colors <- color_mapping[names(color_mapping) %in% unique(data$group)]

p_freq_stage <- ggplot(data, aes(x = stage1_fre, y = stage2_fre, color = group)) +
  geom_point(size = 1.5, alpha = 0.8) + scale_color_manual(values = all_colors) +
  labs(x = "Stage1 Frequency", y = "Stage2 Frequency") +
  ggtitle("Stage1 vs Stage2 Frequency") + base_theme

# 深度比较 + 相关系数
x_pos_depth <- min(data$stage1_mean_depth, na.rm = TRUE) * 1.2
y_pos_depth <- max(data$stage2_mean_depth, na.rm = TRUE) * 0.9
p_depth_stage <- ggplot(data, aes(x = stage1_mean_depth, y = stage2_mean_depth, color = group)) +
  geom_point(size = 1.5, alpha = 0.8) + scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = all_colors) +
  labs(x = "Stage1 Mean Depth (log scale)", y = "Stage2 Mean Depth (log scale)") +
  ggtitle("Stage1 vs Stage2 Depth") + base_theme +
  annotate("text", x = x_pos_depth, y = y_pos_depth,
           label = sprintf("Pf: r = %.2f\nOthers: r = %.2f", r_pf_depth, r_others_depth),
           color = "black", size = 3.5, fontface = "bold", hjust = 0, vjust = 1)

# 覆盖度比较 + 相关系数
x_pos_cov <- min(data$stage1_mean_cov, na.rm = TRUE) * 1.2
y_pos_cov <- max(data$stage2_mean_cov, na.rm = TRUE) * 0.9
p_cov_stage <- ggplot(data, aes(x = stage1_mean_cov, y = stage2_mean_cov, color = group)) +
  geom_point(size = 1.5, alpha = 0.8) + scale_color_manual(values = all_colors) +
  labs(x = "Stage1 Mean Coverage", y = "Stage2 Mean Coverage") +
  ggtitle("Stage1 vs Stage2 Coverage") + base_theme +
  annotate("text", x = x_pos_cov, y = y_pos_cov,
           label = sprintf("Pf: r = %.2f\nOthers: r = %.2f", r_pf_cov, r_others_cov),
           color = "black", size = 3.5, fontface = "bold", hjust = 0, vjust = 1)

# ====================== 11. 组合图形 ======================
combined_plot <- ggarrange(
  p_all, p_step1, p_step2,
  p_comp_cov, p_comp_freq, p_comp_depth,
  p_freq_stage, p_depth_stage, p_cov_stage,
  nrow = 3, ncol = 3,
  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
  font.label = list(size = 14, face = "bold"),
  common.legend = TRUE,
  legend = "bottom"
)

# ====================== 12. 保存 ======================
ggsave("S_F12.png", plot = combined_plot, width = 18, height = 15, dpi = 300, bg = "white")
ggsave("S_F12.pdf", plot = combined_plot, width = 18, height = 15, dpi = 300, bg = "white")

