library(ggsci)
library(ggpubr)
library(ggplot2)

# 1. 数据读取和预处理
data <- read.table("PRJ.71_complete.annotation.txt",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

data$fragment_clean <- sapply(strsplit(data$fragment, " "), function(x) x[1])
data$annotation_clean <- "Others"
data$annotation_clean[grepl("Rhodococcus", data$annotation, ignore.case = TRUE)] <- "Rhodococcus"
data$annotation_clean[grepl("Influenza A virus", data$annotation, ignore.case = TRUE)] <- "Influenza A virus"
data$annotation_clean[grepl("Severe acute respiratory syndrome coronavirus 2", data$annotation, ignore.case = TRUE)] <- "SARS-CoV-2"

result <- data.frame(
  fragment = data$fragment_clean,
  annotation = data$annotation_clean,
  stringsAsFactors = FALSE
)

# 2. 读取主数据并进行标注
Main <- read.table("PRJ.71.FFPE_IAV.V1.expression.table",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

contig_to_annotation <- setNames(result$annotation, result$fragment)
Main$annotation_new <- Main$annotation

for(i in 1:nrow(Main)) {
  contig_name <- Main$contig[i]
  if(contig_name %in% names(contig_to_annotation)) {
    Main$annotation_new[i] <- contig_to_annotation[contig_name]
  } else if(Main$annotation[i] == "Y") {
    Main$annotation_new[i] <- "Influenza A virus"
  } else {
    Main$annotation_new[i] <- "Others"
  }
}

Main_final <- Main
Main_final$annotation <- Main_final$annotation_new
Main_final$annotation_new <- NULL
Main <- Main_final

# 3. 设置颜色映射
nejm_colors <- pal_nejm("default")(8)
color_mapping <- c(
  "Influenza A virus" = "#BC3C29FF",
  "Rhodococcus"       = "#0072B5FF",
  "SARS-CoV-2"        = nejm_colors[3],
  "Others"            = "grey50"
)

# 4. 第一个散点图（保留图例，供底部共享）
Main$log_meanDepth <- log(Main$case.mean.depth)
n_left <- nrow(Main)
iav_count_left <- sum(Main$annotation == "Influenza A virus")
sr_left <- round(sum(Main[Main$annotation == "Influenza A virus", "length"]) / sum(Main$length) * 100, digits = 2)

f1_left <- ggscatter(Main, "gc", "log_meanDepth",
                     xlab = "GC content",
                     color = "annotation",
                     alpha = 0.9,
                     palette = color_mapping,
                     size = "length",
                     title = "All contigs") +
  annotate("text",
           x = max(Main$gc, na.rm = TRUE) * 0.95,
           y = max(Main$log_meanDepth, na.rm = TRUE) * 0.95,
           label = paste0("n = ", n_left, "\nIAV = ", iav_count_left, "\nSR = ", sr_left, "%"),
           hjust = 1, vjust = 1,
           size = 3.5, color = "black", fontface = "bold") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title = element_text(size = 13),   # 增大
        axis.text = element_text(size = 12))    # 增大

# 5. 第二个散点图（无图例）
Main_m2 <- Main[Main$case.fre >= 2, ]
n_mid <- nrow(Main_m2)
iav_count_mid <- sum(Main_m2$annotation == "Influenza A virus")
sr_mid <- round(sum(Main_m2[Main_m2$annotation == "Influenza A virus", "length"]) / sum(Main_m2$length) * 100, digits = 2)

f1_middle <- ggscatter(Main_m2, "gc", "log_meanDepth",
                       xlab = "GC content",
                       color = "annotation",
                       alpha = 0.9,
                       legend = "none",
                       palette = color_mapping,
                       size = "length",
                       title = "Case freq ≥ 2") +
  annotate("text",
           x = max(Main_m2$gc, na.rm = TRUE) * 0.95,
           y = max(Main_m2$log_meanDepth, na.rm = TRUE) * 0.95,
           label = paste0("n = ", n_mid, "\nIAV = ", iav_count_mid, "\nSR = ", sr_mid, "%"),
           hjust = 1, vjust = 1,
           size = 3.5, color = "black", fontface = "bold") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.position = "none")

# 6. 第三个散点图（无图例）
Main_m2_final <- Main_m2[Main_m2$control.fre <= 1, ]
n_right <- nrow(Main_m2_final)
iav_count_right <- sum(Main_m2_final$annotation == "Influenza A virus")
sr_right <- round(sum(Main_m2_final[Main_m2_final$annotation == "Influenza A virus", "length"]) / sum(Main_m2_final$length) * 100, digits = 2)

f1_right <- ggscatter(Main_m2_final, "gc", "log_meanDepth",
                      xlab = "GC content",
                      color = "annotation",
                      alpha = 0.9,
                      legend = "none",
                      palette = color_mapping,
                      size = "length",
                      title = "Case freq ≥2 & Control freq ≤1") +
  annotate("text",
           x = max(Main_m2_final$gc, na.rm = TRUE) * 0.95,
           y = max(Main_m2_final$log_meanDepth, na.rm = TRUE) * 0.95,
           label = paste0("n = ", n_right, "\nIAV = ", iav_count_right, "\nSR = ", sr_right, "%"),
           hjust = 1, vjust = 1,
           size = 3.5, color = "black", fontface = "bold") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.position = "none")

# 7. 三个比较图（使用 Main_m2 数据）
Main_comp <- Main_m2
colors_comp <- color_mapping[unique(Main_comp$annotation)]

p1 <- ggplot(Main_comp, aes(x = case.mean.cov, y = control.mean.cov, color = annotation)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = colors_comp) +
  labs(x = "Patients Mean Coverage", y = "Controls Mean Coverage") +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  ggtitle("Coverage Comparison")

p2 <- ggplot(Main_comp, aes(x = annotation, y = control.fre, fill = annotation)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.6, aes(color = annotation)) +
  scale_fill_manual(values = colors_comp) +
  scale_color_manual(values = colors_comp) +
  labs(x = "Species", y = "Control Frequency") +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),  # 稍微增大刻度标签
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  ggtitle("Control Frequency by Species")

p3 <- ggplot(Main_comp, aes(x = log(case.mean.depth), y = log(control.mean.depth), color = annotation)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = colors_comp) +
  labs(x = "log(Patients Mean Depth)", y = "log(Controls Mean Depth)") +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  ggtitle("Depth Comparison")

# 8. 组合图形：底部共享图例（从 f1_left 提取）
all_plots_combined <- ggarrange(
  f1_left, f1_middle, f1_right,
  p1, p2, p3,
  nrow = 2, ncol = 3,
  labels = c("A", "B", "C", "D", "E", "F"),
  font.label = list(size = 14, face = "bold"),
  common.legend = TRUE,
  legend = "bottom"
)

# 9. 保存
ggsave("S_F11.png", plot = all_plots_combined, width = 16, height = 10, dpi = 300, bg = "white")
ggsave("S_F11.pdf", plot = all_plots_combined, width = 16, height = 10, dpi = 300, bg = "white")
