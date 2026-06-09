# 加载必要的包
library(ggsci)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(cowplot)

# 定义通用颜色
colorsQ <- c("grey50", pal_nejm("default")(1))
names(colorsQ) <- c("N", "Y")

# 辅助函数：获取非空颜色
get_colors <- function(data) {
  cols <- colorsQ[names(colorsQ) %in% unique(data$annotation)]
  return(cols)
}

# ==================== 单病例数据集（PRJ.39）====================
Main1 <- read.table("PRJ.39.sars.V1.expresson.table", head = TRUE, sep = "\t")
Main1$log_meanDepth <- log(Main1$case.mean.depth)

# 图 A：原始数据
nA <- nrow(Main1)
purityA <- round(sum(Main1[Main1$annotation == "Y", "length"]) / sum(Main1$length) * 100, digits = 2)
pA <- ggscatter(Main1, "gc", "log_meanDepth",
                xlab = "GC content", ylab = "log(Length)",
                color = "annotation", alpha = 0.9, legend = "none",
                palette = colorsQ, size = "length",
                title = "A. Single case: all contigs") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.position = "none") +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("n = ", nA, "\nsignal purity = ", purityA, "%"),
           hjust = 1.1, vjust = 1.1, size = 3, color = "gray30")

# 图 B：最终过滤（case.fre ≥1 & control.fre ≤1）
Main1_final <- Main1[Main1$case.fre >= 1 & Main1$control.fre <= 1, ]
nB <- nrow(Main1_final)
purityB <- round(sum(Main1_final[Main1_final$annotation == "Y", "length"]) / sum(Main1_final$length) * 100, digits = 2)
colorsB <- get_colors(Main1_final)
pB <- ggscatter(Main1_final, "gc", "log_meanDepth",
                xlab = "GC content", ylab = "log(Length)",
                color = "annotation", alpha = 0.9, legend = "none",
                palette = colorsB, size = "length",
                title = "B. Single case: case.fre ≥1 & control.fre ≤1") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.position = "none") +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("n = ", nB, "\nsignal purity = ", purityB, "%"),
           hjust = 1.1, vjust = 1.1, size = 3, color = "gray30")

# ==================== 两病例数据集（PRJ.40）====================
Main2 <- read.table("PRJ.40.ljc_sar2_2019.V1.expression.table", head = TRUE, sep = "\t")
Main2$log_meanDepth <- log(Main2$case.mean.depth)

# 图 C：原始数据
nC <- nrow(Main2)
purityC <- round(sum(Main2[Main2$annotation == "Y", "length"]) / sum(Main2$length) * 100, digits = 2)
pC <- ggscatter(Main2, "gc", "log_meanDepth",
                xlab = "GC content", ylab = "log(Length)",
                color = "annotation", alpha = 0.9, legend = "none",
                palette = colorsQ, size = "length",
                title = "C. Two cases: all contigs") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.position = "none") +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("n = ", nC, "\nsignal purity = ", purityC, "%"),
           hjust = 1.1, vjust = 1.1, size = 3, color = "gray30")

# 图 D：case.fre ≥2
Main2_m2 <- Main2[Main2$case.fre >= 2, ]
nD <- nrow(Main2_m2)
purityD <- round(sum(Main2_m2[Main2_m2$annotation == "Y", "length"]) / sum(Main2_m2$length) * 100, digits = 2)
colorsD <- get_colors(Main2_m2)
pD <- ggscatter(Main2_m2, "gc", "log_meanDepth",
                xlab = "GC content", ylab = "log(Length)",
                color = "annotation", alpha = 0.9, legend = "none",
                palette = colorsD, size = "length",
                title = "D. Two cases: case.fre ≥2") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.position = "none") +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("n = ", nD, "\nsignal purity = ", purityD, "%"),
           hjust = 1.1, vjust = 1.1, size = 3, color = "gray30")

# 图 E：最终过滤（case.fre ≥2 & control.fre ≤1）
Main2_final <- Main2_m2[Main2_m2$control.fre <= 1, ]
nE <- nrow(Main2_final)
purityE <- round(sum(Main2_final[Main2_final$annotation == "Y", "length"]) / sum(Main2_final$length) * 100, digits = 2)
colorsE <- get_colors(Main2_final)
pE <- ggscatter(Main2_final, "gc", "log_meanDepth",
                xlab = "GC content", ylab = "log(Length)",
                color = "annotation", alpha = 0.9, legend = "none",
                palette = colorsE, size = "length",
                title = "E. Two cases: case.fre ≥2 & control.fre ≤1") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.position = "none") +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("n = ", nE, "\nsignal purity = ", purityE, "%"),
           hjust = 1.1, vjust = 1.1, size = 3, color = "gray30")

# ==================== 构建共享图例（颜色 + 大小说明） ====================
# 颜色图例
legend_color <- ggplot() +
  geom_point(aes(x = 0, y = 0, color = "N"), size = 2) +
  geom_point(aes(x = 0, y = 0, color = "Y"), size = 2) +
  scale_color_manual(name = "Annotation",
                     values = colorsQ,
                     labels = c("N" = "Not responsible", "Y" = "Responsible")) +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8))
legend_color <- get_legend(legend_color)

# 将颜色图例和大小说明文本放在一起
legend_with_size <- plot_grid(
  legend_color,
  ggdraw() + draw_label("Point size: contig length", size = 8, fontface = "italic"),
  ncol = 1, rel_heights = c(1, 0.4)
)

# ==================== 组合图形 ====================
top_row <- plot_grid(pA, pB, legend_with_size,
                     ncol = 3, nrow = 1,
                     rel_widths = c(1, 1, 0.5),
                     labels = c("", "", ""))

bottom_row <- plot_grid(pC, pD, pE,
                        ncol = 3, nrow = 1,
                        rel_widths = c(1, 1, 1),
                        labels = c("", "", ""))

combined <- plot_grid(top_row, bottom_row, ncol = 1, nrow = 2, rel_heights = c(1, 1))

ggsave("F_S9.png", plot = combined, width = 14, height = 10, dpi = 300, bg = "white")
ggsave("F_S9.pdf", plot = combined, width = 14, height = 10, dpi = 300, bg = "white")

