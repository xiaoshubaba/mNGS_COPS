library(ggsci)
library(ggpubr)
library(ggplot2)

# 读取数据
Main <- read.table("PRJ.15YFL.expression.V3.table", header = TRUE, sep = "\t")
Main$log_meanDepth <- log(Main$case.mean.depth)

# 颜色定义
colorsQ <- c("grey50", pal_nejm("default")(1))

# ---- 图 A：原始数据 ----
sr <- round(sum(Main[Main$annotation == "Y", "length"]) / sum(Main$length) * 100, digits = 2)
f1_left <- ggscatter(Main, "gc", "log_meanDepth",
                     xlab = "GC content", ylab = "log(Mean depth in cases)",
                     color = "annotation", alpha = 0.9, legend = "none",
                     palette = colorsQ, size = "length",
                     title = paste0("Contigs: ", nrow(Main), ", signal purity: ", sr, "%")) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))

# ---- 图 B：case.fre ≥ 2 ----
Main_m2 <- Main[Main$case.fre >= 2, ]
sr <- round(sum(Main_m2[Main_m2$annotation == "Y", "length"]) / sum(Main_m2$length) * 100, digits = 2)
colorsQ_tem <- colorsQ[table(Main_m2$annotation) > 0]
f1_middle <- ggscatter(Main_m2, "gc", "log_meanDepth",
                       xlab = "GC content", ylab = "log(Mean depth in cases)",
                       color = "annotation", alpha = 0.9, legend = "none",
                       palette = colorsQ_tem, size = "length",
                       title = paste0("Contigs: ", nrow(Main_m2), ", signal purity: ", sr, "%")) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))

# ---- 图 C：case.fre ≥ 2 且 control.fre ≤ 1 ----
Main_m2_m2 <- Main_m2[Main_m2$control.fre <= 1, ]
colorsQ_tem <- colorsQ[table(Main_m2_m2$annotation) > 0]
sr <- round(sum(Main_m2_m2[Main_m2_m2$annotation == "Y", "length"]) / sum(Main_m2_m2$length) * 100, digits = 2)
f1_right <- ggscatter(Main_m2_m2, "gc", "log_meanDepth",
                      xlab = "GC content", ylab = "log(Mean depth in cases)",
                      color = "annotation", alpha = 0.9, legend = "none",
                      palette = colorsQ_tem, size = "length",
                      title = paste0("Contigs: ", nrow(Main_m2_m2), ", signal purity: ", sr, "%")) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))

# ---- 图 D：病例 vs 对照 Coverage（使用 Main_m2 数据） ----
Main <- Main_m2
colorsQ_tem <- colorsQ[table(Main$annotation) > 0]
p <- ggplot(Main, aes(x = case.mean.cov, y = control.mean.cov, color = annotation)) +
  geom_point(size = 1.5, alpha = 0.9) +
  scale_color_manual(values = colorsQ_tem) +
  labs(x = "Patients Mean Coverage", y = "Controls Mean Coverage") +
  theme_bw() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.position = "none")

# ---- 图 E：病例 vs 对照 Frequency ----
p1 <- ggplot(Main, aes(x = case.fre, y = control.fre, color = annotation)) +
  geom_point(size = 1.5, alpha = 0.9) +
  scale_color_manual(values = colorsQ_tem) +
  labs(x = "Patients Frequency", y = "Controls Frequency") +
  theme_bw() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.position = "none")

# ---- 图 F：病例 vs 对照 Depth (log) ----
p2 <- ggplot(Main, aes(x = log(case.mean.depth), y = log(control.mean.depth), color = annotation)) +
  geom_point(size = 1.5, alpha = 0.9) +
  scale_color_manual(values = colorsQ_tem) +
  labs(x = "log(Patients Mean Depth)", y = "log(Controls Mean Depth)") +
  theme_bw() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.position = "none")

# 组合图形（2行3列）
pCyto <- ggarrange(f1_left, f1_middle, f1_right, p, p1, p2,
                   labels = c("A", "B", "C", "D", "E", "F"),
                   ncol = 3, nrow = 2,
                   font.label = list(size = 12, face = "bold"),
                   common.legend = FALSE)

# 保存为 PNG 和 PDF
ggsave("F_S7.png", plot = pCyto, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("F_S7.pdf", plot = pCyto, width = 12, height = 8, dpi = 300, bg = "white")

