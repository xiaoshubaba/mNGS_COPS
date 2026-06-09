# 加载必要的包
library(ggplot2)
library(ggpubr)
library(ggsci)

# 颜色设置
colorsQ <- c(pal_nejm("default")(5), "grey50", "grey30")
cR <- rev(pal_nejm("default")(3))

# 读取数据
data <- read.table("Main.table.txt", head = TRUE)
data$round <- factor(data$round, levels = c("First", "Second", "Third"),
                     labels = c("De novo Assemble", "First Filter", "Second Filter"))
data$specy <- factor(data$specy, levels = c("COVID", "IVA"), labels = c("SARS-2-Cov", "IVA"))

MainAB <- read.table("Main.table.top.txt", head = TRUE)
MainAB$IVA_reads_L <- log(MainAB$IVA_reads)

# 图 A
p4 <- ggboxplot(MainAB, x = "cohort", y = "IVA_reads_L", color = "cohort",
                add = "jitter", palette = cR, facet.by = "IVA_dominant_Strain") +
  theme(legend.position = "none", axis.text.x = element_text(size = 8)) +
  labs(x = "", y = "Number of IVA reads (log)") + ggtitle("A")

# 图 B
p5 <- ggboxplot(MainAB, x = "cohort", y = "asmCov", color = "cohort",
                add = "jitter", palette = cR, facet.by = "IVA_dominant_Strain") +
  theme(legend.position = "none", axis.text.x = element_text(size = 8)) +
  labs(x = "", y = "De novo Assembly coverage (%)") + ggtitle("B")

# 图 C
p6 <- ggscatter(MainAB, "IVA_reads_L", "asmCov", facet.by = "cohort",
                color = "cohort", palette = cR, legend = "none",
                add = "loess", conf.int = TRUE) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(legend.position = "none", axis.text.x = element_text(size = 8)) +
  labs(x = "Detected IVA reads (log)", y = "De novo Assembly coverage (%)") +
  ggtitle("C")

# 图 D
p1 <- ggboxplot(data[data$specy == "SARS-2-Cov", ], x = "round", y = "totalNumberContigs",
                color = "specy", add = "jitter", palette = "grey50") +
  theme(legend.position = "none", axis.text.x = element_text(size = 8)) +
  labs(x = "", y = "Number of Contigs") + ggtitle("D")

# 图 E
p2 <- ggboxplot(data, x = "round", y = "coverage", color = "specy",
                add = "jitter", palette = "nejm") +
  facet_wrap(~specy) +
  theme(legend.position = "none", axis.text.x = element_text(size = 8)) +
  labs(x = "", y = "Genomic Coverage (%)") + ggtitle("E")

# 图 F
p3 <- ggboxplot(data, x = "round", y = "signalRatio", color = "specy",
                add = "jitter", palette = "nejm") +
  facet_wrap(~specy) +
  theme(legend.position = "none", axis.text.x = element_text(size = 8)) +
  labs(x = "", y = "Signal Purity (%)") + ggtitle("F")

# 合并所有图形（2行3列）
combined <- ggarrange(p1, p2, p3, p4, p5, p6,
                      ncol = 3, nrow = 2,
                      labels = c("", "", "", "", "", ""),  # 各子图已有标题
                      common.legend = FALSE)

# 保存
ggsave("F_S15.png", plot = combined, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("F_S15.pdf", plot = combined, width = 12, height = 8, dpi = 300, bg = "white")

