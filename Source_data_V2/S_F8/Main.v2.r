# 加载必要的包
library(ggsci)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)

# ==================== 第一部分：散点图（A, B, C）共享一个图例 ====================
Main <- read.table("YFL.Local_sensitive.V2.table.txt", head = TRUE, sep = "\t")
Main$log_meanDepth <- log(Main$length)

# 定义颜色，确保名称与数据中的 annotation 匹配
colorsQ <- c("grey50", pal_nejm("default")(1))
names(colorsQ) <- c("N", "Y")

# 辅助函数：获取非空颜色
get_colors <- function(data) {
  colorsQ[names(colorsQ) %in% unique(data$annotation)]
}

# 图 A：原始数据（无图例，注释内嵌）
nA <- nrow(Main)
purityA <- round(sum(Main[Main$annotation == "Y", "length"]) / sum(Main$length) * 100, digits = 2)
pA <- ggscatter(Main, "gc", "log_meanDepth",
                xlab = "GC content", ylab = "log(Length)",
                color = "annotation", alpha = 0.9, legend = "none",
                palette = colorsQ, size = "length",
                title = "A. All contigs") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.position = "none") +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("n = ", nA, "\nsignal purity = ", purityA, "%"),
           hjust = 1.1, vjust = 1.1, size = 3, color = "gray30")

# 图 B：case.fre ≥2（无图例，注释内嵌）
Main_m2 <- Main[Main$case.fre >= 2, ]
nB <- nrow(Main_m2)
purityB <- round(sum(Main_m2[Main_m2$annotation == "Y", "length"]) / sum(Main_m2$length) * 100, digits = 2)
colorsB <- get_colors(Main_m2)
pB <- ggscatter(Main_m2, "gc", "log_meanDepth",
                xlab = "GC content", ylab = "log(Length)",
                color = "annotation", alpha = 0.9, legend = "none",
                palette = colorsB, size = "length",
                title = "B. case.fre ≥2") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.position = "none") +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("n = ", nB, "\nsignal purity = ", purityB, "%"),
           hjust = 1.1, vjust = 1.1, size = 3, color = "gray30")

# 图 C：case.fre ≥2 且 control.fre ≤1（无图例，注释内嵌）
Main_m2_m2 <- Main_m2[Main_m2$control.fre <= 1, ]
nC <- nrow(Main_m2_m2)
purityC <- round(sum(Main_m2_m2[Main_m2_m2$annotation == "Y", "length"]) / sum(Main_m2_m2$length) * 100, digits = 2)
colorsC <- get_colors(Main_m2_m2)
pC <- ggscatter(Main_m2_m2, "gc", "log_meanDepth",
                xlab = "GC content", ylab = "log(Length)",
                color = "annotation", alpha = 0.9, legend = "none",
                palette = colorsC, size = "length",
                title = "C. + control.fre ≤1") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.position = "none") +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("n = ", nC, "\nsignal purity = ", purityC, "%"),
           hjust = 1.1, vjust = 1.1, size = 3, color = "gray30")

# 构建共享图例（颜色 + 点大小）
# 颜色图例
temp_plot <- ggplot(Main, aes(x = gc, y = log_meanDepth, color = annotation, size = length)) +
  geom_point(alpha = 0.9) +
  scale_color_manual(values = colorsQ, name = "Annotation", labels = c("N", "Y")) +
  guides(size = "none") +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8))
color_legend <- get_legend(temp_plot)

# 点大小图例（contig length）
size_legend_data <- data.frame(
  length = c(500, 1000, 2000),
  x = 1:3,
  y = 1
)
size_legend_plot <- ggplot(size_legend_data, aes(x = x, y = y, size = length)) +
  geom_point(shape = 21, fill = "grey50", color = "black", alpha = 0.9) +
  scale_size_continuous(name = "Contig length (bp)", range = c(2, 8), breaks = c(500, 1000, 2000)) +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8))
size_legend <- get_legend(size_legend_plot)

# 上下组合图例
combined_legend <- plot_grid(color_legend, size_legend, ncol = 1, rel_heights = c(1, 0.8))

# ==================== 第二部分：对齐图（pd, pe, pf）保留各自图例 ====================
# ---- pd: 基于 YFL.k141_396.false.alignment.txt ----
data <- read.table("YFL.k141_396.false.alignment.txt", 
                   header = TRUE, fill = TRUE, stringsAsFactors = FALSE)

data_sorted <- data %>% arrange(fragment, sample)

data_processed <- data_sorted %>%
  mutate(alignment = ifelse(is.na(alignment), "", as.character(alignment)),
         y = rev(row_number())) %>%
  group_by(sample) %>%
  mutate(group_size = n(),
         group_y_min = min(y) - 0.4,
         group_y_max = max(y) + 0.4) %>%
  ungroup()

parse_cigar <- function(cigar, start) {
  if (cigar == "" || is.na(cigar)) return(data.frame())
  ops <- str_match_all(cigar, "(\\d+)([A-Z])")[[1]]
  if (nrow(ops) == 0) return(data.frame())
  df <- data.frame(length = as.numeric(ops[, 2]), type = ops[, 3], stringsAsFactors = FALSE)
  df$ref_start <- cumsum(c(start, df$length))[1:nrow(df)]
  df$ref_end <- df$ref_start + df$length - 1
  df$ref_start[1] <- start
  df$ref_end[1] <- start + df$length[1] - 1
  return(df)
}

alignments <- lapply(1:nrow(data_processed), function(i) {
  df <- parse_cigar(data_processed$alignment[i], data_processed$start[i])
  if (nrow(df) > 0) {
    df$y <- data_processed$y[i]
    df$sample <- data_processed$sample[i]
  }
  return(df)
}) %>% bind_rows()

operation_colors <- c("S" = "#FF9999", "M" = "#99CCFF", "D" = "#99FF99")

x_min <- min(data_processed$start, na.rm = TRUE) - 50
x_max <- max(alignments$ref_end, na.rm = TRUE) + 50

group_lines <- data_processed %>%
  distinct(sample, group_y_min, group_y_max, group_size) %>%
  filter(group_size >= 2) %>%
  mutate(group_x = x_min - 50, group_mid = (group_y_min + group_y_max)/2)

base_line <- data.frame(x = c(x_min, x_max), y = c(0, 0))
position_ticks <- data.frame(x = seq(floor(x_min/50)*50, ceiling(x_max/50)*50, by = 50))

alignments$type <- factor(alignments$type, levels = c("S", "M", "D"))

pd <- ggplot() +
  geom_line(data = base_line, aes(x = x, y = y), color = "gray70", size = 1) +
  geom_segment(data = group_lines, aes(x = group_x, xend = group_x, y = group_y_min, yend = group_y_max), color = "gray70", size = 0.5) +
  geom_text(data = group_lines, aes(x = group_x - 2, y = group_mid, label = paste0("n=", group_size)), size = 3, color = "gray30", hjust = 1) +
  geom_segment(data = alignments, aes(x = ref_start, xend = ref_end, y = y, yend = y, color = type), size = 0.8, lineend = "butt") +
  geom_segment(data = position_ticks, aes(x = x, xend = x, y = 0, yend = -0.15), color = "gray50") +
  geom_text(data = position_ticks, aes(x = x, y = -0.4, label = x), size = 2.5, color = "gray40") +
  scale_color_manual(name = "Operation", values = operation_colors,
                     labels = c("S: Soft Clip", "M: Match", "D: Deletion"), drop = FALSE) +
  labs(x = "Reference Position", y = "", title = "D. False alignment details (k141_396)") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "right", plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.title = element_blank(), legend.text = element_text(size = 8),
        legend.key.size = unit(0.6, "lines")) +
  coord_cartesian(xlim = c(x_min - 60, x_max), ylim = c(-0.7, max(data_processed$y) + 0.5))

# ---- pe: 基于 YFL.alignments.CIGAR.txt ----
alignment_data <- read.table("YFL.alignments.CIGAR.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

fragment_counts <- alignment_data %>% group_by(fragment) %>% summarise(read_count = n())
top_fragments <- fragment_counts %>% filter(read_count > 10)

subset_data <- alignment_data %>%
  inner_join(top_fragments, by = "fragment") %>%
  mutate(fragment_label = paste0(fragment, " (n=", read_count, ")"))

pe <- ggplot(subset_data, aes(x = fragment_label, y = start)) +
  geom_jitter(width = 0.3, height = 0, size = 1.2, alpha = 0.7, color = "steelblue") +
  theme_minimal() +
  labs(x = "Fragment (with read count)", y = "Start Position", title = "E. Fragments with >10 reads") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 9),
        plot.title = element_text(size = 10, face = "bold"))

# ---- pf: 特定片段 k141_80 ----
data_k141 <- alignment_data[alignment_data$fragment == "k141_80", ]
if (nrow(data_k141) > 0) {
  data_sorted_k <- data_k141 %>% arrange(fragment, sample)
  data_processed_k <- data_sorted_k %>%
    mutate(alignment = ifelse(is.na(alignment), "", as.character(alignment)),
           y = rev(row_number())) %>%
    group_by(sample) %>%
    mutate(group_size = n(),
           group_y_min = min(y) - 0.4,
           group_y_max = max(y) + 0.4) %>%
    ungroup()
  
  alignments_k <- lapply(1:nrow(data_processed_k), function(i) {
    df <- parse_cigar(data_processed_k$alignment[i], data_processed_k$start[i])
    if (nrow(df) > 0) {
      df$y <- data_processed_k$y[i]
      df$sample <- data_processed_k$sample[i]
    }
    return(df)
  }) %>% bind_rows()
  
  x_min_k <- min(data_processed_k$start, na.rm = TRUE) - 50
  x_max_k <- max(alignments_k$ref_end, na.rm = TRUE) + 50
  group_lines_k <- data_processed_k %>%
    distinct(sample, group_y_min, group_y_max, group_size) %>%
    filter(group_size >= 2) %>%
    mutate(group_x = x_min_k - 50, group_mid = (group_y_min + group_y_max)/2)
  base_line_k <- data.frame(x = c(x_min_k, x_max_k), y = c(0, 0))
  position_ticks_k <- data.frame(x = seq(floor(x_min_k/50)*50, ceiling(x_max_k/50)*50, by = 50))
  alignments_k$type <- factor(alignments_k$type, levels = c("S", "M", "D"))
  
  pf <- ggplot() +
    geom_line(data = base_line_k, aes(x = x, y = y), color = "gray70", size = 1) +
    geom_segment(data = group_lines_k, aes(x = group_x, xend = group_x, y = group_y_min, yend = group_y_max), color = "gray70", size = 0.5) +
    geom_text(data = group_lines_k, aes(x = group_x - 2, y = group_mid, label = paste0("n=", group_size)), size = 3, color = "gray30", hjust = 1) +
    geom_segment(data = alignments_k, aes(x = ref_start, xend = ref_end, y = y, yend = y, color = type), size = 0.8, lineend = "butt") +
    geom_segment(data = position_ticks_k, aes(x = x, xend = x, y = 0, yend = -0.15), color = "gray50") +
    geom_text(data = position_ticks_k, aes(x = x, y = -0.4, label = x), size = 2.5, color = "gray40") +
    scale_color_manual(name = "Operation", values = operation_colors,
                       labels = c("S: Soft Clip", "M: Match", "D: Deletion"), drop = FALSE) +
    labs(x = "Reference Position", y = "", title = "F. Detailed alignment for fragment k141_80") +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "right", plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
          legend.title = element_blank(), legend.text = element_text(size = 8),
          legend.key.size = unit(0.6, "lines")) +
    coord_cartesian(xlim = c(x_min_k - 60, x_max_k), ylim = c(-0.7, max(data_processed_k$y) + 0.5))
} else {
  pf <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data for fragment k141_80", size = 4) + theme_void()
}

# ==================== 组合图形 ====================
# 2行3列网格，不添加字母标签
grid_plots <- plot_grid(pA, pB, pC, pd, pe, pf, 
                        ncol = 3, nrow = 2,
                        labels = NULL)   # 删除字母标签

# 右侧放置共享图例，宽度比例减小至 0.2 以节约空间
final_plot <- plot_grid(grid_plots, combined_legend, 
                        ncol = 2, 
                        rel_widths = c(1, 0.2))

# 保存
ggsave("F_S8_combined.png", plot = final_plot, width = 16, height = 10, dpi = 300, bg = "white")
ggsave("F_S8_combined.pdf", plot = final_plot, width = 16, height = 10, dpi = 300, bg = "white")

