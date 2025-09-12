 library(ggsci)
library(ggpubr)

colorsQ <- c("grey50", pal_nejm("default")(1))
colorsQ_tem <- c("grey50", pal_nejm("default")(1))

# Patient 0 - raw
Main <- read.table("c0.single.table.txt", header = TRUE, sep = "\t")
Main$log_meanDepth <- log(Main$length)
sr <- round(sum(Main[Main$annotation == "Y", ]$length) / sum(Main$length) * 100, 2)
max_gc <- max(Main$gc, na.rm = TRUE)
max_depth <- max(Main$log_meanDepth, na.rm = TRUE)
p0 <- ggscatter(Main, "gc", "log_meanDepth", xlab = "GC content", color = "annotation", alpha = 0.9, legend = "none", palette = colorsQ, size = "length") +   annotate("text", x = max_gc, y = max_depth, hjust = 1, vjust = 1, size = 3, label = paste0("Patient 0\nContigs: ", nrow(Main), "\nSignal: ", sr, "%"))

# Patient 0 - filtered
Main_y <- Main[Main$control.fre <= 1, ]
sr <- round(sum(Main_y[Main_y$annotation == "Y", ]$length) / sum(Main_y$length) * 100, 2)
max_gc <- max(Main_y$gc, na.rm = TRUE)
max_depth <- max(Main_y$log_meanDepth, na.rm = TRUE)
p0_1 <- ggscatter(Main_y, "gc", "log_meanDepth", xlab = "GC content", color = "annotation", alpha = 0.9, legend = "none", palette = colorsQ_tem, size = "length") +  annotate("text", x = max_gc, y = max_depth, hjust = 1, vjust = 1, size = 3, label = paste0("Patient 0 (filtered)\nContigs: ", nrow(Main_y), "\nSignal: ", sr, "%"))


# Patient 1 - raw
Main <- read.table("c1.express.txt", header = TRUE, sep = "\t")
Main$log_meanDepth <- log(Main$length)
sr <- round(sum(Main[Main$annotation == "Y", ]$length) / sum(Main$length) * 100, 2)
max_gc <- max(Main$gc, na.rm = TRUE)
max_depth <- max(Main$log_meanDepth, na.rm = TRUE)
p1 <- ggscatter(Main, "gc", "log_meanDepth", xlab = "GC content", color = "annotation", alpha = 0.9, legend = "none", palette = colorsQ, size = "length") +  annotate("text", x = max_gc, y = max_depth, hjust = 1, vjust = 1, size = 3, label = paste0("Patient 1\nContigs: ", nrow(Main), "\nSignal: ", sr, "%"))

# Patient 1 - filtered
Main_y <- Main[Main$control.fre <= 1, ]
sr <- round(sum(Main_y[Main_y$annotation == "Y", ]$length) / sum(Main_y$length) * 100, 2)
max_gc <- max(Main_y$gc, na.rm = TRUE)
max_depth <- max(Main_y$log_meanDepth, na.rm = TRUE)
p2 <- ggscatter(Main_y, "gc", "log_meanDepth", xlab = "GC content", color = "annotation", alpha = 0.9, legend = "none", palette = colorsQ_tem, size = "length") + annotate("text", x = max_gc, y = max_depth, hjust = 1, vjust = 1, size = 3, label = paste0("Patient 1 (filtered)\nContigs: ", nrow(Main_y), "\nSignal: ", sr, "%"))

# Patient 2 - raw
Main <- read.table("c2.express.txt", header = TRUE, sep = "\t")
Main$log_meanDepth <- log(Main$length)
sr <- round(sum(Main[Main$annotation == "Y", ]$length) / sum(Main$length) * 100, 2)
max_gc <- max(Main$gc, na.rm = TRUE)
max_depth <- max(Main$log_meanDepth, na.rm = TRUE)
p3 <- ggscatter(Main, "gc", "log_meanDepth", xlab = "GC content", color = "annotation", alpha = 0.9, legend = "none", palette = colorsQ, size = "length") +  annotate("text", x = max_gc, y = max_depth, hjust = 1, vjust = 1, size = 3, label = paste0("Patient 2\nContigs: ", nrow(Main), "\nSignal: ", sr, "%"))

# Patient 2 - filtered
Main_y <- Main[Main$control.fre <= 1, ]
sr <- round(sum(Main_y[Main_y$annotation == "Y", ]$length) / sum(Main_y$length) * 100, 2)
max_gc <- max(Main_y$gc, na.rm = TRUE)
max_depth <- max(Main_y$log_meanDepth, na.rm = TRUE)
p4 <- ggscatter(Main_y, "gc", "log_meanDepth", xlab = "GC content", color = "annotation",  alpha = 0.9, legend = "none", palette = colorsQ_tem, size = "length") +  annotate("text", x = max_gc, y = max_depth, hjust = 1, vjust = 1, size = 3, label = paste0("Patient 2 (filtered)\nContigs: ", nrow(Main_y), "\nSignal: ", sr, "%"))

# Combined patients - raw
Main <- read.table("mix.V2.express.txt", header = TRUE, sep = "\t")
Main$log_meanDepth <- log(Main$length)
sr <- round(sum(Main[Main$annotation == "Y", ]$length) / sum(Main$length) * 100, 2)
max_gc <- max(Main$gc, na.rm = TRUE)
max_depth <- max(Main$log_meanDepth, na.rm = TRUE)
p5 <- ggscatter(Main, "gc", "log_meanDepth", xlab = "GC content", color = "annotation",alpha = 0.9, legend = "none", palette = colorsQ, size = "length") +  annotate("text", x = max_gc, y = max_depth, hjust = 1, vjust = 1, size = 3,label = paste0("Combined\nContigs: ", nrow(Main), "\nSignal: ", sr, "%"))

# Combined - case.fre >= 2
Main_m2 <- Main[which(Main$case.fre >=2), ]
sr <- round(sum(Main_m2[Main_m2$annotation == "Y", ]$length) / sum(Main_m2$length) * 100, 2)
max_gc <- max(Main_m2$gc, na.rm = TRUE)
max_depth <- max(Main_m2$log_meanDepth, na.rm = TRUE)
p6 <- ggscatter(Main_m2, "gc", "log_meanDepth", xlab = "GC content", color = "annotation", alpha = 0.9, legend = "none", palette = colorsQ_tem, size = "length") + annotate("text", x = max_gc, y = max_depth, hjust = 1, vjust = 1, size = 3, label = paste0("Case >= 2\nContigs: ", nrow(Main_m2), "\nSignal: ", sr, "%"))

# Combined - case.fre ≥2 & control.fre ≤1
Main_m2_m2 <- Main_m2[which(Main_m2$control.fre <= 1), ]
sr <- round(sum(Main_m2_m2[Main_m2_m2$annotation == "Y", ]$length) / sum(Main_m2_m2$length) * 100, 2)
max_gc <- max(Main_m2_m2$gc, na.rm = TRUE)
max_depth <- max(Main_m2_m2$log_meanDepth, na.rm = TRUE)
p7 <- ggscatter(Main_m2_m2, "gc", "log_meanDepth", xlab = "GC content", color = "annotation", alpha = 0.9, legend = "none", palette = colorsQ_tem, size = "length") +  annotate("text", x = max_gc, y = max_depth, hjust = 1, vjust = 1, size = 3,           label = paste0("Final\nContigs: ", nrow(Main_m2_m2), "\nSignal: ", sr, "%"))

# Arrange all plots
ggarrange(
  ggarrange(p0, p0_1, ncol = 2),		
  ggarrange(p1, p2, ncol = 2),
  ggarrange(p3, p4, ncol = 2),
  ggarrange(p5, p6, p7, ncol = 3),
  nrow = 4
)
