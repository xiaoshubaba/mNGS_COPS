library(ggplot2)
library(ggsci)
colorsQ = c(pal_nejm("default")(7),"grey50")
base_length <- 30000

gap <- 0.5

rect_data <- read.table("RNA_k75_9508.Rect.box.base.txt",head=T)
rect_data$xstart = as.numeric(rect_data$xstart)
rect_data$xend = as.numeric(rect_data$xend)
rect_data$ystart = as.numeric(rect_data$ystart)
rect_data$yend = as.numeric(rect_data$yend)
# SNPs
rect_snp_hash <- read.table("RNA_k75_9508.Rect.snp.txt",head=T)
# ypos    xpos
rect_snp_hash$ypos =  as.numeric(rect_snp_hash$ypos)
rect_snp_hash$xpos =  as.numeric(rect_snp_hash$xpos)
line_data <- data.frame(
  x = 0,
  xend = base_length,
  y = -1,
  yend = -1
)

annotation_data <- data.frame(
  start = c(302, 21590, 25405, 26257, 26535, 27214, 27406, 27768, 27906, 28286, 29561),
  end = c(21582, 25396, 26232, 26484, 27203, 27399, 27771, 27899, 28271, 29536, 29677),
  gene = c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
)

dash_line_data <- data.frame(
  x = c(annotation_data$start,annotation_data$end),
  y = -1
)

L <- read.table("RNA_k75_9508.Rect.labels_cov.txt",head=T)
labels <- c("",L$Label)

ggplot() +
  geom_rect(data = rect_data, aes(xmin = xstart, xmax = xend, ymin = ystart, ymax = yend), fill = rect_data$fill, color = "white", alpha = 0.7) +
  geom_segment(data= rect_snp_hash,aes(x=xpos,xend=xpos,y=ypos-1,yend=ypos+1,color=type),linetype="solid") + scale_color_manual(values=colorsQ[1:2]) +
  geom_text(data = annotation_data, aes(x = (start + end) / 2, y = -2, label = gene), color = "black", size = 3, angle = 45, vjust = 0.5, hjust = 0) +
  geom_segment(data = line_data, aes(x = x, xend = xend, y = y, yend = yend), color = "grey", size = 1) +
  geom_segment(data = dash_line_data, aes(x = x, xend = x, y = -1.2, yend = -0.8),linetype = "dashed", color = "grey50") +
  scale_y_continuous(breaks = seq(0, 12) * (1 + gap) + gap / 2 -1.2, labels = labels) + 
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank()) + 
  labs(x = "Reference genome position COPs::RNA_k75_9508", y = "Patients in the first week")
