library(ggplot2)
library(ggsci)
base_length <- 1900
gap <- 0.5
rect_data <- read.table("DNA_k75_448.Rect.box.base.txt",head=T)
rect_data$xstart = as.numeric(rect_data$xstart)
rect_data$xend = as.numeric(rect_data$xend)
rect_data$ystart = as.numeric(rect_data$ystart)
rect_data$yend = as.numeric(rect_data$yend)
line_data <- data.frame(
  x = 0,
  xend = base_length,
  y = -1,
  yend = -1
)


L <- read.table("DNA_k75_448.Rect.labels_cov.txt",head=T)
labels <- c("",L$Label)
p = ggplot() +
  geom_rect(data = rect_data, aes(xmin = xstart, xmax = xend, ymin = ystart, ymax = yend), fill = rect_data$fill, color = "white", alpha = 0.7) +
  geom_segment(data = line_data, aes(x = x, xend = xend, y = y, yend = yend), color = "grey", size = 1) +
  scale_y_continuous(breaks = seq(0, 12) * (1 + gap) + gap / 2 -1.2, labels = labels) + 
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank()) + 
  labs(x = "Reference genome position COPs::DNA_k75_448", y = "Patients in the first week")

base_length <- 900
gap <- 0.5
rect_data <- read.table("DNA_k75_10730.Rect.box.base.txt",head=T)
rect_data$xstart = as.numeric(rect_data$xstart)
rect_data$xend = as.numeric(rect_data$xend)
rect_data$ystart = as.numeric(rect_data$ystart)
rect_data$yend = as.numeric(rect_data$yend)
line_data <- data.frame(
  x = 0,
  xend = base_length,
  y = -1,
  yend = -1
)


L <- read.table("DNA_k75_10730.Rect.labels_cov.txt",head=T)
labels <- c("",L$Label)

snps <- read.table("DNA_k75_10730.f1w.snps.txt",head=T)

p1 = ggplot() +
  geom_rect(data = rect_data, aes(xmin = xstart, xmax = xend, ymin = ystart, ymax = yend), fill = rect_data$fill, color = "white", alpha = 0.7) +
  geom_segment(data= snps,aes(x=pos,xend=pos,y=7,yend=9,color=type),linetype="solid") + scale_color_manual(values="#BC3C29FF") +
  geom_segment(data = line_data, aes(x = x, xend = xend, y = y, yend = yend), color = "grey", size = 1) +
  scale_y_continuous(breaks = seq(0, 12) * (1 + gap) + gap / 2 -1.2, labels = labels) + 
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank()) +
  labs(x = "Reference genome position COPs::DNA_k75_10730", y = "Patients in the first week")
