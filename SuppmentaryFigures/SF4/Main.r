library(ggplot2)
library(ggsci)
colorsQ = c(pal_nejm("default")(7),"grey50")
base_length <- 30000

# SNPs
rect_snp_hash <- read.table("covid.snp.rec.txt",head=T)
# ypos    xpos
rect_snp_hash$xpos =  as.numeric(rect_snp_hash$xpos)
line_data <- data.frame(
  x = c(0,0),
  xend = c(base_length,base_length),
  y = c(-1,1),
  yend = c(-1,1)
)

annotation_data <- data.frame(
  start = c(266, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 29558),
  end = c(21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 29674),
	gene = c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")
)

dash_line_data <- data.frame(
  x = c(annotation_data$start,annotation_data$end),
  y = -1
)


ggplot() +
  geom_segment(data= rect_snp_hash,aes(x=xpos,xend=xpos,y=-1,yend=1,color=type),linetype="solid") + scale_color_manual(values=colorsQ[1:2]) +
  geom_text(data = annotation_data, aes(x = (start + end) / 2, y = -2, label = gene), color = "black", size = 3, angle = 45, vjust = 0.5, hjust = 0) +
  geom_segment(data = line_data, aes(x = x, xend = xend, y = y, yend = yend), color = "grey", size = 1) +
  geom_segment(data = dash_line_data, aes(x = x, xend = x, y = -1, yend = -1.2),linetype = "solid", color = "grey50") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank()) +
  labs(x = "Reference genome position NC045512.2",y="")
