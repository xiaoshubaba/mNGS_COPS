library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)

Main = read.table("Single.gatk.reads.txt",head=T,sep="\t")
colorsQ = pal_nejm("default")(2)
Main$gatkCov = Main$gatkCov * 100
p <- ggplot(Main, aes(x = gatkReads)) +
  geom_point(aes(y = asmCov, color = "Denovo assemble coverage"), alpha = 0.6,size=1) +
  geom_point(aes(y = Signal, color = "Signal Ratio"), alpha = 0.6) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) + scale_color_manual(values = c("Denovo assemble coverage" = colorsQ[1], "Signal Ratio" = colorsQ[2])) +
  labs(
    y = "Coverage/Signal",
    x = "Number of detected pathogen reads (log10 scale)",
    color = "Metrics"
  ) +
  theme(axis.title = element_text(size = 10), # Adjust axis label sizes
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),legend.position="none") + geom_smooth(data = Main, aes(y = asmCov, color = "Denovo assemble coverage"), method = "gam", se = TRUE, linetype = "solid") +
  geom_smooth(data = Main, aes(y = Signal, color = "Signal Ratio"), method = "gam", se = TRUE, linetype = "solid") +   ylim(0, 100)
