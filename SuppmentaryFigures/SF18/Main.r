library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)
setwd("/Users/liang/Desktop/project/mNGS_COVID/Figures/S_F29")
Main = read.table("S.txt",head=T)
Main$prj = factor(Main$prj,levels=c("p1","p2","HIV","LASV","Ebola"))
Main$round = factor(Main$round,levels=c(1,2,3,4))

Main_long <- Main %>%
  pivot_longer(cols = c(cov, signal), names_to = "Metric", values_to = "Value")

# Create the prj_Metric column after pivoting
Main_long <- Main_long %>%
  mutate(prj_Metric = interaction(Metric, prj, sep = " - "))


# Create the plot
p1 = ggplot(Main_long, aes(x = round, y = Value, color = Metric)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.7) +
  facet_wrap(~ prj_Metric, scales = "free_y", strip.position = "top", ncol = length(unique(Main_long$Metric))) +
  scale_color_manual(values = c(nCtg = "grey20", "cov" = "#BC3C29FF", "signal" = "grey50")) +
  labs(x = "Round", y = "Value", color = "Metric") +
  theme(axis.title = element_text(size = 10),
        legend.position = "none",
        strip.placement = "outside",
        strip.text = element_text(size = 8),
        strip.text.x = element_text(angle = 0),
        strip.text.y = element_text(angle = 0))
