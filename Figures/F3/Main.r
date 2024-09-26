library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
Main = read.table("Main.5.txt",head=T,sep="\t")
Main = Main[,c(1,2,3,4,5,6,7)]
Main$prj = factor(Main$prj,levels=c("SARS-Cov-2","SARS-Cov-2 Study 1","SARS-Cov-2 Study 2","HIV","LASV","Zika","Ebola","CHIKV","DENG"))
Main$mixNum = factor(Main$mixNum,levels=c(5,10))
Main$round = factor(Main$round,levels=c(1,2,3))

Main_long <- Main %>%
  pivot_longer(cols = c(nCtg,cov, signal), names_to = "Metric", values_to = "Value")
# Create the prj_Metric column after pivoting
#Main_long <- Main_long %>%
#  mutate(prj_Metric = interaction(Metric, prj, sep = " - "))

# Filter the data
Main_long_5 = Main_long[Main_long$mixNum == 5, ]
# order
Main_long_5$Metric = factor(Main_long_5$Metric ,levels=c("nCtg","cov","signal"))

# Create the plot
p1 = ggplot(Main_long_5, aes(x = round, y = Value, color = Metric)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.7) +
  facet_grid(rows = vars(prj), cols = vars(Metric), scales = "free_y", switch = "y",labeller = labeller(Metric = c("cov" = "Denovo assembled genomic coverage","signal" = "Signal ratio","nCtg" = "Number of genomic fragments")))+
  scale_color_manual(values = c(nCtg = "grey20", "cov" = "#BC3C29FF", "signal" = "grey50")) +
  labs(x = "Filter Steps", y = "Value", color = "Metric") +
  theme(axis.title = element_text(size = 10),
        legend.position = "none",
        strip.placement = "outside",
        strip.text = element_text(size = 8),
        strip.text.x = element_text(angle = 0),
        strip.text.y = element_text(angle = 0)) + guides(color = FALSE)  # 1) Hide the legend


#p2
    Main = read.table("Main.5.txt",head=T,sep="\t")
    Main = Main[,c(1,2,3,4,6,7)]
    Main$prj = factor(Main$prj,levels=c("SARS-Cov-2","SARS-Cov-2 Study 1","SARS-Cov-2 Study 2","HIV","LASV","Zika","Ebola","CHIKV","DENG"))
    Main$mixNum = factor(Main$mixNum,levels=c(5,10))
    Main$round = factor(Main$round,levels=c(1,2,3))

    Main_long <- Main %>%
      pivot_longer(cols = c(cov, signal), names_to = "Metric", values_to = "Value")
    # Create the prj_Metric column after pivoting
    #Main_long <- Main_long %>%
    #  mutate(prj_Metric = interaction(Metric, prj, sep = " - "))

    # Filter the data
    Main_long_5 = Main_long[Main_long$mixNum == 5, ]
    # order
    Main_long_5$Metric = factor(Main_long_5$Metric ,levels=c("cov","signal"))

    # Create the plot
    p2 = ggplot(Main_long_5, aes(x = round, y = Value, color = Metric)) +
      geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
      geom_jitter(position = position_dodge(width = 0.75), alpha = 0.7) +
      facet_grid(rows = vars(prj), cols = vars(Metric), scales = "free_y", switch = "y",labeller = labeller(Metric = c("cov" = "Denovo assembled genomic coverage","signal" = "Signal ratio","nCtg" = "Number of genomic fragments")))+
      scale_color_manual(values = c(nCtg = "grey20", "cov" = "#BC3C29FF", "signal" = "grey50")) +
      labs(x = "Filter Steps", y = "Value", color = "Metric") +
      theme(axis.title = element_text(size = 10),
            legend.position = "none",
            strip.placement = "outside",
            strip.text = element_text(size = 8),
            strip.text.x = element_text(angle = 0),
            strip.text.y = element_text(angle = 0)) + guides(color = FALSE)  # 1) Hide the legend
