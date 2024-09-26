library(ggpubr)
library(ggplot2)
Main = read.table("./bootstra.txt",head=T)
Main$MixNum = factor(Main$MixNum)
mean_vals <- aggregate(Common ~ Prj + MixNum, data = Main, mean)
sd_vals <- aggregate(Common ~ Prj + MixNum, data = Main, sd)
mean_sd_vals <- merge(mean_vals, sd_vals, by = c("Prj", "MixNum"), suffixes = c("_mean", "_sd"))

# Set y_position based on the MixNum condition
mean_sd_vals$y_position <- ifelse(mean_sd_vals$MixNum == 10, 1000, 0)

# Create the histogram with annotations
p <- gghistogram(Main, x = "Common", add = "mean",
                 facet.by = "Prj", nrow=5,color = "MixNum", scales = "free",
                 fill = "MixNum", alpha = 0.2,legend="none") +
     labs(fill = "Number of randomly\nChosen samples", color = "Number of randomly\nChosen samples") +
     geom_text(data = mean_sd_vals, aes(x = Common_mean, y = y_position,
                                        label = paste0("Mean: ", round(Common_mean, 2))),
               vjust = -1.5, hjust = 0.5, color = "black", size = 3) +
     geom_text(data = mean_sd_vals, aes(x = Common_mean, y = y_position,
                                        label = paste0("SD: ", round(Common_sd, 2))),
               vjust = -2.5, hjust = 0.5, color = "red", size = 3)(base)
