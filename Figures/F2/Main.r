library(ggplot2)
library(dplyr)
library(ggsci)
colorsQ = c(pal_nejm("default")(2),"grey50")
# Read the data
df <- read.table("single.sample.asm.txt", header = TRUE)

# Multiply gatk_cov by 100
df$gatk_cov <- df$gatk_cov * 100
df_high <- df %>% filter(gatk_Reads > 2034)
df_low <- df %>% filter(gatk_Reads <= 2034)

# Calculate Pearson correlation coefficients and p-values
pearson_low_cov <- cor.test(df_low$gatk_Reads, df_low$cov, method = "pearson")
pearson_high_cov <- cor.test(df_high$gatk_Reads, df_high$cov, method = "pearson")

pearson_low_gatk_cov <- cor.test(df_low$gatk_Reads, df_low$gatk_cov, method = "pearson")
pearson_high_gatk_cov <- cor.test(df_high$gatk_Reads, df_high$gatk_cov, method = "pearson")

pearson_low_signal_ratio <- cor.test(df_low$gatk_Reads, df_low$signal_ratio, method = "pearson")
pearson_high_signal_ratio <- cor.test(df_high$gatk_Reads, df_high$signal_ratio, method = "pearson")

p <- ggplot(df, aes(x = gatk_Reads)) +
  geom_point(aes(y = cov, color = "Denovo assemble genomic coverage"), alpha = 0.6) +
  geom_point(aes(y = gatk_cov, color = "Mapped genomic coverage"), alpha = 0.6) +
  geom_point(aes(y = signal_ratio, color = "Signal Ratio"), alpha = 0.6) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "Detected Pathogen Reads (log10 scale)", y = "Value") +
  scale_color_manual(values = c("Denovo assemble genomic coverage" = colorsQ[1], "Mapped genomic coverage" = colorsQ[2], "Signal Ratio" = colorsQ[3])) +
  geom_smooth(data = df_low, aes(y = cov, color = "Denovo assemble genomic coverage"), method = "gam", se = TRUE, linetype = "solid") +
  geom_smooth(data = df_high, aes(y = cov, color = "Denovo assemble genomic coverage"), method = "gam", se = TRUE, linetype = "solid") +
  geom_smooth(data = df_low, aes(y = gatk_cov, color = "Mapped genomic coverage"), method = "gam", se = TRUE, linetype = "solid") +
  geom_smooth(data = df_high, aes(y = gatk_cov, color = "Mapped genomic coverage"), method = "gam", se = TRUE, linetype = "solid") +
  geom_smooth(data = df_low, aes(y = signal_ratio, color = "Signal Ratio"), method = "gam", se = TRUE, linetype = "solid") +
  geom_smooth(data = df_high, aes(y = signal_ratio, color = "Signal Ratio"), method = "gam", se = TRUE, linetype = "solid") +
  geom_vline(xintercept = 2034, linetype = "dashed", color = "grey") +
  annotate("text", x = 2034, y = -1, label = "x = 2034", hjust = -0.1, color = "grey50") +
  annotate("text", x = 0, y = 100, label = paste0("Low (Denovo coverage): ", round(pearson_low_cov$estimate, 2), ", p = ", round(pearson_low_cov$p.value, 2)), hjust = -0.1, color = colorsQ[1], size = 4) +
  annotate("text", x = 10^4, y = 0, label = paste0("High (Denovo coverage): ", round(pearson_high_cov$estimate, 2), ", p = ", round(pearson_high_cov$p.value, 2)), hjust = -0.1, color = colorsQ[1], size = 4) +
  annotate("text", x = 0, y = 90, label = paste0("Low (Mapped coverage): ", round(pearson_low_gatk_cov$estimate, 2), ", p = ", round(pearson_low_gatk_cov$p.value, 2)), hjust = -0.1, color = colorsQ[2], size = 4) +
  annotate("text", x = 10^4, y = 10, label = paste0("High (Mapped coverage): ", round(pearson_high_gatk_cov$estimate, 2), ", p = ", round(pearson_high_gatk_cov$p.value, 2)), hjust = -0.1, color = colorsQ[2], size = 4) +
  annotate("text", x = 0, y = 80, label = paste0("Low (Signal Ratio): ", round(pearson_low_signal_ratio$estimate, 2), ", p = ", round(pearson_low_signal_ratio$p.value, 2)), hjust = -0.1, color = colorsQ[3], size = 4) +
  annotate("text", x = 10^4, y = 20, label = paste0("High (Signal Ratio): ", round(pearson_high_signal_ratio$estimate, 2), ", p = ", round(pearson_high_signal_ratio$p.value, 2)), hjust = -0.1, color = colorsQ[3], size = 4) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.box = "horizontal")
