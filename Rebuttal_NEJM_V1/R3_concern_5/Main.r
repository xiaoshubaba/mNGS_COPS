library(tidyverse)
library(ggridges)
library(patchwork)

# File list
files <- c("prj.01.control.fre.txt", "prj.02.control.txt", "prj.04.control.fre.txt",
           "prj.06.ebola.fre.txt", "prj.07.lasv.control.fre.txt", "prj.08.den.fre.txt",
           "prj.17.iva.control.fre.txt", "prj.18.ivb.control.fre.txt", "prj.19.rha.control.fre.txt",
           "prj.20.hmb.control.fre.txt", "prj.21.hrsv.control.fre.txt","prj.05.zika.control.fre.txt","prj.09.chiv.fre.txt")

# Virus name mapping
virus_map <- tibble(
  prj = c("prj.01", "prj.02", "prj.04", "prj.06", "prj.07", "prj.08",
          "prj.17", "prj.18", "prj.19", "prj.20", "prj.21","prj.05","prj.09"),
  virus = c("SARS-CoV-2", "SARS-CoV-2", "HIV", "Ebola", "LASV", "Dengue virus",
            "Influenza A", "Influenza B", "Rhinovirus", "Adenovirus", "Parainfluenza virus","Zika","chikv")
)

# 1. Compute key stats --------------------------------------------------------
results <- data.frame(file = character(),
                      pct_90_controls = numeric(),
                      zero_align_controls = numeric())

for (file in files) {
  data <- read.delim(file, header = TRUE)
  pct_90 <- data %>%
    filter(ratio >= 0.9) %>%
    slice(1) %>%
    pull(order)
  zero_align <- mean(data$fre == 0) * 100
  results <- rbind(results, data.frame(file = file,
                                       pct_90_controls = pct_90,
                                       zero_align_controls = zero_align))
}

# Summary stats
mean_pct_90 <- mean(results$pct_90_controls, na.rm = TRUE)
sd_pct_90 <- sd(results$pct_90_controls, na.rm = TRUE)
mean_zero_align <- mean(results$zero_align_controls, na.rm = TRUE)
sd_zero_align <- sd(results$zero_align_controls, na.rm = TRUE)

# 2. Rarefaction curves -------------------------------------------------------
all_data <- list()
for (file in files) {
  prj_id <- str_extract(file, "prj\\.[0-9]+")
  df <- read.delim(file, header = TRUE) %>%
    mutate(percentage = order / max(order) * 100,
           prj = prj_id)
  df <- left_join(df, virus_map, by = "prj")
  all_data[[file]] <- df
}
combined <- bind_rows(all_data)

# Colors by virus
virus_colors <- scales::hue_pal()(length(unique(combined$virus)))
names(virus_colors) <- unique(combined$virus)

# Plot A: Rarefaction curves
pA <- ggplot(combined, aes(x = percentage, y = ratio, color = virus)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
  scale_color_manual(values = virus_colors) +
  labs(title = "A", x = "Control Subset (%)", y = "Cumulative Filtering Ratio", color = "Virus") +
  theme_minimal() +
  theme(legend.position = c(0.8, 0.4),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing.y = unit(0.2, "cm"),
        legend.key.size = unit(0.3, "cm"),
        legend.direction = "vertical",
        plot.title = element_text(face = "bold", size = 14, hjust = -0.1)) +
  guides(color = guide_legend(ncol = 3, byrow = TRUE))

# Plot B
pB <- ggplot(results, aes(x = pct_90_controls)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  geom_vline(xintercept = mean_pct_90, color = "red", linetype = "dashed") +
  annotate("text", x = mean_pct_90 + 5, y = 0.05,
           label = sprintf("Mean: %.1f", mean_pct_90), color = "red") +
  labs(title = "B", x = "Control Order Index", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = -0.1))

# Plot C
pC <- ggplot(results, aes(x = zero_align_controls)) +
  geom_density(fill = "forestgreen", alpha = 0.6) +
  geom_vline(xintercept = mean_zero_align, color = "red", linetype = "dashed") +
  annotate("text", x = mean_zero_align + 5, y = 0.05,
           label = sprintf("Mean: %.1f%%", mean_zero_align), color = "red") +
  labs(title = "C", x = "Zero-Aligned Control (%)", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = -0.1))

# Combine and save
final_plot <- (pA) / (pB + pC) +
  plot_annotation(title = "Control Filtering Efficiency Analysis",
                  subtitle = "Performance assessment under gradient-based filtering design")

print(final_plot)
ggsave("rarefaction_analysis_labeled_by_virus.png", final_plot, width = 12, height = 10, dpi = 300)
