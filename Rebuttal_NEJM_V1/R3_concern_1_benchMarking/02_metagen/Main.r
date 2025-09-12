library(tidyverse)
library(ggpubr)
library(patchwork)

# Load data
data <- read.delim("TABS/whole.cluser.txt", header = TRUE, sep = "\t")

# Plot A: Distribution of clustering metrics
p1 <- data %>%
  pivot_longer(c(homogeneityPathogenScore, PathogenPurity_count, PathogenPurity_length),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         "homogeneityPathogenScore" = "Homogeneity",
                         "PathogenPurity_count" = "Count Purity",
                         "PathogenPurity_length" = "Length Purity")) %>%
  ggplot(aes(value, fill = metric)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~metric, scales = "free") +
  labs(title = "Distribution of Clustering Characteristics",
       x = "Metric Value", y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")

# Plot B: Homogeneity vs. Count Purity by project
p2 <- data %>%
  ggplot(aes(homogeneityPathogenScore, PathogenPurity_count, color = prj)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Homogeneity vs. Count Purity by Project",
       x = "Homogeneity", y = "Count Purity") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 3))

# Plot C: Count vs. Length Purity colored by homogeneity
p3 <- data %>%
  ggplot(aes(PathogenPurity_count, PathogenPurity_length)) +
  geom_point(aes(color = homogeneityPathogenScore), alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_color_viridis_c(option = "plasma") +
  labs(title = "Count vs. Length Purity",
       x = "Count Purity",
       y = "Length Purity",
       color = "Homogeneity") +
  theme_minimal()

# Plot D1: Number of clusters per project
p4 <- ggboxplot(data, x = "prj", y = "n_cluster", 
                fill = "prj", palette = "Set3",
                xlab = FALSE, ylab = "Number of Clusters",
                title = "Cluster Count Distribution by Project") +
  rotate_x_text(45) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8))

# Plot D2: Clustering metrics by project
metrics_data <- data %>%
  dplyr::select(prj, 
         Homogeneity = homogeneityPathogenScore,
         `Count Purity` = PathogenPurity_count,
         `Length Purity` = PathogenPurity_length) %>%
  pivot_longer(-prj, names_to = "metric", values_to = "value")

p5 <- ggplot(metrics_data, aes(x = prj, y = value, fill = metric)) +
  geom_boxplot() +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  labs(title = "Clustering Metrics by Project",
       x = "Project", y = "Value") +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Combine all plots and label with A–D
final_plot <- (p1 | p3) / (p2) / (p4 | p5) +
  plot_annotation(title = "Clustering Characteristics Across 833 Outbreak Simulations",
                  tag_levels = 'A') &
  theme(plot.title = element_text(face = "bold", size = 12))

# Show and save
print(final_plot)
ggsave("clustering_characteristics.png", final_plot, width = 14, height = 12, dpi = 300)
