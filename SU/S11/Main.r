library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggsci)
library(scales)

data <- read.table("Main.6.V2.txt", header = TRUE, sep = "")
data$aR <- data$annotatedRatio / (100 - data$signal)
filtered_data <- data %>% filter(round == 3, signal < 100)

plot_data <- filtered_data %>%
  mutate(
    signal_prop = signal / 100,
    annotated_prop = annotatedRatio / 100,
    unannotated_prop = (100 - signal - annotatedRatio) / 100
  ) %>%
  select(prj, signal_prop, annotated_prop, unannotated_prop)

project_means <- plot_data %>%
  group_by(prj) %>%
  summarise(across(c(signal_prop, annotated_prop, unannotated_prop), mean, na.rm = TRUE)) %>%
  mutate(signal_rank = signal_prop)

long_data_A <- project_means %>%
  pivot_longer(
    cols = c(signal_prop, annotated_prop, unannotated_prop),
    names_to = "category",
    values_to = "proportion"
  ) %>%
  mutate(
    signal_rank = rep(project_means$signal_rank, each = 3),
    category = factor(
      category,
      levels = c("signal_prop", "annotated_prop", "unannotated_prop"),
      labels = c("Signal", "Annotated", "Unannotated")
    )
  )

long_data_B <- project_means %>%
  pivot_longer(
    cols = c(annotated_prop, unannotated_prop),
    names_to = "category",
    values_to = "proportion"
  ) %>%
  mutate(
    signal_rank = rep(project_means$signal_rank, each = 2),
    category = factor(
      category,
      levels = c("annotated_prop", "unannotated_prop"),
      labels = c("Annotated", "Unannotated")
    )
  )

nejm_colors <- pal_nejm("default")(2)
custom_palette_A <- c("Signal" = nejm_colors[1], "Annotated" = nejm_colors[2], "Unannotated" = "grey50")
custom_palette_B <- c("Annotated" = nejm_colors[2], "Unannotated" = "grey50")
alpha_value <- 0.6

p_A <- ggplot(long_data_A, aes(x = proportion, y = reorder(prj, signal_rank), fill = category)) +
  geom_col(width = 0.85, position = position_stack(reverse = TRUE), alpha = alpha_value) +
  scale_fill_manual(values = custom_palette_A, guide = guide_legend(reverse = TRUE)) +
  geom_text(
    aes(label = ifelse(proportion > 0.1, percent(proportion, accuracy = 1), "")),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    color = "white", size = 3.5
  ) +
  labs(x = "Proportion", y = "Project") +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_text(margin = margin(r = 15)),
    plot.margin = margin(10, 10, 10, 10)
  )
