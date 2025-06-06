library(ggpubr)
library(dplyr)
library(tidyr)

# ---- Load Data ----
data <- read.table("Multiple.source.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ---- Reshape to Long Format ----
long_data <- data %>%
  select(project, Signal, SignalCombined) %>%
  pivot_longer(cols = c(Signal, SignalCombined),
               names_to = "type", values_to = "value")

# ---- Boxplot ----
ggboxplot(long_data,
          x = "type",
          y = "value",
          color = "black",   # outline only
          fill = NA,         # no fill
          outlier.shape = NA,
          facet.by = "project",
          add = "none") +
  stat_compare_means(
    method = "wilcox.test",
    paired = TRUE,
    label = "p.format",  # show actual p-value
    label.y.npc = "top"
  ) +
  labs(title = "Signal vs SignalCombined per Project",
       x = "",
       y = "Signal Value") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
