# =============================================================================
# Paired Wilcoxon test boxplots for simulated outbreak data with mock controls
# Data: S6.Simulated.outbreaks.pops.mock.txt
# =============================================================================

# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggsci)

# NEJM color palette (first two colors)
nejm_cols <- pal_nejm("default")(4)
col_base <- nejm_cols[1]   # color for base group (without mock)
col_mock <- nejm_cols[2]   # color for mock group (with mock)

# Read data
data <- read.table("S6.Simulated.outbreaks.pops.mock.txt",
                   header = TRUE,
                   sep = "\t",
                   stringsAsFactors = FALSE)

# Check required columns
if (!all(c("sample", "prj", "signal") %in% colnames(data))) {
  stop("Missing required columns: sample, prj, signal")
}

# Remove rows where signal == 0
data <- data[data$signal != 0, ]
cat("Rows after removing signal=0:", nrow(data), "\n")

# Remove duplicate sample-prj pairs
data <- data %>% distinct(sample, prj, .keep_all = TRUE)

# Identify mock projects and their base counterparts
prj_levels <- unique(data$prj)
cat("Unique prj levels:\n")
print(prj_levels)

mock_prjs <- grep("_mock$", prj_levels, value = TRUE)
base_prjs <- gsub("_mock$", "", mock_prjs)

valid_pairs <- data.frame(base = character(), mock = character(), stringsAsFactors = FALSE)
for (i in seq_along(mock_prjs)) {
  base <- base_prjs[i]
  mock <- mock_prjs[i]
  if (base %in% prj_levels) {
    valid_pairs <- rbind(valid_pairs, data.frame(base = base, mock = mock))
  } else {
    warning("Base project '", base, "' not found for mock '", mock, "'")
  }
}

if (nrow(valid_pairs) == 0) {
  stop("No valid base-mock pairs found.")
}

# Function to create paired plot for one base-mock pair
plot_pair <- function(df, base_group, mock_group) {
  # Filter data for the two groups
  sub_df <- df %>% filter(prj %in% c(base_group, mock_group))

  # Convert to wide format: one row per sample
  wide_df <- sub_df %>%
    select(sample, prj, signal) %>%
    pivot_wider(id_cols = sample, names_from = prj, values_from = signal)

  # Keep only samples present in both groups
  wide_df_complete <- wide_df %>% drop_na()
  n_samples <- nrow(wide_df_complete)
  if (n_samples == 0) {
    warning("No common samples between ", base_group, " and ", mock_group)
    return(NULL)
  }

  # Compute improvement percentage
  wide_df_complete$improvement_pct <- ((wide_df_complete[[mock_group]] - wide_df_complete[[base_group]]) /
                                         wide_df_complete[[base_group]]) * 100
  avg_improvement <- mean(wide_df_complete$improvement_pct, na.rm = TRUE)

  # Paired Wilcoxon test
  wilcox_test <- wilcox.test(wide_df_complete[[mock_group]],
                              wide_df_complete[[base_group]],
                              paired = TRUE)

  # Reshape to long format for plotting
  plot_data <- wide_df_complete %>%
    pivot_longer(cols = c(base_group, mock_group),
                 names_to = "Group",
                 values_to = "Value") %>%
    mutate(Group = factor(Group, levels = c(base_group, mock_group)))

  # Title without "samples"
  title <- paste0(base_group, " (n = ", n_samples, ")")
  subtitle <- paste0("Average improvement: ", round(avg_improvement, 1),
                     "% (p = ", format.pval(wilcox_test$p.value, digits = 2), ")")

  # Generate the plot
  p <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA) +
    geom_line(aes(group = sample), color = "gray50", alpha = 0.5, linewidth = 0.4) +
    geom_point(aes(group = sample), size = 1.8, alpha = 0.7) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = "Signal"
    ) +
    scale_fill_manual(values = setNames(c(col_base, col_mock), c(base_group, mock_group))) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 9),
      panel.grid.minor = element_blank()
    )

  return(p)
}

# Generate plots for all valid pairs (limit to first 4 for 2x2 grid)
plot_list <- list()
for (i in 1:nrow(valid_pairs)) {
  p <- plot_pair(data, valid_pairs$base[i], valid_pairs$mock[i])
  if (!is.null(p)) {
    plot_list[[length(plot_list) + 1]] <- p
  }
  if (length(plot_list) >= 4) break
}

# Combine and save plots in 2x2 layout
if (length(plot_list) > 0) {
  # Ensure exactly 4 elements (pad with NULL if needed)
  if (length(plot_list) < 4) {
    plot_list <- c(plot_list, rep(list(NULL), 4 - length(plot_list)))
  }
  combined <- ggarrange(plotlist = plot_list, ncol = 2, nrow = 2)
  print(combined)

  ggsave("S_F14.png", plot = combined, width = 8, height = 8, dpi = 300, bg = "white")
  ggsave("S_F14.pdf", plot = combined, width = 8, height = 8, dpi = 300, bg = "white")
} else {
  cat("No valid plots generated.\n")
}
