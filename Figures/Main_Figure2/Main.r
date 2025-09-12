# === Three-panel figure with 0–100% y-axes everywhere ===
# p1: Reads vs Signal
# p2: Pooling Sample Number Optimization
# p3: Filtering Parameter Optimization
# Only p3 shows the legend; panels tagged A/B/C

library(tidyverse)
library(scales)
library(ggsci)
library(patchwork)
library(readr)
library(mgcv)

# ----- unified palette/labels -----
metric_colors <- pal_nejm("default")(4)
names(metric_colors) <- c("Genomic coverage","Signal purity","Success rate","Final score")
metric_levels <- names(metric_colors)
title_small <- theme(plot.title = element_text(size = 10, face = "bold"))

# ===================== p1 =====================
main <- read_tsv("Table.5.single.F3.V3.txt")

p1 <- ggplot(main, aes(x = gatkReads)) +
  geom_point(aes(y = asmCov, color = "Genomic coverage"), alpha = 0.6, size = 1) +
  geom_point(aes(y = Signal,  color = "Signal purity"),    alpha = 0.6, size = 1) +
  geom_smooth(aes(y = asmCov, color = "Genomic coverage"), method = "gam", se = TRUE, show.legend = FALSE) +
  geom_smooth(aes(y = Signal, color = "Signal purity"),    method = "gam", se = TRUE, show.legend = FALSE) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = metric_colors, breaks = metric_levels, name = NULL) +
  # Y in 0–100% with percent labels
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20),
                     labels = function(x) paste0(x, "%")) +
  labs(title = "Reads vs Signal",
       x = "Detected pathogen reads (log10)",
       y = "Genomic coverage / Signal purity (%)") +
  theme_minimal(base_size = 11) +
  title_small +
  theme(legend.position = "none")

# ===================== p2 =====================
mix <- read_tsv("Table.5.mix.F3.txt") %>%
  mutate(Mixnumber = as.integer(Mixnumber),
         success   = as.integer(cov > 0))

score_df <- mix %>%
  group_by(Mixnumber) %>%
  summarise(mean_cov = mean(cov),
            mean_signal = mean(signal),
            success_rate = mean(success), .groups = "drop") %>%
  mutate(norm_cov = rescale(mean_cov),
         norm_signal = rescale(mean_signal),
         norm_success = rescale(success_rate),
         score = (norm_cov + norm_signal + norm_success) / 3) %>%
  # convert to percent for plotting 0–100%
  mutate(across(c(norm_cov, norm_signal, norm_success, score), ~ .x * 100))

plot_df2 <- score_df %>%
  select(Mixnumber, norm_cov, norm_signal, norm_success, score) %>%
  pivot_longer(-Mixnumber, names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric,
                         levels = c("norm_cov","norm_signal","norm_success","score"),
                         labels = c("Genomic coverage","Signal purity","Success rate","Final score")))

p2 <- ggplot(plot_df2, aes(Mixnumber, Value, color = Metric)) +
  geom_line(size = 1.1, alpha = 0.8) +
  geom_point(size = 2,   alpha = 0.8) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "red", linewidth = 0.6) +
  scale_color_manual(values = metric_colors, breaks = metric_levels, name = NULL) +
  scale_x_continuous(breaks = c(1,3,5,7,9)) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20),
                     labels = function(x) paste0(x, "%")) +
  labs(title = "Pooling Sample Number Optimization",
       x = "Number of pooling samples",
       y = "Normalized value / score (%)") +
  theme_minimal(base_size = 11) +
  title_small +
  theme(legend.position = "none")

# ===================== p3 =====================
grid <- read_table("global_param_scored.txt") %>%
  arrange(case.cutoff, control.cutoff) %>%
  mutate(combo = factor(paste0("X", case.cutoff, "Y", control.cutoff),
                        levels = unique(paste0("X", case.cutoff, "Y", control.cutoff))))
best_row   <- which.max(grid$tradeoff_score)
best_combo <- as.character(grid$combo[best_row])
best_x     <- which(levels(grid$combo) == best_combo)

plot_df3 <- grid %>%
  mutate(across(c(asmCov_norm, Signal_norm, success_norm, tradeoff_score), ~ .x * 100)) %>%
  select(combo, asmCov_norm, Signal_norm, success_norm, tradeoff_score) %>%
  pivot_longer(-combo, names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         asmCov_norm   = "Genomic coverage",
                         Signal_norm   = "Signal purity",
                         success_norm  = "Success rate",
                         tradeoff_score= "Final score"),
         Metric = factor(Metric, levels = metric_levels))

p3 <- ggplot(plot_df3, aes(combo, Value, color = Metric, group = Metric)) +
  geom_line(size = 1) +
  geom_vline(xintercept = best_x, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 105, hjust = 0,
           label = paste0("Best cutoff: ", best_combo), size = 3.5) +
  scale_color_manual(values = metric_colors, breaks = metric_levels, name = NULL) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20),
                     labels = function(x) paste0(x, "%"),
                     name = "Normalized value / score (%)") +
  labs(title = "Filtering Parameter Optimization",
       x = "Cutoff combination (X = case, Y = control)") +
  theme_minimal(base_size = 11) +
  title_small +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

# ===================== assemble =====================
final_plot <- (p1 + p2) / p3 +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

print(final_plot)
# ggsave("POPS_three_panel.png", final_plot, width = 11, height = 7, dpi = 300)

