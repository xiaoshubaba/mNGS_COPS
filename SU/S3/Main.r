library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(readr)

# ==== 1. read data ====
data <- read.table("Main.6.V2.txt", header = TRUE, sep = "\t")

# ==== 2. 3  rounds ====
rounds <- c(1, 2, 3)
completed_data <- data %>%
  select(prj, sample) %>% distinct() %>%
  crossing(round = rounds) %>%
  left_join(data, by = c("prj", "sample", "round")) %>%
  mutate(
    cov = ifelse(is.na(cov), 0, cov),
    signal = ifelse(is.na(signal), 0, signal)
  )

# ==== 3. round_label  ====
completed_data <- completed_data %>%
  mutate(round_label = factor(round, levels = 1:3,
                              labels = c("Co-assembly", "Forward Selection", "Reverse Selection")))

# ==== 4. mean ====
mean_summary <- completed_data %>%
  group_by(prj, round_label) %>%
  summarise(
    `Genomic Coverage (%)` = mean(cov),
    `Signal Intensity (%)` = mean(signal),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(`Genomic Coverage (%)`, `Signal Intensity (%)`),
               names_to = "metric", values_to = "value")

# ==== 5. statistics ====
compare_stats <- function(df, metric_col) {
  df %>%
    group_by(prj) %>%
    group_modify(~ {
      r1 <- .x %>% filter(round == 1)
      r2 <- .x %>% filter(round == 2)
      r3 <- .x %>% filter(round == 3)
      tibble(
        prj = .x$prj[1],
        p_2vs1 = if (nrow(r1) > 0 & nrow(r2) > 0) wilcox.test(r1[[metric_col]], r2[[metric_col]])$p.value else NA,
        p_3vs1 = if (nrow(r1) > 0 & nrow(r3) > 0) wilcox.test(r1[[metric_col]], r3[[metric_col]])$p.value else NA,
        effect_3vs1 = median(r3[[metric_col]]) - median(r1[[metric_col]])
      )
    }) %>%
    ungroup()
}

# ==== 6. statistics txt ====
cov_stats <- compare_stats(completed_data, "cov")
signal_stats <- compare_stats(completed_data, "signal")

# ==== 7. merge test ====
stats_all <- cov_stats %>%
  left_join(signal_stats, by = "prj", suffix = c("_cov", "_signal"))

# ==== 8. heatmap add ** ====
heatmap_data <- mean_summary %>%
  left_join(stats_all, by = "prj") %>%
  mutate(
    label_value = sprintf("%.1f", value),
    label = case_when(
      metric == "Genomic Coverage (%)" & round_label == "Forward Selection" & p_2vs1_cov < 0.05 ~ paste0(label_value, "*"),
      metric == "Genomic Coverage (%)" & round_label == "Reverse Selection" & p_3vs1_cov < 0.05 ~ paste0(label_value, "*"),
      metric == "Signal Intensity (%)" & round_label == "Forward Selection" & p_2vs1_signal < 0.05 ~ paste0(label_value, "*"),
      metric == "Signal Intensity (%)" & round_label == "Reverse Selection" & p_3vs1_signal < 0.05 ~ paste0(label_value, "*"),
      TRUE ~ label_value
    )
  )

# ==== 9.  n1 and n3  ====
n_by_round <- data %>%
  filter(round %in% c(1, 3)) %>%
  group_by(prj, round) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = round, values_from = n, names_prefix = "n")  # n1, n3

# ==== 10.  ====
heatmap_data <- heatmap_data %>%
  left_join(n_by_round, by = "prj") %>%
  mutate(project_with_n = paste0(prj, " (n1=", n1, "; n3=", n3, ")"))

# ==== 11. sort ====
sorted_projects <- stats_all %>%
  arrange(desc(effect_3vs1_cov)) %>%
  pull(prj)

heatmap_data$project_with_n <- factor(
  heatmap_data$project_with_n,
  levels = heatmap_data %>%
    filter(prj %in% sorted_projects) %>%
    distinct(prj, project_with_n) %>%
    mutate(order = match(prj, sorted_projects)) %>%
    arrange(order) %>%
    pull(project_with_n)
)

# ==== 12. vis ====
heatmap_data <- heatmap_data %>%
  mutate(
    prj = case_when(
      prj == "LA_HMB_2024" ~ "LA_AV_2024",
      prj == "LA_HRSV_2024" ~ "LA_PRIV_2024",
      TRUE ~ prj
    ),
    project_with_n = gsub("_", ".", project_with_n)
  )

p = ggplot(heatmap_data, aes(x = round_label, y = project_with_n, fill = value)) +
  geom_tile(color = "white", size = 0.3) +
  geom_text(aes(label = label), size = 3.5, fontface = "bold", color = "white") +
  scale_fill_gradientn(
    colours = c("#4393C3", "#92C5DE", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
    values = scales::rescale(c(0, 40, 60, 70, 80, 100)),
    limits = c(0, 100),
    name = "Value (%)",
    guide = guide_colorbar(
      barwidth = 0.8,
      barheight = 10,
      title.position = "top",
      title.hjust = 0.5,
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  facet_wrap(~metric, scales = "free_x",
             labeller = labeller(metric = c(
               "Genomic Coverage (%)" = "Genomic Coverage (%)",
               "Signal Intensity (%)" = "Signal Ratio (%)"  
             ))) +
  labs(
    x = "POPS Processing Step",
    y = "Project (Sample n1 and n3)",
    caption = "* p<0.05 vs. co-assembly step"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10, color = "black"),
    legend.text = element_text(size = 9, color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(face = "bold", size = 10, color = "black"),
    strip.text = element_text(face = "bold", size = 10, color = "black"),
    panel.grid = element_blank(),
    plot.caption = element_text(size = 8, hjust = 0, face = "italic", color = "black"),
    panel.spacing = unit(1, "lines"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank(),
    strip.background = element_blank()
  )
