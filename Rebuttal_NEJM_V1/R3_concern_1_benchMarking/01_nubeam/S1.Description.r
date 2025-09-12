library(dplyr)
library(ggplot2)
library(ggpubr)

# Load all files
files <- list.files(pattern = "\\.PairSingle\\.distance\\.txt$", recursive = TRUE, full.names = TRUE)
result_list <- list()
plot_data <- data.frame()

dataset_colors <- c(
  "CHIKV" = "#1f77b4",
  "Ebola" = "#ff7f0e",
  "LASV" = "#2ca02c",
  "SARS02" = "#d62728",
  "SARS01" = "#9467bd"
)

for (file in files) {
  df <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  if (!"distance" %in% colnames(df)) next
  label <- basename(dirname(file))
  label <- gsub("[^a-zA-Z0-9\\-]", "_", label)
  total <- nrow(df)
  full_dist <- sum(df$distance == 1.0, na.rm = TRUE)
  result_list[[label]] <- data.frame(
    Dataset = label,
    Total = total,
    Distance_1.0 = full_dist,
    Proportion = round(full_dist / total, 4)
  )
  df$Dataset <- label
  plot_data <- bind_rows(plot_data, df)
}

summary_df <- bind_rows(result_list)
plot_data$PairType <- ifelse(plot_data$PAIR == "Y", "Intra-sample (PAIR)", "Inter-sample")
plot_data$Dataset <- as.character(plot_data$Dataset)

# --- p1: Density plot ---
annotation_text <- paste(
  apply(summary_df, 1, function(row) {
    paste0(row["Dataset"], ": ", row["Distance_1.0"], "/", row["Total"],
           " (", as.numeric(row["Proportion"])*100, "%)")
  }),
  collapse = "\n"
)

plot1 <- ggplot(plot_data, aes(x = distance, fill = Dataset, color = Dataset)) +
  geom_density(alpha = 0.3, adjust = 1.2) +
  scale_fill_manual(values = dataset_colors) +
  scale_color_manual(values = dataset_colors) +
  labs(
    title = "A:Distribution of Pairwise Hellinger Distances",
    x = "Distance", y = "Density"
  ) +
  annotate("text", x = 0.2, y = Inf, label = annotation_text, hjust = 0, vjust = 1.5, size = 2.5) +
  theme_minimal(base_size = 14)

# --- p2: Intra vs Inter-sample comparison (no p-values shown) ---

	# ---- Add p-values to plot2 facet titles ----
# Calculate per-dataset Wilcoxon p-values
p2_labels <- plot_data %>%
  group_by(Dataset) %>%
  filter(PairType %in% c("Intra-sample (PAIR)", "Inter-sample")) %>%
  summarise(
    p_value = tryCatch({
      wilcox.test(distance[PairType == "Intra-sample (PAIR)"], distance[PairType == "Inter-sample"])$p.value
    }, error = function(e) NA)
  ) %>%
  mutate(
    label = paste0(Dataset, "(p=", formatC(p_value, format = "e", digits = 2), ")")
  )

# Relabel facet titles using labeller
label_map <- setNames(p2_labels$label, p2_labels$Dataset)
custom_labeller <- as_labeller(label_map)

# --- Updated plot2 with p-values in titles ---
plot2 <- ggboxplot(
  plot_data,
  x = "PairType",
  y = "distance",
  color = "PairType",
  palette = c("Intra-sample (PAIR)" = "#1f77b4", "Inter-sample" = "#ff7f0e"),
  ylab = "Hellinger Distance",
  xlab = ""
) +
  facet_wrap(~ Dataset, labeller = custom_labeller) +
  labs(title = "B:Intra- vs. Inter-sample Distance") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.98, 0.02),
    legend.justification = c("right", "bottom"),
    legend.title = element_blank()
  )
# --- p3: Inter-sample distances by cohort (no p-values) ---
plot_data_subset <- plot_data %>% filter(PAIR == "N")

plot3 <- ggboxplot(
  plot_data_subset,
  x = "Cohort",
  y = "distance",
  color = "Cohort",
  ylab = "Hellinger Distance",
  xlab = ""
) +
  facet_wrap(~ Dataset) +
  labs(title = "C:Inter-sample Distance by Cohort") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.98, 0.02),
    legend.justification = c("right", "bottom"),
    legend.title = element_blank()
  )

# --- Pairwise test for p3 and final result table ---
get_enrichment_label <- function(group1, group2, label1, label2, pval) {
  if (pval >= 0.05) return("NA")
  m1 <- mean(group1)
  m2 <- mean(group2)
  if (m1 > m2) return(paste("Enriched in", label1))
  else if (m1 < m2) return(paste("Enriched in", label2))
  else return("NA")
}

pairwise_results <- data.frame()

for (dataset in unique(plot_data_subset$Dataset)) {
  df <- plot_data_subset %>% filter(Dataset == dataset)
  if (all(c("Between", "Within_Case", "Within_Control") %in% unique(df$Cohort))) {
    g_between <- df %>% filter(Cohort == "Between") %>% pull(distance)
    g_case <- df %>% filter(Cohort == "Within_Case") %>% pull(distance)
    g_control <- df %>% filter(Cohort == "Within_Control") %>% pull(distance)
    p1 <- wilcox.test(g_between, g_case)$p.value
    p2 <- wilcox.test(g_between, g_control)$p.value
    result_row <- data.frame(
      Dataset = dataset,
      p_Between_vs_Within_Case = formatC(p1, format = "e", digits = 2),
      Direction_Between_vs_Within_Case = get_enrichment_label(g_between, g_case, "Between", "Within_Case", p1),
      p_Between_vs_Within_Control = formatC(p2, format = "e", digits = 2),
      Direction_Between_vs_Within_Control = get_enrichment_label(g_between, g_control, "Between", "Within_Control", p2)
    )
    pairwise_results <- rbind(pairwise_results, result_row)
  }
}

# Final p3 result table
print(pairwise_results) # pvalue added to plot2 titles
