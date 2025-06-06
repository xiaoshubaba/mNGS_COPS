library(ggplot2)
library(dplyr)
library(ggsci)
library(patchwork)

# ---- Read Data ----
hiv_data <- read.table("PRJ.04.HIV.annotation.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
zika_data <- read.table("PRJ.05.zika.annotation.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ---- Get Top 10 Annotations by Frequency ----
get_top10 <- function(df) {
  df %>%
    group_by(annotation, status) %>%
    summarise(total_freq = sum(frquency, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_freq)) %>%
    slice_max(total_freq, n = 10) %>%
    mutate(annotation = reorder(annotation, total_freq))
}

hiv_top10 <- get_top10(hiv_data)
zika_top10 <- get_top10(zika_data)

# ---- Define Colors ----
nejm_colors <- pal_nejm("default")(2)
fill_colors <- c("Annotated" = nejm_colors[2], "Unannotated" = "grey50")

# ---- Plot Function ----
make_plot <- function(df, title) {
  ggplot(df, aes(x = total_freq, y = annotation, fill = status)) +
    geom_col(alpha = 0.85) +
    scale_fill_manual(values = fill_colors) +
    labs(title = title, x = "Fragment Frequency", y = "Annotation") +
    theme_minimal(base_size = 12) +
    theme(legend.position="none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.y = element_text(margin = margin(r = 10)),
      legend.title = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ---- Create Plots ----
plot_hiv <- make_plot(hiv_top10, "HIV")
plot_zika <- make_plot(zika_top10, "Zika") +
  theme(axis.title.y = element_blank())

# ---- Combine ----
combined_plot <- plot_hiv + plot_zika +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "A")
