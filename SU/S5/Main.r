 library(ggplot2)
library(dplyr)
library(ggsci)

# ---- Load Data ----
yfl_data <- read.table("PRJ.15.YFL.annotation.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ---- Clean and Filter ----
yfl_data$annotation[is.na(yfl_data$annotation)] <- "Unannotated"
yfl_data <- yfl_data %>% filter(status != "Target")

# ---- Get Top 10 Annotations ----
top10_annots <- table(yfl_data$annotation) %>%
  sort(decreasing = TRUE) %>%
  head(10)

# ---- Prepare Plot Data ----
plot_data <- data.frame(
  annotation = names(top10_annots),
  count = as.numeric(top10_annots)
)

# ---- Color Setup ----
nejm_color <- pal_nejm("default")(2)[2]
plot_data$fill <- ifelse(plot_data$annotation == "Unannotated", "grey50", nejm_color)

# ---- Plot ----
ggplot(plot_data, aes(x = reorder(annotation, count), y = count, fill = fill)) +
  geom_col(alpha = 0.85) +
  coord_flip() +
  scale_fill_identity() +
  labs(title = "Top 10 Annotations (Status ≠ Target)", x = "Annotation", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(10, 10, 10, 10)
  )
