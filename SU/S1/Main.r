library(ggplot2)
library(dplyr)

# ---- Read tab-delimited file safely ----
jd_data <- read.table("jd.list.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = NULL)

# ---- Rename projects for clarity ----
jd_data <- jd_data %>%
  mutate(prj = recode(prj,
                      "LA_HRSV_2024" = "LA_PRIV_2024",
                      "LA_HMB_2024" = "LA_Adenovirus_2024"))

# ---- Plot: Jaccard distance by project ----
ggplot(jd_data, aes(x = prj, y = jd)) +
  geom_boxplot(outlier.shape = NA, fill = "#4E79A7", alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  labs(x = "Project", y = "Jaccard Distance", title = "Jaccard Distance Across Projects") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
mean_jd <- jd_data %>%
  group_by(prj) %>%
  summarise(mean_jd = mean(jd, na.rm = TRUE)) %>%
  arrange(desc(mean_jd))

# ---- Print the result ----
print(mean_jd)
