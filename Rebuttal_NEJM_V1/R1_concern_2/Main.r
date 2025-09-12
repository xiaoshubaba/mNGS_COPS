library(ggplot2)
library(dplyr)
library(lubridate) 
colorsQ = c("#BC3C29FF" ,"#0072B5FF" ,"#E18727FF")
df <- read.table("Contamination.txt",head=T)
df$date <- as.Date(as.character(df$date), format = "%Y%m%d")
df <- df %>%
  mutate(month = floor_date(date, unit = "month")) %>% # Create monthly period column
  group_by(cohort, month) %>%
  summarise(count = n(), .groups = 'drop')
stage_date <- as.Date("20221201", format = "%Y%m%d")
df$cohort = factor(df$cohort,levels=c("COVID","Control","Contamination"))
plot <- ggplot(df, aes(x = month, y = count, fill = cohort)) +
  geom_bar(stat = "identity", position = "stack", width = 20) + # Adjust width for better visibility
  scale_fill_manual(values = colorsQ) +  # Specify colors for cohorts
  theme_minimal() +  # Use minimal theme for a clean look
  labs(x = "Month", y = "Number of Samples", title = "Distribution of Samples by Cohort") +
  geom_vline(xintercept = as.numeric(stage_date), linetype = "dashed", color = "grey", size = 0.5) +  # Vertical line for stage division
  annotate("text", x = stage_date, y = max(df$count) - 1, label = "Stage 2", color = "red", hjust = -0.1, size = 3) +  # Annotation for stage 2
  annotate("text", x = stage_date, y = max(df$count) - 5, label = "Stage 1", color = "red", hjust = 1.1, size = 3) +  # Annotation for stage 1
  scale_x_date(date_breaks = "1 month", date_labels = "%Y-%m") +  # Date format for x-axis
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility


  # Load required libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Load your data
df <- read.table("SARS2.old.txt", header=TRUE, sep="\t")

# Summarize coverage and signal across POPS steps (round)
df_summary <- df %>%
  group_by(round) %>%
  summarise(mean_cov = mean(cov, na.rm = TRUE),
            mean_signal = mean(signal, na.rm = TRUE))

# Melt the data for plotting
df_melt <- melt(df, id.vars = c("sample", "prj", "mixNum", "round"),
                measure.vars = c("cov", "signal"),
                variable.name = "Metric", value.name = "Value")

# Plot boxplots for coverage and signal purity across rounds
p2 = ggplot(df_melt, aes(x = factor(round), y = Value, fill = Metric)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8) +
  facet_wrap(~Metric, scales = "free_y", ncol = 1, labeller = labeller(Metric = c(cov = "Genomic Coverage", signal = "Signal Purity"))) +
  labs(title = "Changes in Genomic Coverage and Signal Purity across POPS Steps",
       x = "POPS Step (Round)",
       y = "Value") +
  scale_fill_manual(values = c("cov" = "#6CA6CD", "signal" = "#FF7F50")) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "gray90"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

