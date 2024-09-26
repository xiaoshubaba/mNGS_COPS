library(ggplot2)
library(tidyr)
library(ggsci)
colorsQ = c(pal_nejm("default")(5),"grey50")
#visual
df <- read.table("COVID.patient.n219.singleAsm.status.m200.txt",head=T)
df_long <- pivot_longer(df, cols = c(fragment_num, mapped_fragment_num, mapped_ref_bases, BreakPoint, cov, snp, indel, signal_ratio), names_to = "feature", values_to = "value")
ggplot(df_long, aes(x = type, y = value, color = type)) +
  geom_boxplot(outlier.shape = NA) + scale_color_manual(values= colorsQ)+ # Do not show outliers
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +  # Add jitter to show individual points
  facet_wrap(~feature, scales = "free_y") +  # Faceting by feature, with free y scales
  labs(title = "Comparison of Genome Assembly Methods Across Features",
       x = NULL,  # Remove x-axis label
       y = "Value",
       color = "Method") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis names
        axis.ticks.x = element_blank())  # Remove x-axis names
#comparision
library(tidyverse)

# Assuming df is already loaded into the environment
# Example: df <- read.csv("your_data.csv")

# Filtering out pairs of methods to compare
unique_methods <- unique(df$type)
comparisons <- combn(unique_methods, 2, simplify = FALSE)

# Function to perform Wilcoxon rank-sum tests and return formatted results
perform_wilcox_test <- function(data, feature, pair) {
  method1_data <- data[data$type == pair[1], ]
  method2_data <- data[data$type == pair[2], ]

  wilcox_result <- wilcox.test(method1_data[[feature]], method2_data[[feature]],
                               alternative = "greater")

  # Conditionally determine which method is better based on p-value
  if (wilcox_result$p.value < 0.1) {

    better_method <- ifelse(median(method1_data[[feature]]) > median(method2_data[[feature]]), pair[1], pair[2])
  } else {
    better_method <- "No significant difference"
  }

  # Return a data frame with the results
  data.frame(
    Method1 = pair[1],
    Method2 = pair[2],
    Feature = feature,
    P_Value = wilcox_result$p.value,
    Better_Method = better_method
  )
}

# Features to analyze
features <- c("cov", "signal_ratio", "BreakPoint")

# Applying the function to each feature and binding rows
results <- do.call(rbind, lapply(features, function(feature) {
  do.call(rbind, lapply(comparisons, function(pair) perform_wilcox_test(df, feature, pair)))
}))

# Print the combined results
print(results)
