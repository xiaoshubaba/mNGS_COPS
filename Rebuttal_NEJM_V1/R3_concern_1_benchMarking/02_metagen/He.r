# Load required package
if (!requireNamespace("entropy", quietly = TRUE)) install.packages("entropy", repos = "https://cloud.r-project.org")
library(entropy)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript h.r input_file output_file")
}

input_file <- args[1]
output_file <- args[2]

# Load data
TAB <- read.table(input_file, header = TRUE, stringsAsFactors = FALSE)

# Filter to pathogen fragments (annotation ≠ "N")
TAB_pathogen <- TAB[TAB$annotation != "N", ]

# Compute homogeneity score
if (nrow(TAB_pathogen) == 0) {
  homogeneity <- NA
} else {
  cluster_counts <- table(TAB_pathogen$cluster)
  if (length(cluster_counts) <= 1) {
    homogeneity <- 1
  } else {
    H_C <- entropy.empirical(cluster_counts, unit = "log2")
    H_C_max <- log2(length(cluster_counts))
    homogeneity <- ifelse(H_C_max == 0, 1, 1 - (H_C / H_C_max))
  }
}

# Write numeric result only
write(homogeneity, file = output_file)
