 library(dplyr)

# Get all project directories (each folder = one project)
project_dirs <- list.dirs(".", recursive = FALSE)

# Initialize result container
summary_list <- list()

for (proj_dir in project_dirs) {
  proj_name <- basename(proj_dir)
  cohort_file <- file.path(proj_dir, paste0(proj_name, ".sample.redundant.cohort.txt"))
  matrix_file <- file.path(proj_dir, paste0(proj_name, ".sample.redundant.ditance.txt"))

  if (!file.exists(cohort_file) || !file.exists(matrix_file)) next

  # Load cohort and distance matrix
  cohort <- read.table(cohort_file, header = TRUE, stringsAsFactors = FALSE)
  dist_matrix <- read.table(matrix_file, header = TRUE, row.names = 1, check.names = FALSE)

  # Filter cohort: only "Case" or "Control"
  cohort <- cohort[cohort$cohort %in% c("Case", "Control"), ]
  valid_samples <- cohort$sample

  # Filter dist_matrix to retain only valid samples (both rows and columns)
  dist_matrix <- dist_matrix[valid_samples, valid_samples]

  # Re-align sample names
  common_samples <- intersect(rownames(dist_matrix), cohort$sample)
  cohort <- cohort[cohort$sample %in% common_samples, ]
  dist_matrix <- dist_matrix[common_samples, common_samples]

  for (group in unique(cohort$cohort)) {
    group_samples <- cohort$sample[cohort$cohort == group]
    raw_n <- length(group_samples)

    # Filter: ≥3 distances < 1.0 (excluding self and distance = 1.0)
    filtered_samples <- group_samples[sapply(group_samples, function(s) {
      dists <- dist_matrix[s, -which(colnames(dist_matrix) == s)]
      sum(dists < 1.0 & dists > 0, na.rm = TRUE) >= 3
    })]

    if (length(filtered_samples) == 0) {
      summary_list[[paste0(proj_name, "_", group)]] <- data.frame(
        Dataset = proj_name,
        Cohort = group,
        Raw_Sample_Number = raw_n,
        Filtered_Sample_Number = 0,
        Correct_Predictions = 0,
        Accuracy = NA
      )
      next
    }

    # Classify based on top 3 neighbors
    classification <- sapply(filtered_samples, function(s) {
      neighbors <- sort(dist_matrix[s, -which(colnames(dist_matrix) == s)])
      top3 <- names(neighbors)[1:3]
      top3_cohorts <- cohort$cohort[match(top3, cohort$sample)]
      sum(top3_cohorts == "Case") >= 2
    })

    # Get true labels
    true_labels <- cohort$cohort[match(filtered_samples, cohort$sample)]
    correct_preds <- sum(classification == (true_labels == "Case"))
    total_preds <- length(filtered_samples)

    summary_list[[paste0(proj_name, "_", group)]] <- data.frame(
      Dataset = proj_name,
      Cohort = group,
      Raw_Sample_Number = raw_n,
      Filtered_Sample_Number = total_preds,
      Correct_Predictions = correct_preds,
      Accuracy = round(correct_preds / total_preds, 4)
    )
  }
}

# Combine all summaries into a single table
final_summary <- bind_rows(summary_list)
print(final_summary)
