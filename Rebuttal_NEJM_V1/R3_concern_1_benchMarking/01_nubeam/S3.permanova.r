library(vegan)

# Get all project directories (each folder = one project)
project_dirs <- list.dirs(".", recursive = FALSE)

# Initialize result table
results <- data.frame(Project = character(), Term = character(), R2 = numeric(), P = numeric(), stringsAsFactors = FALSE)

# Loop through projects
for (proj_dir in project_dirs) {
  proj_name <- basename(proj_dir)
  cohort_file <- file.path(proj_dir, paste0(proj_name, ".sample.redundant.cohort.txt"))
  matrix_file <- file.path(proj_dir, paste0(proj_name, ".sample.redundant.ditance.txt"))
  if (!file.exists(cohort_file) || !file.exists(matrix_file)) next
  metaData <- read.table(cohort_file, header = TRUE, stringsAsFactors = FALSE)
  DISTANCE <- as.dist(read.table(matrix_file, header = TRUE, row.names = 1, check.names = FALSE))
  # Match samples
  common_samples <- intersect(rownames(as.matrix(DISTANCE)), metaData$sample)
  if (length(common_samples) < 3) next
  metaData <- metaData[metaData$sample %in% common_samples, ]
  DISTANCE <- as.dist(as.matrix(DISTANCE)[common_samples, common_samples])
  # Replace cohort labeled as string "NA" with "Other"
  metaData$cohort[which(metaData$cohort == "NA")] <- "Other"
  # Drop real NA rows in cohort or sequencer
  metaData <- metaData[!is.na(metaData$cohort), ]
  if ("sequencer" %in% colnames(metaData)) {
    metaData <- metaData[!is.na(metaData$sequencer), ]
  }

  # Run adonis2
  if ("sequencer" %in% colnames(metaData) && length(unique(metaData$sequencer)) > 1) {
    obj <- adonis2(DISTANCE ~ cohort + sequencer, data = metaData, permutations = 2000)
    results <- rbind(results, 
                     data.frame(Project = proj_name, Term = "cohort", R2 = obj$R2[1], P = obj$Pr[1]),
                     data.frame(Project = proj_name, Term = "sequencer", R2 = obj$R2[2], P = obj$Pr[2]))
  } else {
    obj <- adonis2(DISTANCE ~ cohort, data = metaData, permutations = 2000)
    results <- rbind(results, 
                     data.frame(Project = proj_name, Term = "cohort", R2 = obj$R2[1], P = obj$Pr[1]))
  }
}

# Output
print(results)
write.table(results,"t3.txt",sep="\t",quote=F,row.names=FALSE)
