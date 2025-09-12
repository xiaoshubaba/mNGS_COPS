library(ggsci)
library(ggpubr)
Main <- read.table("YFL.Local_sensitive.V2.table.txt",head=T,sep="\t")
Main$log_meanDepth = log(Main$length)
#right
colorsQ = c("grey50",pal_nejm("default")(1))
sr = round(sum(Main[which(Main$annotation=="Y"),]$length)/sum(Main$length)*100,digits=2)
f1_left = ggscatter(Main,"gc","log_meanDepth",xlab="GC content",color="annotation",alpha=I(0.9),legend="none",palette=colorsQ,size="length",title=paste0("contigs:",dim(Main)[1], ", signal ratio:",sr,"%"))

Main_m2 = Main[which(Main$case.fre>=2),]
sr = round(sum(Main_m2[which(Main_m2$annotation=="Y"),]$length)/sum(Main_m2$length)*100,digits=2)
colorsQ_tem = c()
for (j in 1:dim(table(Main_m2$annotation))){
         if (table(Main_m2$annotation)[j] > 0){
         colorsQ_tem = c(colorsQ_tem,colorsQ[j])
         }
}

f1_middle = ggscatter(Main_m2,"gc",xlab="GC content","log_meanDepth",color="annotation",alpha=I(0.9),legend="none",palette=colorsQ_tem,size="length",title=paste0("contigs:",dim(Main_m2)[1],", signal ratio:",sr,"%"))
#
Main_m2_m2 = Main_m2[which(Main_m2$control.fre<=1),]
colorsQ_tem = c()
for (j in 1:dim(table(Main_m2_m2$annotation))){
         if (table(Main_m2_m2$annotation)[j] > 0){
         colorsQ_tem = c(colorsQ_tem,colorsQ[j])
         }
}
sr = round(sum(Main_m2_m2[which(Main_m2_m2$annotation=="Y"),]$length)/sum(Main_m2_m2$length)*100,digits=2)
f1_right = ggscatter(Main_m2_m2,"gc",xlab="GC content","log_meanDepth",color="annotation",alpha=I(0.9),legend="none",palette=colorsQ_tem,size="length",title=paste0("contigs:",dim(Main_m2_m2)[1], ", signal ratio:",sr,"%"))

########################
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

# Read data file
data <- read.table("YFL.k141_396.false.alignment.txt", 
                  header = TRUE, 
                  fill = TRUE, 
                  stringsAsFactors = FALSE)

# Sort data: first by fragment, then by sample
data_sorted <- data %>%
  arrange(fragment, sample)

# Preprocess data with correct sorting
data_processed <- data_sorted %>%
  mutate(
    alignment = ifelse(is.na(alignment), "", as.character(alignment)),
    y = rev(row_number())
  ) %>%
  group_by(sample) %>%
  mutate(
    group_size = n(),
    group_y_min = min(y) - 0.4,
    group_y_max = max(y) + 0.4
  ) %>%
  ungroup()

# Parse CIGAR string
parse_cigar <- function(cigar, start) {
  if (cigar == "" || is.na(cigar)) return(data.frame())
  ops <- str_match_all(cigar, "(\\d+)([A-Z])")[[1]]
  if (nrow(ops) == 0) return(data.frame())
  df <- data.frame(
    length = as.numeric(ops[, 2]),
    type = ops[, 3],
    stringsAsFactors = FALSE
  )
  df$ref_start <- cumsum(c(start, df$length))[1:nrow(df)]
  df$ref_end <- df$ref_start + df$length - 1
  df$ref_start[1] <- start
  df$ref_end[1] <- start + df$length[1] - 1
  return(df)
}

# Create alignment data for all reads
alignments <- lapply(1:nrow(data_processed), function(i) {
  df <- parse_cigar(data_processed$alignment[i], data_processed$start[i])
  if (nrow(df) > 0) {
    df$y <- data_processed$y[i]
    df$sample <- data_processed$sample[i]
  }
  return(df)
}) %>% bind_rows()

# Color scheme with light colors
operation_colors <- c(
  "S" = "#FF9999",  # Soft Clip
  "M" = "#99CCFF",  # Match
  "D" = "#99FF99"   # Deletion
)

# Calculate plot range
x_min <- min(data_processed$start, na.rm = TRUE) - 50
x_max <- max(alignments$ref_end, na.rm = TRUE) + 50

# Create group connection lines (only for groups with 2+ sequences)
group_lines <- data_processed %>%
  distinct(sample, group_y_min, group_y_max, group_size) %>%
  filter(group_size >= 2) %>%
  mutate(
    group_x = x_min - 50,
    group_mid = (group_y_min + group_y_max)/2
  )

# Reference baseline
base_line <- data.frame(
  x = c(x_min, x_max),
  y = c(0, 0)
)

# Create position ticks
position_ticks <- data.frame(
  x = seq(floor(x_min/50)*50, ceiling(x_max/50)*50, by = 50)
)

# Ensure correct legend order
legend_order <- c("S", "M", "D")
alignments$type <- factor(alignments$type, levels = legend_order)

# Create the alignment plot
pd <- ggplot() +
  # Gray reference baseline
  geom_line(data = base_line, aes(x = x, y = y), color = "gray70", size = 1) +
  # Sample connection lines (for groups with 2+ sequences)
  geom_segment(
    data = group_lines,
    aes(x = group_x, xend = group_x, y = group_y_min, yend = group_y_max),
    color = "gray70", size = 0.5
  ) +
  # Add group size labels
  geom_text(
    data = group_lines,
    aes(x = group_x - 2, y = group_mid, label = paste0("n=", group_size)),
    size = 3, color = "gray30", hjust = 1
  ) +
  # Alignment segments
  geom_segment(
    data = alignments,
    aes(x = ref_start, xend = ref_end, y = y, yend = y, color = type),
    size = 0.8, lineend = "butt"
  ) +
  # Position ticks
  geom_segment(
    data = position_ticks,
    aes(x = x, xend = x, y = 0, yend = -0.15),
    color = "gray50"
  ) +
  # Position labels
  geom_text(
    data = position_ticks,
    aes(x = x, y = -0.4, label = x),
    size = 2.5, color = "gray40"
  ) +
  # Color scale with correct legend order
  scale_color_manual(
    name = "Operation",
    values = operation_colors,
    labels = c("S: Soft Clip", "M: Match", "D: Deletion"),
    drop = FALSE
  ) +
  # Axis labels and title
  labs(
    x = "Reference Position",
    y = ""  ) +
  # Theme settings
  theme_minimal() +
  theme(
    panel.grid = element_blank(),    axis.text.x = element_blank(),

    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.8, "lines")
  ) +
  # Coordinate range
  coord_cartesian(
    xlim = c(x_min - 60, x_max),
    ylim = c(-0.7, max(data_processed$y) + 0.5)
)
library(dplyr)
library(ggplot2)
library(readr)

alignment_data <- read.table("R3.4.data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

fragment_counts <- alignment_data %>%
  group_by(fragment) %>%
  summarise(read_count = n())

top_fragments <- fragment_counts %>%
  filter(read_count > 10)

subset_data <- alignment_data %>%
  inner_join(top_fragments, by = "fragment") %>%
  mutate(fragment_label = paste0(fragment, " (n=", read_count, ")"))

pe <- ggplot(subset_data, aes(x = fragment_label, y = start)) +
  geom_jitter(width = 0.3, height = 0, size = 1.2, alpha = 0.7, color = "steelblue") +
  theme_minimal() +
  labs(
    x = "Fragment (with read count)",
    y = "Start Position"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_blank()
  )

# Sort data: first by fragment, then by sample
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
alignment_data <- read.table("R3.4.data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data <- alignment_data[which(alignment_data$fragment=="k141_80"),]
data_sorted <- data %>%
  arrange(fragment, sample)

# Preprocess data with correct sorting
data_processed <- data_sorted %>%
  mutate(
    alignment = ifelse(is.na(alignment), "", as.character(alignment)),
    y = rev(row_number())
  ) %>%
  group_by(sample) %>%
  mutate(
    group_size = n(),
    group_y_min = min(y) - 0.4,
    group_y_max = max(y) + 0.4
  ) %>%
  ungroup()

# Parse CIGAR string
parse_cigar <- function(cigar, start) {
  if (cigar == "" || is.na(cigar)) return(data.frame())
  ops <- str_match_all(cigar, "(\\d+)([A-Z])")[[1]]
  if (nrow(ops) == 0) return(data.frame())
  df <- data.frame(
    length = as.numeric(ops[, 2]),
    type = ops[, 3],
    stringsAsFactors = FALSE
  )
  df$ref_start <- cumsum(c(start, df$length))[1:nrow(df)]
  df$ref_end <- df$ref_start + df$length - 1
  df$ref_start[1] <- start
  df$ref_end[1] <- start + df$length[1] - 1
  return(df)
}

# Create alignment data for all reads
alignments <- lapply(1:nrow(data_processed), function(i) {
  df <- parse_cigar(data_processed$alignment[i], data_processed$start[i])
  if (nrow(df) > 0) {
    df$y <- data_processed$y[i]
    df$sample <- data_processed$sample[i]
  }
  return(df)
}) %>% bind_rows()

# Color scheme with light colors
operation_colors <- c(
  "S" = "#FF9999",  # Soft Clip
  "M" = "#99CCFF",  # Match
  "D" = "#99FF99"   # Deletion
)

# Calculate plot range
x_min <- min(data_processed$start, na.rm = TRUE) - 50
x_max <- max(alignments$ref_end, na.rm = TRUE) + 50

# Create group connection lines (only for groups with 2+ sequences)
group_lines <- data_processed %>%
  distinct(sample, group_y_min, group_y_max, group_size) %>%
  filter(group_size >= 2) %>%
  mutate(
    group_x = x_min - 50,
    group_mid = (group_y_min + group_y_max)/2
  )

# Reference baseline
base_line <- data.frame(
  x = c(x_min, x_max),
  y = c(0, 0)
)

# Create position ticks
position_ticks <- data.frame(
  x = seq(floor(x_min/50)*50, ceiling(x_max/50)*50, by = 50)
)

# Ensure correct legend order
legend_order <- c("S", "M", "D")
alignments$type <- factor(alignments$type, levels = legend_order)

# Create the alignment plot
pf <- ggplot() +
  # Gray reference baseline
  geom_line(data = base_line, aes(x = x, y = y), color = "gray70", size = 1) +
  # Sample connection lines (for groups with 2+ sequences)
  geom_segment(
    data = group_lines,
    aes(x = group_x, xend = group_x, y = group_y_min, yend = group_y_max),
    color = "gray70", size = 0.5
  ) +
  # Add group size labels
  geom_text(
    data = group_lines,
    aes(x = group_x - 2, y = group_mid, label = paste0("n=", group_size)),
    size = 3, color = "gray30", hjust = 1
  ) +
  # Alignment segments
  geom_segment(
    data = alignments,
    aes(x = ref_start, xend = ref_end, y = y, yend = y, color = type),
    size = 0.8, lineend = "butt"
  ) +
  # Position ticks
  geom_segment(
    data = position_ticks,
    aes(x = x, xend = x, y = 0, yend = -0.15),
    color = "gray50"
  ) +
  # Position labels
  geom_text(
    data = position_ticks,
    aes(x = x, y = -0.4, label = x),
    size = 2.5, color = "gray40"
  ) +
  # Color scale with correct legend order
  scale_color_manual(
    name = "Operation",
    values = operation_colors,
    labels = c("S: Soft Clip", "M: Match", "D: Deletion"),
    drop = FALSE
  ) +
  # Axis labels and title
  labs(
    x = "Reference Position",
    y = ""  ) +
  # Theme settings
  theme_minimal() +
  theme(
    panel.grid = element_blank(),    axis.text.x = element_blank(),

    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.8, "lines")
  ) +
  # Coordinate range
  coord_cartesian(
    xlim = c(x_min - 60, x_max),
    ylim = c(-0.7, max(data_processed$y) + 0.5)
)
