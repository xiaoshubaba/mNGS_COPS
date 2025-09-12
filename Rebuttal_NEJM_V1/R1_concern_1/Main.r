library(tidyverse)
library(scales)
library(ggsci)

# Files + simulated totals
sim_map <- tibble::tribble(
  ~virus,         ~path,                                        ~N_sim,
  "SARS_cov_2",   "SARS_cov_2.mimic.escape.txt",      2985400,
  "HIV",          "HIV.mimic.escape.txt",              902300,
  "IVA",  "IVA.mimic.escape.txt",             1310600
) %>% filter(file.exists(path))

read_escape <- function(virus, path, N_sim) {
  readr::read_tsv(
    path,
    col_types = readr::cols(
      virus = readr::col_character(),
      length = readr::col_integer(),
      Expect_errors = readr::col_integer(),
      Actual_error = readr::col_double(),
      expected_Escape_blast = readr::col_double(),
      expected_Escape_bowtie = readr::col_double()
    )
  ) %>%
    transmute(
      virus = virus,                       # enforce label from sim_map
      length_bp = length,
      expect_errors = Expect_errors,
      miss_hits = expected_Escape_blast,   # absolute "miss" hits (BLAST)
      N_sim = N_sim
    ) %>%
    group_by(virus, length_bp, expect_errors) %>%
    summarise(miss_hits = sum(miss_hits), N_sim = dplyr::first(N_sim), .groups = "drop") %>%
    mutate(miss_ratio = miss_hits / N_sim)
}

# Build data
escape_df <- purrr::pmap_dfr(sim_map, read_escape) %>%
  mutate(
    virus = factor(virus, levels = c("SARS_cov_2", "HIV", "IVA")),
    mut_frac = expect_errors / length_bp,         # x = #mutations / read length
    detect_raw = 1 - miss_ratio,                  # y = detection = 1 - escape
    length_bp = factor(length_bp, levels = c(50, 150))
  )

# NEJM palette (3 colors, fixed mapping to virus levels)
nejm_pal <- ggsci::pal_nejm("default")(3)
names(nejm_pal) <- levels(escape_df$virus)

shapes <- c("50" = 16, "150" = 17)

# Plot
leg_opts <- theme(
  legend.position = c(0.98, 0.02),        # (x, y) in NPC coordinates
  legend.justification = c(1, 0),         # anchor bottom-right
  legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
  legend.key.size = unit(8, "pt"),
  legend.title = element_text(size = 9),
  legend.text  = element_text(size = 8)
)
#
p0 <- ggplot(
  escape_df,
  aes(x = mut_frac, y = detect_raw,
      color = virus, shape = length_bp,
      group = interaction(virus, length_bp))
) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  geom_point(size = 2, alpha = 0.95) +
  scale_color_manual(values = nejm_pal, name = "Virus") +
  scale_shape_manual(values = shapes, name = "Read length (bp)") +
  scale_x_continuous(labels = label_percent(accuracy = 1)) +
  scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1)) +
  labs(
    x = "Mutation fraction per read",
    y = "Missing ratio"
  ) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())
p0 <- p0 + leg_opts          # keep legend here



#
# ---- Part 1 analysis: time-limited reference impact ----
# Inputs: /mnt/data/mis.ratio.txt
# Final output: a single ggarrange figure with p1 (distributions) + p2 (scatter),
#               using the SAME legend placed on the right.

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(ggpubr)
})

# ---------- load & normalize ----------
ratio_raw <- readr::read_tsv("mis.ratio.txt", show_col_types = FALSE)

# normalize column names to lowercase snake_case
names(ratio_raw) <- names(ratio_raw) %>%
  stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
  stringr::str_replace("^_|_$", "") %>%
  tolower()

ratio <- ratio_raw

# helper: pick a required column by regex pattern
col_or_stop <- function(df, pattern) {
  hits <- names(df)[grepl(pattern, names(df), ignore.case = TRUE)]
  if (length(hits) == 0) stop("Missing required column matching: ", pattern)
  hits[1]
}

# harmonize key columns
sample_col <- if ("sample_id" %in% names(ratio)) "sample_id" else if ("sample" %in% names(ratio)) "sample" else names(ratio)[1]
reads_col  <- col_or_stop(ratio, "^n_?reads$")
sars2_col  <- col_or_stop(ratio, "^expected.*without.*sars.*cov.*2$")
bat_col    <- col_or_stop(ratio, "^expected.*if.*without.*clos(e|)t.*bat")

ratio <- ratio %>%
  dplyr::rename(sample_id = !!sample_col,
                n_reads   = !!reads_col) %>%
  dplyr::mutate(
    emr_no_sars2 = .data[[sars2_col]],
    bat_mis      = .data[[bat_col]]
  )

# keep ratios in [0,1]
ratio <- ratio %>%
  dplyr::mutate(
    emr_no_sars2 = dplyr::if_else(emr_no_sars2 >= 0 & emr_no_sars2 <= 1, emr_no_sars2, NA_real_),
    bat_mis      = dplyr::if_else(bat_mis      >= 0 & bat_mis      <= 1, bat_mis,      NA_real_)
  )

# labels & a shared color palette (consistent across plots)
metric_levels <- c("emr_no_sars2", "bat_mis")
metric_labs <- c(
  emr_no_sars2 = "No SARS-CoV-2 in DB",
  bat_mis      = "No SARS-CoV-2 + closest bat (cumulative)"
)
pal <- scales::hue_pal()(length(metric_levels))
names(pal) <- metric_levels

# ---------- p1: distributions (overlaid histograms; legend hidden to share p2's) ----------
dist_long <- ratio %>%
  dplyr::select(all_of(metric_levels)) %>%
  tidyr::pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::mutate(metric = factor(metric, levels = metric_levels))

p1 <- ggplot(dist_long, aes(x = value, fill = metric)) +
  geom_histogram(position = "identity", alpha = 0.45, bins = 50, color = NA) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = pal, labels = metric_labs, breaks = metric_levels) +
  labs(
    x = "Expected missing ratio",
    y = "Samples",
    fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# ---------- p2: n_reads vs missingness (single panel, both metrics) ----------
scatter_long <- ratio %>%
  dplyr::select(n_reads, all_of(metric_levels)) %>%
  tidyr::pivot_longer(all_of(metric_levels),
                      names_to = "metric", values_to = "missing_ratio") %>%
  dplyr::filter(!is.na(n_reads), !is.na(missing_ratio)) %>%
  dplyr::mutate(metric = factor(metric, levels = metric_levels))

# Spearman stats for both metrics
stats <- scatter_long %>%
  dplyr::group_by(metric) %>%
  dplyr::summarise(
    n   = dplyr::n(),
    rho = suppressWarnings(cor(n_reads, missing_ratio, method = "spearman")),
    p   = suppressWarnings(cor.test(n_reads, missing_ratio, method = "spearman", exact = FALSE)$p.value),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = paste0("Spearman \u03C1 = ", sprintf("%.2f", rho),
                   "\nP = ", format.pval(p, digits = 2, eps = 1e-16))
  )

# place the labels in the top-right, staggered to avoid overlap
y_max <- max(scatter_long$n_reads, na.rm = TRUE)
x_max <- max(scatter_long$missing_ratio, na.rm = TRUE)
ann <- stats %>%
  dplyr::arrange(desc(metric)) %>%
  dplyr::mutate(
    x = pmin(x_max - 0.02, 1),
    y = seq(from = y_max * 0.98, to = y_max * 0.86, length.out = dplyr::n())
  )

p2 <- ggplot(scatter_long, aes(x = missing_ratio, y = n_reads, color = metric)) +
  geom_point(alpha = 0.45, size = 1) +
  geom_smooth(se = TRUE) +
  scale_color_manual(values = pal, labels = metric_labs, breaks = metric_levels) +
  scale_y_continuous(labels = label_number(scale_cut = cut_si(""))) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    x = "Expected missing ratio",
    y = "Read depth (n_reads)",
    color = NULL
  ) +
  theme_minimal() +
  geom_text(data = ann, aes(x = x, y = y, label = label, color = metric),
            inherit.aes = FALSE, hjust = 1, vjust = 1, size = 3)

# ---------- arrange with a single shared legend on the right ----------
combined <- ggarrange(
  p2, p1,            # put the legend-carrying plot first
  ncol = 2, widths = c(1.15, 1),
  common.legend = TRUE, legend = "bottom"
)

#p4
library(tidyverse)
library(scales)

# Load and tidy
err <- read.delim("mis.error.frq.txt", check.names = FALSE) %>%
  rename(
    length_bp = length,
    ref       = ref,
    errors    = error,
    freq      = Fre
  ) %>%
  mutate(
    length_bp = as.integer(length_bp),
    errors    = as.integer(errors),
    freq      = as.numeric(freq),
    # Order facets nicely
    length_bp = factor(length_bp, levels = c(50, 100), labels = c("50 bp", "100 bp")),
    ref       = factor(ref, levels = c("SARS-cov-2", "SARS-cov-1", "CloseBat"))
  )

# Optional palette (each panel is a single series; legend not needed)
pal <- c(
  "SARS-cov-2" = "#1f77b4",
  "SARS-cov-1" = "#d62728",
  "CloseBat"   = "#2ca02c"
)

p4 <- ggplot(err, aes(x = errors, y = freq)) +
  geom_col(aes(fill = ref), width = 0.9, show.legend = FALSE) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(
    x = "Number of mismatches (errors)",
    y = "Frequency (reads)"
  ) +
  facet_grid(ref ~ length_bp, scales = "free") +  # free x and y
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor.x = element_blank(),
    strip.text = element_text(face = "bold")
  )


  p4 <- ggplot(err, aes(x = errors, y = freq)) +
  geom_col(aes(fill = ref), width = 0.9, show.legend = FALSE) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(
    x = "Number of mismatches (errors)",
    y = "Frequency (reads)"
  ) +
  facet_grid(ref ~ length_bp, scales = "free") +  # free x and y
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor.x = element_blank(),
    strip.text = element_text(face = "bold")
  )
p4_small <- p4 + theme(
  strip.text.x = element_text(size = 8, face = "bold"),
  strip.text.y = element_text(size = 8, face = "bold")
)

# Use the small-label plot in the final layout
p4_noleg <- p4_small + theme(legend.position = "none")
row_top <- ggarrange(p0, p4_noleg, labels = c("A", "B"), ncol = 2, widths = c(1, 1))
final_fig <- ggarrange(row_top, combined_CD, nrow = 2, heights = c(1, 1.15))
