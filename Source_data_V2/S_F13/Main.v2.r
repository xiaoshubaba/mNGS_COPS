# --- Figure 4: Combined panels A, B, C ---
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(scales)
  library(stringr)
  library(purrr)
  library(patchwork)
  library(grid)
})

# NEJM palette
nejm_cols <- pal_nejm("default")(4)
col_int   <- nejm_cols[1]
col_comb  <- nejm_cols[2]
col_main  <- nejm_cols[1]
fill_alpha <- 0.35
jit_alpha  <- 0.35
jit_col    <- "grey60"

# ================================
# Panel A: Control fragments
# ================================
dfA_raw <- read_tsv(
  "S1.PRJ.Controls.status",
  col_types = cols(
    prj = col_character(),
    pathogen = col_character(),
    control_index = col_double(),
    `num_non-pathogen_segments` = col_double(),
    ratio = col_double()
  )
)

dfA <- dfA_raw %>%
  rename(num_non_pathogen_segments = `num_non-pathogen_segments`)

orderA <- dfA %>%
  group_by(prj) %>%
  summarise(med = median(num_non_pathogen_segments, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(prj)

dfA <- dfA %>%
  mutate(prj = factor(prj, levels = rev(orderA)))

pA <- ggplot(dfA, aes(x = prj, y = num_non_pathogen_segments)) +
  geom_boxplot(
    outlier.shape = NA,
    fill = scales::alpha(col_main, fill_alpha),
    color = col_main,
    linewidth = 0.4
  ) +
  geom_jitter(
    width = 0.15, height = 0,
    color = jit_col, alpha = jit_alpha, size = 0.9
  ) +
  coord_flip() +
  # Apply log10 scale to the numerical values (original y-axis)
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  labs(
    title = "A. Control fragments",
    x = NULL,
    y = "Non-pathogen fragments per control (log10 scale)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "plain"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 15)
  )

# ================================
# Panel B: Samples to 90% saturation
# ================================
# ---------- Helper: project name standardization ----------
fix_project_names <- function(x) {
  dplyr::recode(
    x,
    "p1"     = "LA_SARS2_2024",
    "p2"  = "Admed_SARS2_2020",
    "p4"   = "Grundy_HIV_2023",
    "p5"   = "Sardi_Zika_2016",
    "p6"   = "Tony_Ebola_2014",
    "p7"  = "Judith_LASV_2023",
    "p8" = "Zhang_Dengue_2021",
    "p9"   = "Saha_CHIKV_2019",
    "p17"       = "LA_IVA_2024",
    "p18"       = "LA_IVB_2024",
    "p19"       = "LA_RHV_2024",
    "p20"       = "LA_Adenovirus_2024",
    "p21"      = "LA_PRIV_2024",
	"p51"      = "Sarkar_deng_2023",
	"p52" = "Natasha_HIV_2022",
	"p53" = "Ramachadran_HIV_2022",
	"p62" = "Miftahul_sars2_2022",
	"p64" = "Lia_HRSV_A_2024",
	"p65" = "Lia_HRSV_B_2024",
	"p68" = "Eran_SARS2_2020",
	"p70" = "sfl_sars2_2023",
	"p74" = "Yadav_Deng_2024",
	"p46" = "Gabor_WNV_2025",
	"p76" = "Pavels_SARS2_2021",
	"p77" = "Wendy_li_sars2_2022",
	"p78" = "See_sars_2023",
    .default = x
  )
}

sum_files <- list.files("01_random_txt", pattern = "summary\\.txt$", full.names = TRUE)
read_summary <- function(fp) {
  tb <- read_tsv(
    fp,
    col_types = cols(
      sim = col_double(),
      threshold = col_double(),
      totalContigs = col_double(),
      orderAtThreshold = col_double()
    )
  )
  if (!"prj" %in% names(tb)) {
    stem <- str_remove(basename(fp), "\\.summary\\.txt$")
    tb$prj <- stem
  }
  tb
}

dfB_all <- map_df(sum_files, read_summary)

dfB <- dfB_all %>%
  mutate(prj = fix_project_names(prj)) %>%
  filter(threshold >= 0.9 - 1e-9, threshold <= 0.9 + 1e-9)

orderB <- dfB %>%
  group_by(prj) %>%
  summarise(med = median(orderAtThreshold, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(prj)

dfB <- dfB %>%
  mutate(prj = factor(prj, levels = rev(orderB)))

pB <- ggplot(dfB, aes(x = prj, y = orderAtThreshold)) +
  geom_boxplot(
    outlier.shape = NA,
    fill = scales::alpha(col_main, fill_alpha),
    color = col_main,
    linewidth = 0.4
  ) +
  geom_jitter(
    width = 0.15, height = 0,
    color = jit_col, alpha = jit_alpha, size = 0.9
  ) +
  coord_flip() +
  labs(
    title = "B. Saturation samples",
    x = NULL,
    y = "Samples to 90% cumulative"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "plain"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 15)
  )

# ================================
# Panel C: Saturation curve with 95% CI
# ================================
sim_files <- list.files("01_random_txt", pattern = "^p.*sim.*\\.txt$", full.names = TRUE)

read_one <- function(fp) {
  suppressWarnings(
    read_tsv(
      fp,
      col_types = cols(
        order   = col_double(),
        fre     = col_double(),
        frePool = col_double(),
        ratio   = col_double()
      )
    ) %>%
      arrange(order) %>%
      mutate(
        sim_id       = basename(fp),
        frac_samples = order / max(order, na.rm = TRUE),
        ratio        = pmin(pmax(ratio, 0), 1)
      )
  )
}

# Read all simulation files
sim_list <- map(sim_files, ~ tryCatch(read_one(.x), error = function(e) NULL))
sim_list <- sim_list[!vapply(sim_list, is.null, logical(1))]

# Check if any files were read
if (length(sim_list) == 0) {
  stop("No simulation files found in 'txt' directory with pattern '^p.*sim.*\\.txt$'")
}

# Define grid for interpolation
grid <- seq(0, 1, by = 0.01)

# Function to interpolate one simulation
interp_one <- function(df) {
  df <- df %>% arrange(frac_samples)
  yi <- approx(x = df$frac_samples, y = df$ratio, xout = grid,
               method = "linear", rule = 2, ties = "ordered")$y
  tibble(frac = grid, ratio = yi)
}

# Interpolate all simulations
interp_all <- map(sim_list, interp_one)
interp_df  <- bind_rows(interp_all, .id = "sim_idx")

# Calculate summary statistics
sum_df <- interp_df %>%
  group_by(frac) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    lower      = quantile(ratio, 0.025, na.rm = TRUE),
    upper      = quantile(ratio, 0.975, na.rm = TRUE),
    .groups    = "drop"
  )

# Function to find x at target y
find_x_at_y <- function(df, y_target = 0.90) {
  df <- df %>% arrange(frac)
  approx(x = df$ratio, y = df$frac, xout = y_target,
         method = "linear", rule = 2, ties = "ordered")$y
}

# Calculate x90 for each simulation
x90_results <- interp_df %>%
  group_by(sim_idx) %>%
  summarise(x90 = find_x_at_y(cur_data(), 0.90), .groups = "drop")

x90_vec <- x90_results$x90
x90_med <- median(x90_vec, na.rm = TRUE)
x90_lo  <- quantile(x90_vec, 0.025, na.rm = TRUE)
x90_hi  <- quantile(x90_vec, 0.975, na.rm = TRUE)

# Create the plot
pC <- ggplot(sum_df, aes(x = frac, y = mean_ratio)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = alpha(col_main, 0.20), color = NA) +
  geom_line(color = col_main, linewidth = 1.1) +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = x90_med, linetype = "dashed", color = "grey40") +
  # Add vertical lines for confidence interval
  geom_vline(xintercept = x90_lo, linetype = "dotted", color = "grey60", alpha = 0.6) +
  geom_vline(xintercept = x90_hi, linetype = "dotted", color = "grey60", alpha = 0.6) +
  # Add combined annotation
  annotate("text",
           x = 0.25, y = 0.15,
           label = paste0(round(x90_med * 100, 1), "% of controls needed for 90% removal\n",
                         "95% CI (", round(x90_lo * 100, 1), "-",
                         round(x90_hi * 100, 1), "%)"),
           size = 3.2, vjust = 0, hjust = 0, color = "grey25",
           lineheight = 0.9) +
  # Scales
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 10),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 10),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  # Labels
  labs(
    title = "C. Saturation curve",
    x = "Proportion of controls used",
    y = "Cumulative non-pathogen removal",
    caption = paste0("Based on ", length(sim_list), " simulations\n",
                    "Dashed line: median, dotted lines: 95% CI")
  ) +
  # Theme
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "plain"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 15, 5, 15),
    plot.caption = element_text(size = 8, color = "gray40", hjust = 0)
  )

# Print summary statistics
cat("\n=== Saturation Curve Analysis ===\n")
cat(sprintf("Number of simulations: %d\n", length(sim_list)))
cat(sprintf("Median proportion to reach 90%%: %.1f%%\n", x90_med * 100))
cat(sprintf("95%% Confidence Interval: [%.1f%%, %.1f%%]\n", x90_lo * 100, x90_hi * 100))
cat(sprintf("Range: %.1f%% to %.1f%%\n", min(x90_vec, na.rm = TRUE) * 100,
            max(x90_vec, na.rm = TRUE) * 100))

# ================================
# Combine panels A, B, C with improved layout
# ================================
figure_3 <- (pA + pB) / pC +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(theme = theme(plot.margin = margin(10, 10, 10, 10)))

# Print the complete figure (3 panels only)
print(figure_3)
