# --- Figure 4: Combined panels A, B, C, D ---
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
    .default = x
  )
}

# ================================
# Panel A: Control fragments
# ================================
dfA_raw <- read_tsv(
  "F2A.num.txt",
  col_types = cols(
    prj = col_character(),
    pathogen = col_character(),
    control_index = col_double(),
    `num_non-pathogen_segments` = col_double(),
    ratio = col_double()
  )
)

dfA <- dfA_raw %>%
  mutate(prj = fix_project_names(prj)) %>%
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
  labs(
    title = "A. Control fragments",
    x = NULL,
    y = "Non-pathogen fragments per control"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "plain"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 15)  # Adjust margins
  )

# ================================
# Panel B: Samples to 90% saturation
# ================================
sum_files <- list.files("txt", pattern = "summary\\.txt$", full.names = TRUE)

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
    plot.margin = margin(5, 5, 5, 15)  # Adjust margins
  )

# ================================
# Panel C: Saturation curve with 95% CI
# ================================
sim_files <- list.files("txt", pattern = "^p.*sim.*\\.txt$", full.names = TRUE)

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
      ) %>%
      select(sim_id, frac_samples, ratio)
  )
}

sim_list <- map(sim_files, ~ tryCatch(read_one(.x), error = function(e) NULL))
sim_list <- sim_list[!vapply(sim_list, is.null, logical(1))]

grid <- seq(0, 1, by = 0.01)

interp_one <- function(df) {
  df <- df %>% arrange(frac_samples)
  yi <- approx(x = df$frac_samples, y = df$ratio, xout = grid,
               method = "linear", rule = 2, ties = "ordered")$y
  tibble(frac = grid, ratio = yi)
}

interp_all <- map(sim_list, interp_one)
interp_df  <- bind_rows(interp_all, .id = "sim_idx")

sum_df <- interp_df %>%
  group_by(frac) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    lower      = quantile(ratio, 0.025, na.rm = TRUE),
    upper      = quantile(ratio, 0.975, na.rm = TRUE),
    .groups    = "drop"
  )

find_x_at_y <- function(df, y_target = 0.90) {
  df <- df %>% arrange(frac)
  approx(x = df$ratio, y = df$frac, xout = y_target,
         method = "linear", rule = 2, ties = "ordered")$y
}

x90_vec <- interp_df %>%
  group_by(sim_idx) %>%
  summarise(x90 = find_x_at_y(cur_data(), 0.90), .groups = "drop") %>%
  pull(x90)

x90_med <- median(x90_vec, na.rm = TRUE)
x90_lo  <- quantile(x90_vec, 0.025, na.rm = TRUE)
x90_hi  <- quantile(x90_vec, 0.975, na.rm = TRUE)

pC <- ggplot(sum_df, aes(x = frac, y = mean_ratio)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = alpha(col_main, 0.20), color = NA) +
  geom_line(color = col_main, linewidth = 1.1) +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = x90_med, linetype = "dashed", color = "grey40") +
  annotate("text",
           x = 0.3, y = 0.2,
           label = paste0("~", percent(x90_med, accuracy = 1),
                          " of controls to reach 90%"),
           size = 3.5, vjust = 0, color = "grey25") +
  scale_x_continuous(labels = percent_format(accuracy = 10), limits = c(0, 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 10), limits = c(0, 1)) +
  labs(
    title = "C. Saturation curve",
    x = "Proportion of controls used (%)",
    y = "Cumulative non-pathogen removal (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title        = element_text(face = "plain"),
    panel.grid.minor  = element_blank(),
    plot.margin = margin(5, 15, 5, 15)  # Adjust margins
  )

# ================================
# Panel D: Multiple control datasets improve signal purity
# ================================
df_raw <- read_tsv(
  "Figure4D.Multiple.source.V2.txt",
  col_types = cols(
    outbreak        = col_character(),
    prj             = col_character(),
    Combine_purity  = col_double(),
    Single_purity   = col_double()
  )
)

df_long <- df_raw %>%
  pivot_longer(
    cols = c(Single_purity, Combine_purity),
    names_to = "control_set",
    values_to = "purity"
  ) %>%
  mutate(
    control_set = recode(control_set,
                         "Single_purity"  = "Internal controls",
                         "Combine_purity" = "Multiple control datasets")
  )

counts_map <- tibble::tribble(
  ~project,             ~n_internal, ~n_combined,
  "Admed_SARS2_2020",         31,         871,
  "Grundy_HIV_2023",           9,        2455,
  "Tony_Ebola_2014",          48,        2494,
  "Saha_CHIKV_2019",          91,        2750,
  "Judith_LASV_2023",        119,        2565
)

present_projects <- sort(unique(df_long$prj))
counts_map <- counts_map %>% filter(project %in% present_projects)

label_map <- counts_map %>%
  transmute(
    prj = project,
    `Internal controls`         = paste0("n=", n_internal),
    `Multiple control datasets` = paste0("n=", n_combined)
  ) %>%
  pivot_longer(-prj, names_to = "control_set", values_to = "x_lab")

df_plot <- df_long %>%
  inner_join(label_map, by = c("prj", "control_set")) %>%
  mutate(
    control_set = factor(control_set, levels = c("Internal controls", "Multiple control datasets")),
    x_lab = factor(x_lab)
  ) %>%
  arrange(prj, control_set) %>%
  mutate(x_lab = factor(x_lab, levels = unique(x_lab)))

pval_tbl <- df_long %>%
  group_by(prj) %>%
  summarise(
    p_val = tryCatch(
      wilcox.test(purity ~ control_set)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_lab = ifelse(is.na(p_val), "p = NA", paste0("p = ", format.pval(p_val, digits = 3)))
  )

ypos_tbl <- df_plot %>%
  group_by(prj) %>%
  summarise(y_pos = min(110, max(purity, na.rm = TRUE) * 1.08), .groups = "drop")

anno_df <- counts_map %>%
  transmute(
    prj = project,
    group1 = paste0("n=", n_internal),
    group2 = paste0("n=", n_combined)
  ) %>%
  left_join(ypos_tbl, by = "prj") %>%
  left_join(pval_tbl, by = "prj") %>%
  transmute(
    prj,
    group1,
    group2,
    y.position = y_pos,
    label = p_lab
  )

pD <- ggplot(df_plot, aes(x = x_lab, y = purity, fill = control_set, color = control_set)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35, linewidth = 0.4) +
  geom_jitter(width = 0.12, height = 0, color = "grey60", alpha = 0.35, size = 0.9) +
  facet_wrap(~ prj, nrow = 1, scales = "free_x", strip.position = "bottom") +
  ggpubr::stat_pvalue_manual(
    anno_df,
    label        = "label",
    xmin         = "group1",
    xmax         = "group2",
    y.position   = "y.position",
    tip.length   = 0.01,
    bracket.size = 0.4,
    inherit.aes  = FALSE
  ) +
  scale_fill_manual(values = c("Internal controls" = col_int,
                               "Multiple control datasets" = col_comb),
                    name = NULL) +
  scale_color_manual(values = c("Internal controls" = col_int,
                                "Multiple control datasets" = col_comb),
                     guide = "none") +
  coord_cartesian(ylim = c(0, 110)) +
  labs(
    title = "D. Multiple control datasets improve signal purity",
    x = NULL,
    y = "Signal purity (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title        = element_text(face = "plain"),
    panel.grid.minor  = element_blank(),
    legend.position   = "right",
    strip.placement   = "outside",
    strip.background  = element_rect(fill = NA, color = "grey80", linewidth = 0.4),
    strip.text        = element_text(face = "bold"),
    panel.spacing.x   = unit(16, "pt"),
    panel.border      = element_rect(color = "grey85", fill = NA, linewidth = 0.4),
    plot.margin = margin(5, 15, 5, 15)  # Adjust margins
  )

# ================================
# Combine all panels with improved layout
# ================================
figure_4 <- (pA + pB) / pC / pD + 
  plot_layout(heights = c(1, 1, 1.2)) +
  plot_annotation(theme = theme(plot.margin = margin(10, 10, 10, 10)))

# Print the complete figure
print(figure_4)

# Save with appropriate dimensions
ggsave("figure4.png", figure_4, width = 12, height = 14, dpi = 300)
