suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(grid)      # for unit()
})

# ---------- helpers ----------
standardize_names <- function(x) {
  x %>%
    str_replace_all("-", "_") %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("Closet", "Closest") %>%
    str_replace_all("cov_2", "CoV_2") %>%
    str_replace_all("sars_cov_2", "SARS_CoV_2") %>%
    str_replace_all("sars_cov_1", "SARS_CoV_1")
}

find_col <- function(nm, pattern) {
  hits <- grep(pattern, nm, ignore.case = TRUE, value = TRUE)
  if (length(hits) != 1L)
    stop(sprintf("Column match error for pattern: %s\nFound: %s",
                 pattern, paste(hits, collapse=", ")))
  hits
}

to_num <- function(x) {
  readr::parse_number(x, locale = readr::locale(decimal_mark = ".", grouping_mark = ","))
}

read_and_prepare <- function(path, method_label) {
  df_raw <- readr::read_tsv(
    file = path,
    col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA", "NaN"),
    trim_ws = TRUE,
    guess_max = 1e6,
    progress = FALSE
  )
  names(df_raw) <- standardize_names(names(df_raw))
  sample_col     <- find_col(names(df_raw), "^sample$")
  miss_sars2_col <- find_col(names(df_raw), "Expected.*Missing.*Ratio.*without.*SARS.*CoV.*2")
  miss_bat_col   <- find_col(names(df_raw), "Expected.*Missing.*Ratio.*without.*Closest.*Bat")

  df_raw %>%
    transmute(
      sample = .data[[sample_col]],
      `Without SARS-CoV-2`  = to_num(.data[[miss_sars2_col]]),
      `Without Closest Bat` = to_num(.data[[miss_bat_col]])
    ) %>%
    pivot_longer(
      cols = c(`Without SARS-CoV-2`, `Without Closest Bat`),
      names_to = "Scenario",
      values_to = "MissingRatio"
    ) %>%
    mutate(Method = method_label)
}

# POPS (reference-free); collapse repeats per sample by median
read_and_prepare_pops <- function(path) {
  df_raw <- readr::read_tsv(
    file = path,
    col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA", "NaN"),
    trim_ws = TRUE,
    guess_max = 1e6,
    progress = FALSE
  )
  names(df_raw) <- standardize_names(names(df_raw))
  sample_col <- find_col(names(df_raw), "^sample$")
  pops_col   <- find_col(names(df_raw), "Expected.*Missing.*Ratio.*POPS")

  df_raw %>%
    transmute(sample = .data[[sample_col]],
              MissingRatio = to_num(.data[[pops_col]])) %>%
    group_by(sample) %>%
    summarise(MissingRatio = median(MissingRatio, na.rm = TRUE), .groups = "drop") %>%
    mutate(Method = "POPS", Scenario = "Reference-free")
}

# ---------- read files for p1 ----------
blastx_long <- read_and_prepare("mis.protein.ratio.txt", method_label = "BLASTx")
blastn_long <- read_and_prepare("mis.ratio.txt",         method_label = "BLASTn")
pops_long   <- read_and_prepare_pops("mis.ratio.POPs.txt")

df_long_all <- bind_rows(blastn_long, blastx_long, pops_long) %>%
  filter(!is.na(MissingRatio)) %>%
  mutate(
    Method   = factor(Method,   levels = c("BLASTn", "BLASTx", "POPS")),
    Scenario = if_else(Method == "POPS", "Reference-free", Scenario),
    Group = factor(
      case_when(
        Method == "BLASTn" & Scenario == "Without SARS-CoV-2"   ~ "BLASTn (Without SARS-CoV-2)",
        Method == "BLASTn" & Scenario == "Without Closest Bat"  ~ "BLASTn (Without Closest Bat)",
        Method == "BLASTx" & Scenario == "Without SARS-CoV-2"   ~ "BLASTx (Without SARS-CoV-2)",
        Method == "BLASTx" & Scenario == "Without Closest Bat"  ~ "BLASTx (Without Closest Bat)",
        Method == "POPS"                                        ~ "POPS",
        TRUE ~ NA_character_
      ),
      levels = c(
        "BLASTn (Without SARS-CoV-2)",
        "BLASTn (Without Closest Bat)",
        "BLASTx (Without SARS-CoV-2)",
        "BLASTx (Without Closest Bat)",
        "POPS"
      )
    )
  ) %>%
  filter(!is.na(Group))

# ---------- p1: 5-column boxplot ----------
p1 <- ggplot(df_long_all, aes(x = Group, y = MissingRatio, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.9) +
  geom_jitter(width = 0.12, alpha = 0.45, size = 1.1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Expected Missing Ratio: Homology Methods Under Reference Removal vs POPS (Reference-free)",
    x = NULL, y = "Expected Missing Ratio"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 11),  # smaller title
    axis.text.x = element_text(angle = 18, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.title = element_blank()
  )

# ---------- p2: POPS vs SARS-cov-2 error profile ----------
# Read & clean
df_raw <- readr::read_table2("W.error.txt", col_types = cols())
df <- df_raw %>%
  dplyr::rename(
    length = !!names(df_raw)[[1]],
    ref    = !!names(df_raw)[[2]],
    error  = !!names(df_raw)[[3]],
    Fre    = !!names(df_raw)[[4]]
  ) %>%
  dplyr::mutate(
    ref = stringr::str_to_lower(ref) |>
      stringr::str_replace_all("sars[\\-_ ]?co[vv][\\-_ ]?2", "sars-cov-2") |>
      stringr::str_replace_all("pops", "pops"),
    ref = dplyr::case_when(
      ref == "pops" ~ "POPs",
      ref == "sars-cov-2" ~ "SARS-cov-2",
      TRUE ~ ref
    ),
    length = as.integer(length),
    error  = as.integer(error),
    Fre    = as.numeric(Fre)
  ) %>%
  dplyr::filter(ref %in% c("POPs", "SARS-cov-2"), Fre > 0)

# Weighted helpers
w_mean <- function(x, w) sum(w * x) / sum(w)
w_var_rep <- function(x, w) { n <- sum(w); m <- w_mean(x, w); if (n <= 1) return(NA_real_); sum(w * (x - m)^2) / (n - 1) }

welch_weighted <- function(d_sub) {
  g <- d_sub %>%
    dplyr::group_by(ref) %>%
    dplyr::summarise(
      n_eff = sum(Fre),
      mean_error = w_mean(error, Fre),
      var_rep = w_var_rep(error, Fre),
      .groups = "drop"
    )
  if (nrow(g) != 2) return(NULL)
  g1 <- g %>% dplyr::filter(ref == "POPs")
  g2 <- g %>% dplyr::filter(ref == "SARS-cov-2")
  se_diff <- sqrt((g1$var_rep / g1$n_eff) + (g2$var_rep / g2$n_eff))
  diff    <- g1$mean_error - g2$mean_error
  df_welch <- (se_diff^2)^2 /
    ((g1$var_rep / g1$n_eff)^2 / (g1$n_eff - 1) +
       (g2$var_rep / g2$n_eff)^2 / (g2$n_eff - 1))
  tval <- diff / se_diff
  pval <- 2 * pt(-abs(tval), df = df_welch)
  ci_low  <- diff - qt(0.975, df = df_welch) * se_diff
  ci_high <- diff + qt(0.975, df = df_welch) * se_diff
  sd_pooled <- sqrt((g1$var_rep + g2$var_rep) / 2)
  cohens_d  <- diff / sd_pooled
  dplyr::tibble(
    n_pops = g1$n_eff, n_sars = g2$n_eff,
    mean_pops = g1$mean_error, mean_sars = g2$mean_error,
    diff_mean = diff, ci_low = ci_low, ci_high = ci_high,
    t = tval, df = df_welch, p_value = pval,
    cohen_d = cohens_d,
    pct_reduction_pops_vs_sars = (1 - (g1$mean_error / g2$mean_error)) * 100
  )
}

fmt_p <- function(p) {
  p <- as.numeric(p)
  out <- rep(NA_character_, length(p))
  out[is.na(p)] <- "NA"
  out[p == 0] <- "<1e-16"
  out[p > 0 & p < 1e-4] <- "<1e-4"
  out[p >= 1e-4] <- sprintf("%.3g", p[p >= 1e-4])
  out
}

test_results <- df %>%
  dplyr::group_by(length) %>%
  dplyr::group_modify(~{
    wt <- welch_weighted(.x)
    tab <- xtabs(Fre ~ ref + error, data = .x)
    chi <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 5000))
    dplyr::tibble(
      n_levels_error = ncol(tab),
      chisq_stat = unname(chi$statistic),
      chisq_p_value = chi$p.value
    ) %>% dplyr::bind_cols(wt)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(length)

# Build p2 (smaller title; legend will be removed later)
p2 <- df %>%
  ggplot(aes(x = error, y = Fre, color = ref)) +
  geom_line(aes(group = ref)) +
  geom_point(size = 1) +
  facet_wrap(~ length, scales = "free_y") +
  labs(
    x = "Error count",
    y = "Frequency (Fre)",
    color = "Reference",
    title = "Error count distributions by read length"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 11),   # smaller title
    legend.position = "bottom",
    legend.box = "horizontal",
    strip.text = element_text(size = 10)    # smaller facet labels (optional)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

# Add compact top-right stats (smaller annotation text)
labels_tbl <- test_results %>%
  mutate(
    welch_p_str = fmt_p(p_value),
    chisq_p_str = fmt_p(chisq_p_value),
    label = paste(
      paste0("Mean POPs=", sprintf("%.4f", mean_pops),
             " vs SARS-cov-2=", sprintf("%.4f", mean_sars)),
      paste0("\u0394 (POPs\u2212SARS-cov-2)=",
             sprintf("%.4f", diff_mean), " [",
             sprintf("%.4f", ci_low), ", ", sprintf("%.4f", ci_high), "]"),
      paste0("Welch p=", welch_p_str),
      paste0("Chi-square p=", chisq_p_str),
      paste0("Reduction=", sprintf("%.1f", pct_reduction_pops_vs_sars), "%"),
      sep = "\n"
    )
  ) %>%
  select(length, label)

facet_pos <- df %>%
  group_by(length) %>%
  summarise(
    x_pos = max(error, na.rm = TRUE) * 0.98,
    y_pos = max(Fre,   na.rm = TRUE) * 0.98,
    .groups = "drop"
  )

labels_pos <- left_join(labels_tbl, facet_pos, by = "length")

p2 <- p2 +
  geom_label(
    data = labels_pos,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    hjust = 1, vjust = 1,
    label.size = 0.2,
    label.padding = unit(0.06, "lines"),  # tighter box
    size = 2.0,                           # smaller text
    lineheight = 0.9
  )

# Remove legend from p2 for the combined figure
p2_noleg <- p2 + theme(legend.position = "none")

# ---------- p5–p7 (rowC) ----------
# Fallback palettes if not defined in your env
if (!exists("colorsQ"))     colorsQ     <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728")
if (!exists("colorsQ_tem")) colorsQ_tem <- colorsQ

Main <- read.table("mix.V2.express.txt", header = TRUE, sep = "\t")
Main$log_meanDepth <- log(Main$length)

sr <- round(sum(Main[Main$annotation == "Y", ]$length) / sum(Main$length) * 100, 2)
max_gc <- max(Main$gc, na.rm = TRUE)
max_depth <- max(Main$log_meanDepth, na.rm = TRUE)
p5 <- ggscatter(Main, "gc", "log_meanDepth", xlab = "GC content",
                color = "annotation", alpha = 0.9, legend = "none",
                palette = colorsQ, size = "length") +
  annotate("text", x = max_gc, y = max_depth, hjust = 1, vjust = 1, size = 3,
           label = paste0("Combined\nContigs: ", nrow(Main), "\nSignal: ", sr, "%"))

Main_m2 <- Main[which(Main$case.fre >=2), ]
sr <- round(sum(Main_m2[Main_m2$annotation == "Y", ]$length) / sum(Main_m2$length) * 100, 2)
max_gc <- max(Main_m2$gc, na.rm = TRUE)
max_depth <- max(Main_m2$log_meanDepth, na.rm = TRUE)
p6 <- ggscatter(Main_m2, "gc", "log_meanDepth", xlab = "GC content",
                color = "annotation", alpha = 0.9, legend = "none",
                palette = colorsQ_tem, size = "length") +
  annotate("text", x = max_gc, y = max_depth, hjust = 1, vjust = 1, size = 3,
           label = paste0("Case >= 2\nContigs: ", nrow(Main_m2), "\nSignal: ", sr, "%"))

Main_m2_m2 <- Main_m2[which(Main_m2$control.fre <= 1), ]
sr <- round(sum(Main_m2_m2[Main_m2_m2$annotation == "Y", ]$length) / sum(Main_m2_m2$length) * 100, 2)
max_gc <- max(Main_m2_m2$gc, na.rm = TRUE)
max_depth <- max(Main_m2_m2$log_meanDepth, na.rm = TRUE)
p7 <- ggscatter(Main_m2_m2, "gc", "log_meanDepth", xlab = "GC content",
                color = "annotation", alpha = 0.9, legend = "none",
                palette = colorsQ_tem, size = "length") +
  annotate("text", x = max_gc, y = max_depth, hjust = 1, vjust = 1, size = 3,
           label = paste0("Final\nContigs: ", nrow(Main_m2_m2), "\nSignal: ", sr, "%"))

rowC <- ggarrange(p5, p6, p7, ncol = 3, nrow = 1, align = "hv", labels = NULL)

# ---------- Arrange: Row1 (p1|p2) and Row2 (rowC) ----------
row1 <- ggarrange(
  p1, p2_noleg,
  ncol = 2, nrow = 1,
  labels = c("A", "B"),
  font.label = list(size = 14, face = "bold"),
  label.x = 0.01, label.y = 0.99
)

fig_AB_C <- ggarrange(
  row1, rowC,
  ncol = 1, nrow = 2,
  heights = c(1, 1),
  labels = c("", "C"),
  font.label = list(size = 14, face = "bold"),
  label.x = 0.01, label.y = 0.99
)

print(fig_AB_C)

