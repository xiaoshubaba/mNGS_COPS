library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)
#SF2-A
Main = read.table("Single.gatk.reads.txt",head=T,sep="\t")
colorsQ = pal_nejm("default")(2)
Main$gatkCov = Main$gatkCov * 100
sf2A <- ggplot(Main, aes(x = gatkReads)) +
  geom_point(aes(y = asmCov, color = "Denovo assemble coverage"), alpha = 0.6,size=1) +
  geom_point(aes(y = Signal, color = "Signal Ratio"), alpha = 0.6) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) + scale_color_manual(values = c("Denovo assemble coverage" = colorsQ[1], "Signal Ratio" = colorsQ[2])) +
  labs(
    y = "Coverage/Signal",
    x = "Number of detected pathogen reads (log10 scale)",
    color = "Metrics"
  ) +
  theme(axis.title = element_text(size = 10), # Adjust axis label sizes
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),legend.position="none") + geom_smooth(data = Main, aes(y = asmCov, color = "Denovo assemble coverage"), method = "gam", se = TRUE, linetype = "solid") +
  geom_smooth(data = Main, aes(y = Signal, color = "Signal Ratio"), method = "gam", se = TRUE, linetype = "solid") +   ylim(0, 100)
#SF2-B

library(scales)

# Load data
df <- read_tsv("Table.5.mix.F3.txt") %>%
  mutate(Mixnumber = as.integer(Mixnumber))

# ---------------------------------------
# 1. Define Success: cov > 0 = successful
# ---------------------------------------
df <- df %>%
  mutate(success = ifelse(cov > 0, 1, 0))

# ---------------------------------------
# 2. Multi-Objective Scoring
# ---------------------------------------
score_df <- df %>%
  group_by(Mixnumber) %>%
  summarise(
    mean_cov = mean(cov),
    mean_signal = mean(signal),
    success_rate = mean(success),
    .groups = "drop"
  ) %>%
  mutate(
    norm_cov = rescale(mean_cov),
    norm_signal = rescale(mean_signal),
    norm_success = rescale(success_rate),
    # Weighted multi-objective score
    score = ((norm_cov) +
            (norm_signal) +
            (norm_success))/3
  )
score_df <- score_df %>%
  mutate(
    score = (norm_cov * norm_signal * norm_success)^(1/3)
  )
# ---------------------------------------
# 3. Logistic Model Fit for Coverage
# ---------------------------------------
logistic_fit <- nls(mean_cov ~ SSlogis(Mixnumber, Asym, xmid, scal), data = score_df)
score_df$logistic_cov <- predict(logistic_fit, newdata = score_df)
inflection_mix <- coef(logistic_fit)["xmid"]

# ---------------------------------------
# 4. Export Final Table
# ---------------------------------------
write_csv(score_df, "Mixnumber_Optimization_Summary.csv")

# ---------------------------------------
# 5. Print Formula and Inflection Point
# ---------------------------------------
 plot_df <- score_df %>%
  dplyr::select(Mixnumber, norm_cov, norm_signal, norm_success, score) %>%
  pivot_longer(
    cols = -Mixnumber,
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(Metric,
                    levels = c("norm_cov", "norm_signal", "norm_success", "score"),
                    labels = c("Normalized Coverage", "Normalized Signal", "Normalized Success Rate", "Final Score"))
  )
 nejm_colors = pal_nejm("default")(4)
 p <- ggplot(plot_df, aes(x = Mixnumber, y = Value, color = Metric)) +
  geom_line(alpha = 0.6, size = 1.2) +
  geom_point(alpha = 0.6, size = 2) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "red", linewidth = 0.6) +
  scale_color_manual(values = nejm_colors, name = "Metric") +
  scale_x_continuous(
    breaks = c(1, 3, 5, 7, 9),  # Only show these ticks
    labels = c("1", "3", "5", "7", "9")
  ) +
  labs(
    x = "Mixnumber",
    y = "Normalized Value / Score"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 11)
  )
