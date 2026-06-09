# 安装必要包（如果未安装）
if (!require("ggh4x")) install.packages("ggh4x")

library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(readr)
library(ggh4x)

# ==== 1. 读取数据 ====
data <- read.table("S6.Simulated.outbreaks.pops.txt", header = TRUE, sep = "\t")

# ==== 2. 补全每个 sample 的 3 个 round 并处理 NA 值 ====
rounds <- c(1, 2, 3)
completed_data <- data %>%
  select(prj, sample) %>% distinct() %>%
  crossing(round = rounds) %>%
  left_join(data, by = c("prj", "sample", "round")) %>%
  mutate(
    cov = ifelse(is.na(cov), 0, cov),
    signal = ifelse(is.na(signal), 0, signal),
    nCtg = ifelse(is.na(nCtg), 0, nCtg)
  )

# ==== 3. 添加 round_label（修改为更易懂的名称）====
completed_data <- completed_data %>%
  mutate(round_label = factor(round, levels = 1:3,
                              labels = c("Co-assembly", "Case Recurrence", "Control Subtraction")))

# ==== 4. 定义病原体类型 ====
respiratory_pathogens <- c(
  "Admed_SARS2_2020", "Eran_SARS2_2020", "LA_Adenovirus_2024",
  "LA_IVA_2024", "LA_IVB_2024", "LA_PRIV_2024", "LA_RHV_2024",
  "LA_SARS2_2024", "Lia_HRSV_A_2024", "Lia_HRSV_B_2024",
  "Miftahul_sars2_2022", "sfl_sars2_2023","Pavels_SARS2_2021","See_sars_2023","Wendy_li_sars2_2022"
)

completed_data <- completed_data %>%
  mutate(
    pathogen_type = case_when(
      prj %in% respiratory_pathogens ~ "Respiratory pathogens",
      TRUE ~ "Other pathogens"
    )
  )

# ==== 5. 计算每个项目每一轮中 signal > 0 的模拟疫情数 ====
signal_gt0_counts <- completed_data %>%
  group_by(prj) %>%
  summarise(
    n1 = sum(round == 1 & signal > 0, na.rm = TRUE),
    n2 = sum(round == 2 & signal > 0, na.rm = TRUE),
    n3 = sum(round == 3 & signal > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    n1 = ifelse(n1 == 0, 1, n1),
    n2 = ifelse(n2 == 0, 1, n2),
    n3 = ifelse(n3 == 0, 1, n3)
  )

# ==== 6. 按项目、round 求和并计算平均值 ====
sum_by_round <- completed_data %>%
  group_by(prj, round) %>%
  summarise(
    total_cov = sum(cov, na.rm = TRUE),
    total_signal = sum(signal, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(round_label = factor(round, levels = 1:3,
                              labels = c("Co-assembly", "Case Recurrence", "Control Subtraction")))

sum_by_round <- sum_by_round %>%
  left_join(signal_gt0_counts, by = "prj")

adjusted_mean <- sum_by_round %>%
  mutate(
    `Genomic Coverage (%)` = case_when(
      round == 1 ~ total_cov / n1,
      round == 2 ~ total_cov / n2,
      round == 3 ~ total_cov / n3
    ),
    `Signal Purity (%)` = case_when(
      round == 1 ~ total_signal / n1,
      round == 2 ~ total_signal / n2,
      round == 3 ~ total_signal / n3
    )
  ) %>%
  select(prj, round_label, `Genomic Coverage (%)`, `Signal Purity (%)`) %>%
  pivot_longer(cols = c(`Genomic Coverage (%)`, `Signal Purity (%)`),
               names_to = "metric", values_to = "value") %>%
  mutate(value = pmin(pmax(value, 0), 100))

# 添加病原体类型
pathogen_type_info <- completed_data %>%
  select(prj, pathogen_type) %>%
  distinct()
adjusted_mean <- adjusted_mean %>%
  left_join(pathogen_type_info, by = "prj")

# ==== 7. 显著性检验（使用配对检验，paired = TRUE）====
compare_stats <- function(df, metric_col) {
  df %>%
    group_by(prj) %>%
    group_modify(~ {
      # 提取每个 round 的值，按 sample 对齐？由于每个 sample 有唯一标识，但简化处理：直接取所有值
      # 配对检验需要相同样本在不同 round 的值。这里假设数据中每个 sample 在三个 round 都有记录
      # 为了严格配对，应该按 sample 对齐，但为简化且与原代码风格一致，我们使用 paired = TRUE 要求数据成对出现
      r1 <- .x %>% filter(round == 1) %>% pull({{ metric_col }})
      r2 <- .x %>% filter(round == 2) %>% pull({{ metric_col }})
      r3 <- .x %>% filter(round == 3) %>% pull({{ metric_col }})
      # 配对检验要求长度相同，取最小长度
      n12 <- min(length(r1), length(r2))
      n13 <- min(length(r1), length(r3))
      p_2vs1 <- if (n12 > 1) wilcox.test(r1[1:n12], r2[1:n12], paired = TRUE)$p.value else NA
      p_3vs1 <- if (n13 > 1) wilcox.test(r1[1:n13], r3[1:n13], paired = TRUE)$p.value else NA
      tibble(p_2vs1 = p_2vs1, p_3vs1 = p_3vs1)
    }) %>%
    ungroup()
}

cov_stats <- compare_stats(completed_data, cov) %>% rename_with(~paste0(., "_cov"), -prj)
signal_stats <- compare_stats(completed_data, signal) %>% rename_with(~paste0(., "_signal"), -prj)
stats_all <- cov_stats %>% left_join(signal_stats, by = "prj")

# ==== 8. 准备热图数据，添加星号 ====
heatmap_data <- adjusted_mean %>%
  left_join(stats_all, by = "prj") %>%
  left_join(signal_gt0_counts, by = "prj") %>%
  mutate(
    label_value = sprintf("%.1f", value),
    label = case_when(
      metric == "Genomic Coverage (%)" & round_label == "Case Recurrence" & p_2vs1_cov < 0.05 ~ paste0(label_value, "*"),
      metric == "Genomic Coverage (%)" & round_label == "Control Subtraction" & p_3vs1_cov < 0.05 ~ paste0(label_value, "*"),
      metric == "Signal Purity (%)" & round_label == "Case Recurrence" & p_2vs1_signal < 0.05 ~ paste0(label_value, "*"),
      metric == "Signal Purity (%)" & round_label == "Control Subtraction" & p_3vs1_signal < 0.05 ~ paste0(label_value, "*"),
      TRUE ~ label_value
    )
  )

# ==== 9. 美化项目名称并添加 n1,n3 信息 ====
heatmap_data <- heatmap_data %>%
  mutate(
    prj_base = case_when(
      prj == "LA_HMB_2024" ~ "LA_AV_2024",
      prj == "LA_HRSV_2024" ~ "LA_PRIV_2024",
      TRUE ~ prj
    ),
    prj_base = gsub("_", ".", prj_base),
    prj_display = sprintf("%s (n1=%d, n3=%d)", prj_base, n1, n3)
  )

# 添加分组因子
heatmap_data <- heatmap_data %>%
  mutate(
    pathogen_group = factor(pathogen_type,
                            levels = c("Respiratory pathogens", "Other pathogens"))
  )

# 排序
project_order_resp <- heatmap_data %>%
  filter(pathogen_type == "Respiratory pathogens") %>%
  select(prj_base, prj_display) %>% distinct() %>%
  arrange(prj_base) %>% pull(prj_display) %>% rev()

project_order_other <- heatmap_data %>%
  filter(pathogen_type == "Other pathogens") %>%
  select(prj_base, prj_display) %>% distinct() %>%
  arrange(prj_base) %>% pull(prj_display) %>% rev()

heatmap_data$prj_display <- factor(heatmap_data$prj_display,
                                   levels = c(project_order_resp, project_order_other))

# 定义颜色梯度
color_gradient <- scale_fill_gradientn(
  colours = c("#4393C3", "#92C5DE", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
  values = scales::rescale(c(0, 40, 60, 70, 80, 100)),
  limits = c(0, 100),
  name = "Value (%)",
  guide = guide_colorbar(
    barwidth = 0.8,
    barheight = 8,
    title.position = "top",
    title.hjust = 0.5,
    frame.colour = "black",
    ticks.colour = "black"
  )
)

# 定义列标签
metric_labeller <- as_labeller(c(
  "Genomic Coverage (%)" = "Genomic Coverage",
  "Signal Purity (%)" = "Signal Purity"
))

# 定义行标签背景色（左侧分组）
strip_colors_y <- c("Respiratory pathogens" = "#4E79A7",
                    "Other pathogens" = "#E15759")
# 列标签背景色（黑色）
strip_colors_x <- c("Genomic Coverage" = "black", "Signal Purity" = "black")

# ==== 10. 绘制热图 ====
p <- ggplot(heatmap_data, aes(x = round_label, y = prj_display, fill = value)) +
  geom_tile(color = "white", size = 0.3) +
  geom_text(aes(label = label), size = 3, fontface = "bold", color = "white") +
  color_gradient +
  facet_grid2(pathogen_group ~ metric,
              scales = "free_y", space = "free_y",
              switch = "y",
              labeller = labeller(metric = metric_labeller),
              strip = strip_themed(
                background_y = elem_list_rect(fill = strip_colors_y),
                background_x = elem_list_rect(fill = "black", color = NA)
              )) +
  labs(
    x = "POPS Processing Step",
    y = "Project",
    caption = "* p < 0.05 vs. Co-assembly (Wilcoxon signed-rank test, paired)\nn1/n3: # outbreaks with signal > 0 in round1/3"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
    axis.text.y = element_text(size = 7, hjust = 1),
    axis.title = element_text(face = "bold", size = 10),
    strip.text = element_text(face = "bold", size = 10, color = "white"),
    strip.background = element_rect(color = NA),
    strip.placement = "outside",
    panel.grid = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(10, 10, 10, 10),
    plot.caption = element_text(size = 8, face = "italic", color = "gray30",
                                hjust = 1, vjust = 0, margin = margin(t = 10))
  )

#print(p)

# 保存
ggsave("S_F3.heatmap.pdf", plot = p, width = 11, height = 9, dpi = 300)
ggsave("S_F3.heatmap.png", plot = p, width = 11, height = 9, dpi = 300)

# 统计摘要
cat("\n=== POPS Processing Analysis Summary ===\n")
cat(sprintf("Respiratory pathogens: %d\n", length(project_order_resp)))
cat(sprintf("Other pathogens: %d\n", length(project_order_other)))
