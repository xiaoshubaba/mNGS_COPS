library(ggplot2)
library(dplyr)
library(tidyr)

# 读取数据
data <- read.table("mNGS.rna.lib.com.V1.txt", header = TRUE)
data_percentage <- data
df_percentage <-data %>%
  rowwise() %>%
  mutate(across(2:7, ~ . / sum(c_across(3:8)) * 100))

# 将数据转换为长格式
data_long <- df_percentage %>%
  pivot_longer(cols = starts_with("others"):covid, names_to = "variable", values_to = "value")
data_long <- data_long %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate(covid_unclassified = covid + unclassified) %>%
  pivot_longer(cols = c(others, otherVirus, fungi, unclassified, bacteria, covid, covid_unclassified),
               names_to = "variable",
               values_to = "value")


# 定义颜色
colorsQ = c(pal_nejm("default")(5),"grey50","grey30")
my_colors <- colorsQ

# 修正变量名称的拼写并确保顺序
summary_data$variable <- factor(summary_data$variable, levels = c("unclassified","covid","bacteria", "fungi", "otherVirus","others","covid_unclassified"))

p = ggboxplot(data_long, x = "cohort", y = "value", color = "variable",
          add = "jitter", palette = colorsQ, legend = "none") +
  labs(title = "Distribution of Variables by Cohort",
       x = "Cohort",
       y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(method = "t.test", label = "p.signif", hide.ns = TRUE) +
  facet_wrap(~ variable, nrow = 1)
