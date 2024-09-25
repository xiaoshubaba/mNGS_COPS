library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)
data <- read.table("samples.rna.reads.number.txt", header = TRUE)
data_percentage <- data
df_percentage <-data %>%
  rowwise() %>%
  mutate(across(4:6, ~ . / sum(c_across(4:6)) * 100))

data_long <- df_percentage %>%
  pivot_longer(cols = starts_with("host"):clean, names_to = "variable", values_to = "value")

colorsQ = rev(pal_nejm("default")(2))
my_colors <- colorsQ

p = ggboxplot(data_long, x = "cohort", y = "value", color = "cohort",
          add = "jitter", palette = colorsQ, legend = "none") +
  labs(title = "Number of sequences by Cohort",
       y = "Number of sequences (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(method = "t.test", label = "p.signif", hide.ns = TRUE) +
  facet_wrap(~ variable, nrow = 1)

p1 = ggboxplot(data,x="cohort",y="raw",color="cohort",add="jitter", palette= colorsQ,legend = "none") + labs(x="",y="Number of Raw sequences") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +stat_compare_means(method = "t.test", label = "p.signif", hide.ns = TRUE)
f = ggarrage(p,p1,nrow=1)
