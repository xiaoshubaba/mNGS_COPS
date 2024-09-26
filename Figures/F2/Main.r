library(ggplot2)
library(dplyr)
library(ggsci)
library(tidyr)

# panelA
colorsQ = c(pal_nejm("default")(2),"grey50")
# Read the data
df <- read.table("single.sample.asm.txt", header = TRUE)

# Multiply gatk_cov by 100
df$gatk_cov <- df$gatk_cov * 100
df_high <- df %>% filter(gatk_Reads > 2034)
df_low <- df %>% filter(gatk_Reads <= 2034)

# Calculate Pearson correlation coefficients and p-values
pearson_low_cov <- cor.test(df_low$gatk_Reads, df_low$cov, method = "pearson")
pearson_high_cov <- cor.test(df_high$gatk_Reads, df_high$cov, method = "pearson")

pearson_low_gatk_cov <- cor.test(df_low$gatk_Reads, df_low$gatk_cov, method = "pearson")
pearson_high_gatk_cov <- cor.test(df_high$gatk_Reads, df_high$gatk_cov, method = "pearson")

pearson_low_signal_ratio <- cor.test(df_low$gatk_Reads, df_low$signal_ratio, method = "pearson")
pearson_high_signal_ratio <- cor.test(df_high$gatk_Reads, df_high$signal_ratio, method = "pearson")

p <- ggplot(df, aes(x = gatk_Reads)) +
  geom_point(aes(y = cov, color = "Denovo assemble genomic coverage"), alpha = 0.6) +
  geom_point(aes(y = gatk_cov, color = "Mapped genomic coverage"), alpha = 0.6) +
  geom_point(aes(y = signal_ratio, color = "Signal Ratio"), alpha = 0.6) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "Detected Pathogen Reads (log10 scale)", y = "Value") +
  scale_color_manual(values = c("Denovo assemble genomic coverage" = colorsQ[1], "Mapped genomic coverage" = colorsQ[2], "Signal Ratio" = colorsQ[3])) +
  geom_smooth(data = df_low, aes(y = cov, color = "Denovo assemble genomic coverage"), method = "gam", se = TRUE, linetype = "solid") +
  geom_smooth(data = df_high, aes(y = cov, color = "Denovo assemble genomic coverage"), method = "gam", se = TRUE, linetype = "solid") +
  geom_smooth(data = df_low, aes(y = gatk_cov, color = "Mapped genomic coverage"), method = "gam", se = TRUE, linetype = "solid") +
  geom_smooth(data = df_high, aes(y = gatk_cov, color = "Mapped genomic coverage"), method = "gam", se = TRUE, linetype = "solid") +
  geom_smooth(data = df_low, aes(y = signal_ratio, color = "Signal Ratio"), method = "gam", se = TRUE, linetype = "solid") +
  geom_smooth(data = df_high, aes(y = signal_ratio, color = "Signal Ratio"), method = "gam", se = TRUE, linetype = "solid") +
  geom_vline(xintercept = 2034, linetype = "dashed", color = "grey") +
  annotate("text", x = 2034, y = -1, label = "x = 2034", hjust = -0.1, color = "grey50") +
  annotate("text", x = 0, y = 100, label = paste0("Low (Denovo coverage): ", round(pearson_low_cov$estimate, 2), ", p = ", round(pearson_low_cov$p.value, 2)), hjust = -0.1, color = colorsQ[1], size = 4) +
  annotate("text", x = 10^4, y = 0, label = paste0("High (Denovo coverage): ", round(pearson_high_cov$estimate, 2), ", p = ", round(pearson_high_cov$p.value, 2)), hjust = -0.1, color = colorsQ[1], size = 4) +
  annotate("text", x = 0, y = 90, label = paste0("Low (Mapped coverage): ", round(pearson_low_gatk_cov$estimate, 2), ", p = ", round(pearson_low_gatk_cov$p.value, 2)), hjust = -0.1, color = colorsQ[2], size = 4) +
  annotate("text", x = 10^4, y = 10, label = paste0("High (Mapped coverage): ", round(pearson_high_gatk_cov$estimate, 2), ", p = ", round(pearson_high_gatk_cov$p.value, 2)), hjust = -0.1, color = colorsQ[2], size = 4) +
  annotate("text", x = 0, y = 80, label = paste0("Low (Signal Ratio): ", round(pearson_low_signal_ratio$estimate, 2), ", p = ", round(pearson_low_signal_ratio$p.value, 2)), hjust = -0.1, color = colorsQ[3], size = 4) +
  annotate("text", x = 10^4, y = 20, label = paste0("High (Signal Ratio): ", round(pearson_high_signal_ratio$estimate, 2), ", p = ", round(pearson_high_signal_ratio$p.value, 2)), hjust = -0.1, color = colorsQ[3], size = 4) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.box = "horizontal")

#panelB
getLong <- function(data){
names(data) <- trimws(names(data))
data$cov <- as.numeric(as.character(data$cov))
data$signal_ratio <- as.numeric(as.character(data$signal_ratio))
data$sampleNum <- as.factor(data$sampleNum)
data$completeNess <- as.numeric(data$completeNess)
long_data <- data %>%
  pivot_longer(cols = c(cov, signal_ratio), names_to = "variable", values_to = "value", values_transform = list(value = as.numeric))
# means
summary_data <- long_data %>%
  group_by(sampleNum, variable) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.25, na.rm = TRUE),
    upper = quantile(value, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )
summary_data$sampleNum = factor(summary_data$sampleNum,levels=1:20) # last number
return(summary_data)
}

# Function to read data and assign a sample number
getF <- function(fileName, number) {
  R <- read.table(fileName, header = TRUE)

  if (dim(R)[1] < 200){
  fill_R = c()
  len = 200 - dim(R)[1]
        for (i in 1:len){
        fill_tem = data.frame(sample=paste0("R.",number,".","blank"),fragment_num=0,mapped_fragment_num=0,mapped_ref_bases=0,BreakPoint=0,cov=0,snp=0,indel=0,signal_ratio=NA)
        fill_R = rbind(fill_R,fill_tem)
        }
  R = rbind(R,fill_R)
  }
  R$sampleNum <- number
  R$completeNess <- R$cov>95 & R$mapped_fragment_num==1
  return(R)
}

bootstrap_stat <- function(data, indices) {
  d <- data[indices]
  mean(d, na.rm = TRUE)
}


# Reading the data from two files and binding them
R1 <- getF("R1.assess.txt", "1")
R2 <- getF("R2.assess.txt", "2")
R3 <- getF("R3.assess.txt", "3")
R4 <- getF("R4.assess.txt", "4")
R5 <- getF("R5.assess.txt", "5")
R6 <- getF("R6.assess.txt", "6")
R7 <- getF("R7.assess.txt", "7")
R8 <- getF("R8.assess.txt", "8")
R9 <- getF("R9.assess.txt", "9")
R10 <- getF("R10.assess.txt", "10")
R11 <- getF("R11.assess.txt", "11")
R12 <- getF("R12.assess.txt", "12")
R13 <- getF("R13.assess.txt", "13")
R14 <- getF("R14.assess.txt", "14")
R15 <- getF("R15.assess.txt", "15")
R20 <- getF("R20.assess.txt", "20")
#
RL1 <- getF("RL1.assess.txt", "1")
RL2 <- getF("RL2.assess.txt", "2")
RL3 <- getF("RL3.assess.txt", "3")
RL4 <- getF("RL4.assess.txt", "4")
RL5 <- getF("RL5.assess.txt", "5")
RL6 <- getF("RL6.assess.txt", "6")
RL7 <- getF("RL7.assess.txt", "7")
RL8 <- getF("RL8.assess.txt", "8")
RL9 <- getF("RL9.assess.txt", "9")
RL10 <- getF("RL10.assess.txt", "10")
RL11 <- getF("RL11.assess.txt", "11")
RL12 <- getF("RL12.assess.txt", "12")
RL13 <- getF("RL13.assess.txt", "13")
RL14 <- getF("RL14.assess.txt", "14")
RL15 <- getF("RL15.assess.txt", "15")
RL20 <- getF("RL20.assess.txt", "20")
data <- rbind(R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15,R20)
dataL <- rbind(RL1,RL2,RL3,RL4,RL5,RL6,RL7,RL8,RL9,RL10,RL11,RL12,RL13,RL14,RL15,RL20)
data$project = "All samples"
dataL$project = "low COVID samples"
data <- rbind(data,dataL)
data$sampleNum = factor(data$sampleNum,levels=c(1:15,20))
RL1 <- getF("RL1.assess.txt", "1")
RL2 <- getF("RL2.assess.txt", "2")
RL3 <- getF("RL3.assess.txt", "3")
RL4 <- getF("RL4.assess.txt", "4")
RL5 <- getF("RL5.assess.txt", "5")
RL6 <- getF("RL6.assess.txt", "6")
RL7 <- getF("RL7.assess.txt", "7")
RL8 <- getF("RL8.assess.txt", "8")
RL9 <- getF("RL9.assess.txt", "9")
RL10 <- getF("RL10.assess.txt", "10")
RL11 <- getF("RL11.assess.txt", "11")
RL12 <- getF("RL12.assess.txt", "12")
RL13 <- getF("RL13.assess.txt", "13")
RL14 <- getF("RL14.assess.txt", "14")
RL15 <- getF("RL15.assess.txt", "15")
RL20 <- getF("RL20.assess.txt", "20")
data <- rbind(R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15,R20)
dataL <- rbind(RL1,RL2,RL3,RL4,RL5,RL6,RL7,RL8,RL9,RL10,RL11,RL12,RL13,RL14,RL1
5,RL20)
data$project = "All samples"
dataL$project = "low COVID samples"
data <- rbind(data,dataL)
data$sampleNum = factor(data$sampleNum,levels=c(1:15,20))

data_long <-data %>%  pivot_longer(cols = c(cov, signal_ratio), names_to = "Metric", values_to = "Value")
p <- ggplot(data_long, aes(x = as.factor(sampleNum), y = Value, color = Metric)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.4, size = 1) +
  facet_grid(rows = vars(project), cols = vars(Metric), scales = "free_x", switch = "y",
             labeller = labeller(Metric = c("cov" = "De novo assembled genomic coverage",
                                            "signal_ratio" = "Signal ratio"))) +
  scale_color_manual(values = c("cov" = "#BC3C29FF", "signal" = "grey50")) +
  labs(x = "Number of pooled samples", y = "Value", color = "Metric") +
  theme(axis.title = element_text(size = 10), strip.placement = "outside") +
  stat_summary(
    aes(group = Metric),
    fun = mean,
    geom = "line",
    position = position_dodge(width = 0.75),
    size = 1,
    linetype = "solid"
  ) +
  guides(color = FALSE)  # 1) Hide the legend
#panelC
Main = read.table("Table.5.single.F3.V2.txt",head=T,sep="\t")
Main$project = factor(Main$project,levels=c("SARS-Cov-2 study 1","SARS-Cov-2 study 2","HIV","Zika","Ebola","LASV","Dengue","CHIKV"))
colorsQ = c(pal_nejm("default")(2),"grey50")
Main$gatkCov = Main$gatkCov * 100
p <- ggplot(Main, aes(x = gatkReads)) +
  geom_point(aes(y = asmCov, color = "Denovo assemble coverage"), alpha = 0.6,size=1) +
  geom_point(aes(y = gatkCov, color = "Mapped coverage"), alpha = 0.6) +
  geom_point(aes(y = Signal, color = "Signal Ratio"), alpha = 0.6) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) + scale_color_manual(values = c("Denovo assemble coverage" = colorsQ[1], "Mapped coverage" = colorsQ[2], "Signal Ratio" = colorsQ[3])) +
  facet_wrap(~ project, ncol = 1,nrow=8,scales = "free_x",strip.position="left") +  # Split into 3 plots based on 'project' variable
  labs(
    y = "Coverage/Signal",
    x = "Number of detected pathogen reads (log10 scale)",
    color = "Metrics"
  ) +
  theme(axis.title = element_text(size = 10), # Adjust axis label sizes
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),legend.position="none") + geom_smooth(data = Main, aes(y = asmCov, color = "Denovo assemble coverage"), method = "gam", se = TRUE, linetype = "solid") +
  geom_smooth(data = Main, aes(y = gatkCov, color = "Mapped coverage"), method = "gam", se = TRUE, linetype = "solid") +
  geom_smooth(data = Main, aes(y = Signal, color = "Signal Ratio"), method = "gam", se = TRUE, linetype = "solid") +   ylim(0, 100)
