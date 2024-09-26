library(ggplot2)
library(dplyr)
 library(patchwork)
df = read.table("table.txt",head=T)
df$sampleType = "COVID"
	 spearman_raw <- cor.test(log10(df$gatk_Reads), log10(df$raw), method = "spearman")
	 spearman_clean <- cor.test(log10(df$gatk_Reads), log10(df$clean), method = "spearman")
	 p1 <- ggplot(df, aes(x = log10(gatk_Reads), y = log10(raw))) +
	   geom_point(color = "#BC3C29FF") +
	   geom_smooth(method = "loess", color = "#BC3C29FF") +
	   labs(
	        x = "COVID reads (log10)",
	        y = "Raw reads (log10)") + annotate("text", x = 2, y = max(log10(df$raw), na.rm = TRUE) - 0.5, label= paste0("Spearman: ", round(spearman_raw$estimate, 2), "\nP-value: ", format(spearman_raw$p.value, digits = 3)), color = "#BC3C29FF")
	 p2 <- ggplot(df, aes(x = log10(gatk_Reads), y = log10(clean))) +
	   geom_point(color = "grey50") +
	   geom_smooth(method = "loess", color = "grey50") +
	   labs(
			x = "COVID reads (log10)",
			y = "Clean reads (log10)") +
annotate("text", x = 2, y = max(log10(df$raw), na.rm = TRUE) - 0.5, label= paste0("Spearman: ", round(spearman_raw$estimate, 2), "\nP-value: ", format(spearman_clean$p.value, digits = 3)), color = "grey50")

df$L_gatk_reads = log10(df$gatk_Reads)
 df$sampleType=factor(df$sampleType,levels=c("BALF","sputum","nasopharyngeal","Lung"))
#sample  fragment_num    mapped_fragment_num     mapped_ref_bases        BreakPoint      cov     snp     indel   signal_ratio    gatk_length     gatk_cov        gatk_Depth      gatk_Segments   gatk_Reads      sampleType      raw     clean
p = ggboxplot(df,"sampleType","L_gatk_reads",add="jitter",xlab="")
p.1 = ggboxplot(df,"sampleType","cov",add="jitter",xlab="")
p.2 = ggboxplot(df,"sampleType","gatk_cov",add="jitter",xlab="")
p.3 = ggboxplot(df,"sampleType","signal_ratio",add="jitter",xlab="")
p12 = ggarrange(p1,p2,nrow=1,ncol=2)
p0123 = ggarrange(p,p.1,p.2,p.3,nrow=1,ncol=4)
pf = ggarrange(p0123,p12,nrow=2,ncol=1)
