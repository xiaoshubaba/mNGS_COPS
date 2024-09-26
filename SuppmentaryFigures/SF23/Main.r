library(ggsci)
library(ggplot2)
colorsQ = c(pal_nejm("default")(3),"grey50")
my_colors <- colorsQ

Main = read.table("main.txt",head=T)
Main$col =  paste(Main$prj,Main$cohort,sep=".")
Main$col = factor(Main$col,levels=c("p1.COVID","p1.Control","p2.COVID","p2.Control","","p3.COVID"))
Main$raw_L = log(Main$raw)
Main$clean_L = log(Main$clean)
Main$gatkReads_L = log(Main$gatk_reads)
Main$gatkReads_rpm = Main$gatk_reads/Main$clean*1e6
Main$asm_cov = Main$asm_cov/100
p1 = ggboxplot(Main, x = "col", y = "raw_L", color = "cohort",add = "jitter",add.params = list(size = 1), outlier.shape = NA, palette = colorsQ, legend = "none") +labs(x = "", y = "Raw Sequences(log)") +  theme(axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8),axis.title.y = element_text(size = 8))
p2 = ggboxplot(Main, x = "col", y = "clean_L", color = "cohort",add = "jitter",add.params = list(size = 1), outlier.shape = NA, palette = colorsQ, legend = "none") +labs(x = "", y = "Clean Sequences(log)") +  theme(axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8),axis.title.y = element_text(size = 8))
p3 = ggboxplot(Main, x = "col", y = "gatkReads_L", color = "cohort",add = "jitter", add.params = list(size = 1), outlier.shape = NA,palette = colorsQ, legend = "none") +labs(x = "", y = "Detected COVID Reads(log)") +  theme(axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8),axis.title.y = element_text(size = 8))
p4 = ggboxplot(Main, x = "col", y = "gatkReads_rpm", color = "cohort",add = "jitter", add.params = list(size = 1), outlier.shape = NA,palette = colorsQ, legend = "none") +labs(x = "", y = "Detected COVID Reads(RPM)") +  theme(axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8),axis.title.y = element_text(size = 8))
p5 = ggboxplot(Main[which(Main$cohort=="COVID"),], x = "col", y = "gatk_cov", color = "cohort",add = "jitter", add.params = list(size = 1), outlier.shape = NA, palette = colorsQ, legend = "none") +labs(x = "", y = "Mapped Reads Genomic Coverage") +  theme(axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8),axis.title.y = element_text(size = 8))
p6 = ggboxplot(Main[which(Main$cohort=="COVID"),], x = "col", y = "asm_cov", color = "cohort",add = "jitter", add.params = list(size = 1), outlier.shape = NA, palette = colorsQ, legend = "none") +labs(x = "", y = "De novo Assembled Genomic Coverage") +  theme(axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8),axis.title.y = element_text(size = 8))
ggarrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2,labels=LETTERS[1:6])
