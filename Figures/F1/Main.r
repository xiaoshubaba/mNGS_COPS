library(ggsci)
library(ggpubr)

Main = read.table("Main.F1.txt",head=T,sep="\t")
Main$class = NA
Main$class[which(Main$specy== "Severe acute respiratory syndrome-related coronavirus")] = "COVID";Main$class[which(Main$family=="Enterobacteriaceae")] = "Bacteria::Enterobacteriaceae";Main$class[which(Main$family=="Pseudomonadaceae")] = "Bacteria::Pseudomonadaceae";Main$class[which(Main$family=="Burkholderiaceae")] = "Bacteria::Burkholderiaceae";Main$class[which(Main$family=="Neisseriaceae")] = "Bacteria::Neisseriaceae";Main$class[which(Main$family=="Streptococcaceae")] = "Bacteria::Streptococcaceae"
i = 1
while (i<=dim(Main)[1]){
        if( (is.na(Main$superkingdom[i]) == "FALSE") && (is.na(Main$class[i])) ){
        Main$class[i] = "Others"
        }
i = i + 1
}
i = 1
while (i<=dim(Main)[1]){
        if (is.na(Main$class[i])){
        Main$class[i] = "Unannotated"
        }
i = i + 1
}
Main$class = factor(Main$class,levels=c("COVID","Bacteria::Enterobacteriaceae","Bacteria::Pseudomonadaceae","Bacteria::Burkholderiaceae","Bacteria::Neisseriaceae","Bacteria::Streptococcaceae","Others","Unannotated"))
colorsQ = c(pal_nejm("default")(7),"grey50")
f1_left = ggscatter(Main,"gc","log_meanDepth",xlab="GC content",color="class",alpha=I(0.9),legend="none",palette=colorsQ,size="length",title=paste0("contigs:",dim(Main)[1], ", signal ratio:0.157%"))
Main_m6 = Main[which(Main$frequency>=6),]
Main_m6$class = factor(Main_m6$class,levels=c("COVID","Bacteria::Enterobacteriaceae","Bacteria::Pseudomonadaceae","Bacteria::Burkholderiaceae","Bacteria::Neisseriaceae","Bacteria::Streptococcaceae","Others","Unannotated"))
colorsQ_tem = c()
for (j in 1:dim(table(Main_m6$class))){
         if (table(Main_m6$class)[j] > 0){
         colorsQ_tem = c(colorsQ_tem,colorsQ[j])
         }
}
f1_middle = ggscatter(Main_m6,"gc",xlab="GC content","log_meanDepth",color="class",alpha=I(0.9),legend="none",palette=colorsQ_tem,size="length",title=paste0("contigs:",dim(Main_m6)[1],", signal ratio:0.76%"))
Main_m6_m10 = Main_m6[which(Main_m6$CS1Fre<=10),]
Main_m6_m10$class = factor(Main_m6_m10$class,levels=c("COVID","Bacteria::Enterobacteriaceae","Bacteria::Pseudomonadaceae","Bacteria::Burkholderiaceae","Bacteria::Neisseriaceae","Bacteria::Streptococcaceae","Others","Unannotated"))
colorsQ_tem = c()
for (j in 1:dim(table(Main_m6_m10$class))){
         if (table(Main_m6_m10$class)[j] > 0){
         colorsQ_tem = c(colorsQ_tem,colorsQ[j])
         }
}
f1_right = ggscatter(Main_m6_m10,"gc",xlab="GC content","log_meanDepth",color="class",alpha=I(0.9),legend="none",palette=colorsQ_tem,size="length",title=paste0("contigs:",dim(Main_m6_m10)[1], ", signal ratio:91.51%"))
legend = ggscatter(Main,"gc","log_meanDepth",color="class",alpha=I(0.9),legend="right",palette=colorsQ,size="length",title=paste0("contigs:",dim(Main)[1]))
########################
Main = Main_m6
Main$class = factor(Main$class,levels=c("COVID","Bacteria::Enterobacteriaceae","Bacteria::Pseudomonadaceae","Bacteria::Burkholderiaceae","Bacteria::Neisseriaceae","Bacteria::Streptococcaceae","Others","Unannotated"))
colorsQ_tem = c()
for (j in 1:dim(table(Main$class))){
         if (table(Main$class)[j] > 0){
         colorsQ_tem = c(colorsQ_tem,colorsQ[j])
         }
}
p = ggplot(Main, aes(x = meanCov, y = CS1MeanCov, color = class) )+
	geom_point(size = 1.5, alpha = 0.9) +
	scale_color_manual(values = colorsQ_tem) +
	labs(x = "First Week patients Mean Coverage", y = "Stage 1 Controls Mean Coverage") +
	theme(
		axis.title.x = element_text(size = 12),
		axis.title.y = element_text(size = 12),legend.position="none"
	)
#
contigs_to_tag <- c("DNA_k75_10730", "DNA_k75_448", "RNA_k75_9508")
for (contig_name in contigs_to_tag) {
	p <- p + geom_segment(data = Main[Main$contig == contig_name, ],
										aes(x = meanCov, y = CS1MeanCov,
										xend = meanCov + 0.05, yend = CS1MeanCov + 0.05),
										arrow = arrow(length = unit(0.15, "inches")), color = "black") +
						geom_text(data = Main[Main$contig == contig_name, ],
										aes(x = meanCov + 0.05, y = CS1MeanCov + 0.05, label = contig),
										hjust = 0, vjust = 0, size = 3)
}
left = p

p1 = ggplot(Main, aes(x = frequency, y = CS1Fre, color = class) )+
geom_point(size = 1.5, alpha = 0.9) +
scale_color_manual(values = colorsQ_tem) +
labs(x = "First Week patients Frequency", y = "Stage 1 Controls Frequency") +
theme(
		axis.title.x = element_text(size = 12),
		axis.title.y = element_text(size = 12),legend.position="none"
	)
for (contig_name in contigs_to_tag) {
	p1 <- p1 + geom_segment(data = Main[Main$contig == contig_name, ],
										aes(x = frequency, y = CS1Fre,
										xend = frequency + 0.5, yend = CS1Fre + 0.5),
										arrow = arrow(length = unit(0.15, "inches")), color = "black") +
						geom_text(data = Main[Main$contig == contig_name, ],
										aes(x = frequency + 0.5, y = CS1Fre + 0.5, label = contig),
										hjust = 0, vjust = 0, size = 3)
}


p2 = ggplot(Main, aes(x = log(meanDepth), y = log(CS1MeanDep), color = class) )+
	geom_point(size = 1.5, alpha = 0.9) +
	scale_color_manual(values = colorsQ_tem) +
	labs(x = "First Week patients Mean Depth", y = "Stage 1 Controls Mean Depth") +
	theme(
		axis.title.x = element_text(size = 12),
		axis.title.y = element_text(size = 12),legend.position="none"
	)

for (contig_name in contigs_to_tag) {
	p2 <- p2 + geom_segment(data = Main[Main$contig == contig_name, ],
										aes(x =log(meanDepth), y = log(CS1MeanDep),
										xend = log(meanDepth)+0.5 , yend = log(CS1MeanDep) + 0.5),
										arrow = arrow(length = unit(0.15, "inches")), color = "black") +
						geom_text(data = Main[Main$contig == contig_name, ],
										aes(x = log(meanDepth) + 0.5, y = log(CS1MeanDep) + 0.5, label = contig),
										hjust = 0, vjust = 0, size = 3)
}
left  = p
middle = p1
right = p2
ggarragne(f1_left,f1_middle,f1_right,left,middle,right)
