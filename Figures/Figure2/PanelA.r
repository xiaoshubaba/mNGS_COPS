library(ggsci)
library(ggpubr)
Main <- read.table("F3d.V4.table.txt",head=T,sep="\t")
Main$log_meanDepth = log(Main$length)
#right
i = 1
Main$class = NA
Main$class[which(Main$Specy== "Severe acute respiratory syndrome coronavirus 2")] = "Sars-Cov-2";
Main$class[which(Main$Specy=="Candida albicans SC5314")] = "Candida albicans";
Main$class[which(Main$Specy=="Influenza A virus")] = "Influenza A virus";
Main$class[which(Main$Specy=="Candida albicans")] = "Candida albicans";
Main$class[which(Main$Specy=="Klebsiella pneumoniae")] = "Klebsiella pneumoniae"
Main$class[which(Main$Specy=="Human respiratory syncytial virus B")] = "HRSV"
i = 1
while (i<=dim(Main)[1]){
        if( is.na(Main$class[i]) ){
        Main$class[i] = "Others"
        }
i = i + 1
}
i = 1
while (i<=dim(Main)[1]){
        if (is.na(Main$Specy[i])){
        Main$class[i] = "Unannotated"
        }
i = i + 1
}
Main$class = factor(Main$class,levels=c("Sars-Cov-2","Candida albicans","Influenza A virus","Klebsiella pneumoniae","HRSV","Others","Unannotated"))
colorsQ = c(pal_nejm("default")(6),"grey50")
sr = round(sum(Main[which(Main$annotation=="Y"),]$length)/sum(Main$length)*100,digits=2)
f1_left = ggscatter(Main,"gc","log_meanDepth",xlab="GC content",color="class",alpha=I(0.9),legend="none",palette=colorsQ,size="length",title=paste0("contigs:",dim(Main)[1], ", signal ratio:",sr,"%"))

Main_m2 = Main[which(Main$case.fre>=2),]
Main_m2$class = factor(Main_m2$class,levels=c("Sars-Cov-2","Candida albicans","Influenza A virus","Klebsiella pneumoniae","HRSV","Others","Unannotated"))
sr = round(sum(Main_m2[which(Main_m2$annotation=="Y"),]$length)/sum(Main_m2$length)*100,digits=2)
colorsQ_tem = c()
for (j in 1:dim(table(Main_m2$class))){
         if (table(Main_m2$class)[j] > 0){
         colorsQ_tem = c(colorsQ_tem,colorsQ[j])
         }
}

f1_middle = ggscatter(Main_m2,"gc",xlab="GC content","log_meanDepth",color="class",alpha=I(0.9),legend="none",palette=colorsQ_tem,size="length",title=paste0("contigs:",dim(Main_m2)[1],", signal ratio:",sr,"%"))
#
Main_m2_m2 = Main_m2[which(Main_m2$control.fre<=2),]
Main_m2_m2$class = factor(Main_m2_m2$class,levels=c("Sars-Cov-2","Candida albicans","Influenza A virus","Klebsiella pneumoniae","HRSV","Others","Unannotated"))
colorsQ_tem = c()
for (j in 1:dim(table(Main_m2_m2$class))){
         if (table(Main_m2_m2$class)[j] > 0){
         colorsQ_tem = c(colorsQ_tem,colorsQ[j])
         }
}
sr = round(sum(Main_m2_m2[which(Main_m2_m2$annotation=="Y"),]$length)/sum(Main_m2_m2$length)*100,digits=2)
f1_right = ggscatter(Main_m2_m2,"gc",xlab="GC content","log_meanDepth",color="class",alpha=I(0.9),legend="none",palette=colorsQ_tem,size="length",title=paste0("contigs:",dim(Main_m2_m2)[1], ", signal ratio:",sr,"%"))

legend = ggscatter(Main,"gc","log_meanDepth",color="class",alpha=I(0.9),legend="right",palette=colorsQ,size="length",title=paste0("contigs:",dim(Main)[1]))
########################
Main = Main_m2
Main$class = factor(Main$class,levels=c("Sars-Cov-2","Candida albicans","Influenza A virus","Klebsiella pneumoniae","HRSV","Others","Unannotated"))
colorsQ_tem = c()
for (j in 1:dim(table(Main$class))){
         if (table(Main$class)[j] > 0){
         colorsQ_tem = c(colorsQ_tem,colorsQ[j])
         }
}
p = ggplot(Main, aes(x = case.mean.cov, y = control.mean.cov, color = class) )+
	geom_point(size = 1.5, alpha = 0.9) +
	scale_color_manual(values = colorsQ_tem) +
	labs(x = "First Week patients Mean Coverage", y = "Controls Mean Coverage") +
	theme(
		axis.title.x = element_text(size = 12),
		axis.title.y = element_text(size = 12),legend.position="none"
	)
#

p1 = ggplot(Main, aes(x = case.fre, y = control.fre, color = class) )+
	geom_point(size = 1.5, alpha = 0.9) +
	scale_color_manual(values = colorsQ_tem) +
	labs(x = "First Week patients Frequency", y = "Controls Frequency") +
	theme(
		axis.title.x = element_text(size = 12),
		axis.title.y = element_text(size = 12),legend.position="none"
	)

p2 = ggplot(Main, aes(x = log(case.mean.depth), y = log(control.mean.depth), color = class) )+
	geom_point(size = 1.5, alpha = 0.9) +
	scale_color_manual(values = colorsQ_tem) +
	labs(x = "First Week patients Mean Depth", y = "Controls Mean Depth") +
	theme(
		axis.title.x = element_text(size = 12),
		axis.title.y = element_text(size = 12),legend.position="none"
	)
left  = p
middle = p1
right = p2
