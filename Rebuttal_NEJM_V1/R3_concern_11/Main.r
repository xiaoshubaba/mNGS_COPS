library(ggsci)
library(ggpubr)
Main <- read.table("PRJ.11.cyto.2patients.txt",head=T,sep="\t")
Main$log_meanDepth = log(Main$length)
#right
colorsQ = c("grey50",pal_nejm("default")(1))
sr = round(sum(Main[which(Main$annotation=="Y"),]$length)/sum(Main$length)*100,digits=2)
f1_left = ggscatter(Main,"gc","log_meanDepth",xlab="GC content",color="annotation",alpha=I(0.9),legend="none",palette=colorsQ,size="length",title=paste0("contigs:",dim(Main)[1], ", signal ratio:",sr,"%"))

Main_m2 = Main[which(Main$case.fre>=2),]
sr = round(sum(Main_m2[which(Main_m2$annotation=="Y"),]$length)/sum(Main_m2$length)*100,digits=2)

colorsQ_tem = pal_nejm("default")(1)
f1_middle = ggscatter(Main_m2,"gc",xlab="GC content","log_meanDepth",color="annotation",alpha=I(0.9),legend="none",palette=colorsQ_tem,size="length",title=paste0("contigs:",dim(Main_m2)[1],", signal ratio:",sr,"%"))
#
Main_m2_m2 = Main_m2[which(Main_m2$control.fre<=1),]
colorsQ_tem = pal_nejm("default")(1)
sr = round(sum(Main_m2_m2[which(Main_m2_m2$annotation=="Y"),]$length)/sum(Main_m2_m2$length)*100,digits=2)
f1_right = ggscatter(Main_m2_m2,"gc",xlab="GC content","log_meanDepth",color="annotation",alpha=I(0.9),legend="none",palette=colorsQ_tem,size="length",title=paste0("contigs:",dim(Main_m2_m2)[1], ", signal ratio:",sr,"%"))

########################
Main = Main_m2
colorsQ_tem = pal_nejm("default")(1)
p = ggplot(Main, aes(x = case.mean.cov, y = control.mean.cov, color = annotation) )+
	geom_point(size = 1.5, alpha = 0.9) +
	scale_color_manual(values = colorsQ_tem) +
	labs(x = "patients Mean Coverage", y = "Controls Mean Coverage") +
	theme(
		axis.title.x = element_text(size = 12),
		axis.title.y = element_text(size = 12),legend.position="none"
	)
#

p1 = ggplot(Main, aes(x = case.fre, y = control.fre, color = annotation) )+
	geom_point(size = 1.5, alpha = 0.9) +
	scale_color_manual(values = colorsQ_tem) +
	labs(x = "patients Frequency", y = "Controls Frequency") +
	theme(
		axis.title.x = element_text(size = 12),
		axis.title.y = element_text(size = 12),legend.position="none"
	)

p2 = ggplot(Main, aes(x = log(case.mean.depth), y = log(control.mean.depth), color = annotation) )+
	geom_point(size = 1.5, alpha = 0.9) +
	scale_color_manual(values = colorsQ_tem) +
	labs(x = "patients Mean Depth", y = "Controls Mean Depth") +
	theme(
		axis.title.x = element_text(size = 12),
		axis.title.y = element_text(size = 12),legend.position="none"
	)
left  = p
middle = p1
right = p2
ggarrange(f1_left,f1_middle,f1_right,left,middle,right,labels=c("A","B","C","D","E","F"))
