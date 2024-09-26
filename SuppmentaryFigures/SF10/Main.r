library(ggsci)
library(ggpubr)
CalSR <- function(data){
R = sum(data[which(data$Specy==1),]$length)/sum(data$length)
return(round(R,digits=4))
}
colorsQ = c("grey50",pal_nejm("default")(1))

#top
#left
Main <- read.table("first_5P.main.information.txt",head=T,sep="\t")
Main$Specy = factor(Main$Specy)
Main$log_meanDepth = log(Main$meanDepth)
Main$S1_log_meanDepth = log(Main$S1_meanDep)
#right
f1_left = ggscatter(Main,"gc","log_meanDepth",color="Specy",alpha=I(0.9),legend="none",palette=colorsQ,size="length",title=paste0("contigs:",dim(Main)[1]," Signal Ratio:",CalSR(Main)))
Main_m6 = Main[which(Main$frequency>3),]
f1_middle = ggscatter(Main_m6,"gc","log_meanDepth",color="Specy",alpha=I(0.9),legend="none",palette=colorsQ,size="length",title=paste0("contigs:",dim(Main_m6)[1]," Signal Ratio:",CalSR(Main_m6)))
Main_m6_m10 = Main_m6[which(Main_m6$S1_fre<10),]
f1_right = ggscatter(Main_m6_m10,"gc","log_meanDepth",color="Specy",alpha=I(0.9),legend="none",palette=colorsQ[2],size="length",title=paste0("contigs:",dim(Main_m6_m10)[1]," Signal Ratio:",CalSR(Main_m6_m10)))
legend = ggscatter(Main,"gc","log_meanDepth",color="Specy",alpha=I(0.9),legend="right",palette=colorsQ,size="length",title=paste0("contigs:",dim(Main)[1]))
# bottom
p1 <- ggscatter(Main_m6,"log_meanDepth","S1_log_meanDepth",color="Specy",fill="Specy",palette=colorsQ,size="length",legend="none")
p2 <- ggscatter(Main_m6,"frequency","S1_fre",color="Specy",fill="Specy",palette=colorsQ,size="length",legend="none")
p3 <- ggscatter(Main_m6,"meanCov","S1_meanCov",color="Specy",fill="Specy",palette=colorsQ,size="length",legend="none")
