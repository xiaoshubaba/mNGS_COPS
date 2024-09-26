library(ggsci)
library(ggpubr)
CalSR <- function(data){
R = sum(data[which(data$Specy==1),]$length)/sum(data$length)
return(round(R,digits=4))
}
C = pal_nejm("default")(2)

#top
#left
Main <- read.table("hiv.V2.feature.txt",head=T,sep="\t")
Main$Specy = factor(Main$Specy)
Main$log_meanDepth = log(Main$meanDepth)
Main$S1_log_meanDepth = log(Main$S1_meanDep)
Main$ifCops = paste(Main$frequency>=4,Main$S1_fre<10,Main$Specy)

#right
colorsQ = c("grey50","grey50",C[1],"grey50",C[2],C[1])

f1_left = ggscatter(Main,"gc","log_meanDepth",color="ifCops",alpha=I(0.9),legend="none",palette=colorsQ,size="length",title=paste0("contigs:",dim(Main)[1]," Signal Ratio:",CalSR(Main)))
Main_m6 = Main[which(Main$frequency>=4),]
colorsQ = c("grey50",C[2],C[1])
f1_middle = ggscatter(Main_m6,"gc","log_meanDepth",color="ifCops",alpha=I(0.9),legend="none",palette=colorsQ,size="length",title=paste0("contigs:",dim(Main_m6)[1]," Signal Ratio:",CalSR(Main_m6)))
Main_m6_m10 = Main_m6[which(Main_m6$S1_fre<10),]
f1_right = ggscatter(Main_m6_m10,"gc","log_meanDepth",color="ifCops",alpha=I(0.9),legend="none",palette=c(C[2],C[1]),size="length",title=paste0("contigs:",dim(Main_m6_m10)[1]," Signal Ratio:",CalSR(Main_m6_m10)))
# bottom
p1 <- ggscatter(Main_m6,"log_meanDepth","S1_log_meanDepth",color="ifCops",palette=colorsQ,size="length",legend="none")
p2 <- ggscatter(Main_m6,"frequency","S1_fre",color="ifCops",palette=colorsQ,size="length",legend="none")
p3 <- ggscatter(Main_m6,"meanCov","S1_meanCov",color="ifCops",palette=colorsQ,size="length",legend="none")
ggarrange(f1_left,f1_middle,f1_right,p1,p2,p3)
