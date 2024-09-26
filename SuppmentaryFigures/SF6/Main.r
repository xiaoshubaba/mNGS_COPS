library(ggplot2)
library(ggsci)
library(pROC)
colorsQ = c(pal_nejm("default")(7),"grey50")
#panalA
validation = read.table("roc.discovery.S1.txt",head=T)
roc_obj1 <- roc(validation$cohort, validation[,2])
roc_obj2 <- roc(validation$cohort, validation[,3])
roc_obj3 <- roc(validation$cohort, validation[,4])
auc1 <- auc(roc_obj1)
ci1 <- ci.auc(roc_obj1)
auc2 <- auc(roc_obj2)
ci2 <- ci.auc(roc_obj2)
auc3 <- auc(roc_obj3)
ci3 <- ci.auc(roc_obj3)

roc_data1 <- data.frame(
  FPR = 1 - roc_obj1$specificities,
  TPR = roc_obj1$sensitivities,
  Model = "DNA_k75_448"
)
roc_data2 <- data.frame(
  FPR = 1 - roc_obj2$specificities,
  TPR = roc_obj2$sensitivities,
  Model = "RNA_k75_9508"
)

roc_data3 <- data.frame(
  FPR = 1 - roc_obj3$specificities,
  TPR = roc_obj3$sensitivities,
  Model = "DNA_k75_10730"
)

roc_data <- rbind(roc_data1, roc_data2, roc_data3)
custom_colors <- c("DNA_k75_448" = "grey50", "RNA_k75_9508" = "#BC3C29FF", "DNA_k75_10730" = "grey50")
roc_plot <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Model) )+
geom_line(size = 1) + scale_color_manual(values = custom_colors) +
geom_abline(intercept = 0, slope = 1, linetype = "dashed",color = "grey") +
labs(
    title = "ROC Curves of COPs contigs in f1w and ContolS1",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
) +
theme_minimal() +
annotate("text", x = 0.6, y = 0.4, label = paste("Model 1: AUC =", round(auc1, 2), "\n95% CI:", round(ci1[1], 2), "-", round(ci1[3], 2)), color = "grey50") +
annotate("text", x = 0.6, y = 0.3, label = paste("Model 2: AUC =", round(auc2, 2), "\n95% CI:", round(ci2[1], 2), "-", round(ci2[3], 2)), color = "#BC3C29FF") +
annotate("text", x = 0.6, y = 0.2, label = paste("Model 3: AUC =", round(auc3, 2), "\n95% CI:", round(ci3[1], 2), "-", round(ci3[3], 2)), color = "grey50")

#panelB
validation = read.table("roc.validation.S2.txt",head=T)
roc_obj1 <- roc(validation$cohort, validation[,2])
roc_obj2 <- roc(validation$cohort, validation[,3])
roc_obj3 <- roc(validation$cohort, validation[,4])
auc1 <- auc(roc_obj1)
ci1 <- ci.auc(roc_obj1)
auc2 <- auc(roc_obj2)
ci2 <- ci.auc(roc_obj2)
auc3 <- auc(roc_obj3)
ci3 <- ci.auc(roc_obj3)

roc_data1 <- data.frame(
  FPR = 1 - roc_obj1$specificities,
  TPR = roc_obj1$sensitivities,
  Model = "DNA_k75_448"
)
roc_data2 <- data.frame(
  FPR = 1 - roc_obj2$specificities,
  TPR = roc_obj2$sensitivities,
  Model = "RNA_k75_9508"
)

roc_data3 <- data.frame(
  FPR = 1 - roc_obj3$specificities,
  TPR = roc_obj3$sensitivities,
  Model = "DNA_k75_10730"
)

roc_data <- rbind(roc_data1, roc_data2, roc_data3)
custom_colors <- c("DNA_k75_448" = "grey50", "RNA_k75_9508" = "#BC3C29FF", "DNA_k75_10730" = "grey50")
roc_plot <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Model) )+
geom_line(size = 1) + scale_color_manual(values = custom_colors) +
geom_abline(intercept = 0, slope = 1, linetype = "dashed",color = "grey") +
labs(
    title = "ROC Curves with AUC Values",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
) +
theme_minimal() +
annotate("text", x = 0.6, y = 0.4, label = paste("Model 1: AUC =", round(auc1, 2), "\n95% CI:", round(ci1[1], 2), "-", round(ci1[3], 2)), color = "grey50") +
annotate("text", x = 0.6, y = 0.3, label = paste("Model 2: AUC =", round(auc2, 2), "\n95% CI:", round(ci2[1], 2), "-", round(ci2[3], 2)), color = "#BC3C29FF") +
annotate("text", x = 0.6, y = 0.2, label = paste("Model 3: AUC =", round(auc3, 2), "\n95% CI:", round(ci3[1], 2), "-", round(ci3[3], 2)), color = "grey50")

#panel C
cutoff_pool <- as.vector(read.table("roc.validation.cutoff.S2.txt",head=T))[[1]]
tp = c(); fp = c(); fn = c(); tn = c()
i = 1
for (i in 1:(length(cutoff_pool)-1)){
	j = 1
	pool_diagnosis = c()
	while (j <= dim(validation)[1]){
	value = validation[j,]
	diagnosis = as.numeric(value[2:4]) > cutoff_pool[i]
		if (TRUE %in% diagnosis){
		pool_diagnosis[j] = "Vali"
		}
		else{
		pool_diagnosis[j] = "Control"
		}
	j = j + 1
	}
	pool_tem_diagnosis = data.frame(validation,cohort=validation$cohort,pool_diagnosis=pool_diagnosis)
	tp[i] = length(which(pool_tem_diagnosis$cohort=="Vali" & pool_tem_diagnosis$pool_diagnosis=="Vali"))
	tn[i] = length( which(pool_tem_diagnosis$cohort=="Control" & pool_tem_diagnosis$pool_diagnosis=="Control"))
	fp[i] = length( which(pool_tem_diagnosis$cohort=="Control" & pool_tem_diagnosis$pool_diagnosis=="Vali"))
	fn[i] = length( which(pool_tem_diagnosis$cohort=="Vali" & pool_tem_diagnosis$pool_diagnosis=="Control"))
}
Precision=tp/(tp+fp);Recall=tp/(tp+fn);
f1s = 2*((Precision*Recall)/(Precision+Recall))
output = data.frame( cutoff=cutoff_pool[1:149],tp=tp,tn=tn,fp=fp,fn=fn,pre=Precision,rec=Recall,f1s=f1s)
df <- data.frame(cutoff=output$cutoff,pre=output$pre,rec=output$rec,f1s=output$f1s)
dfm <- melt(df, id = "cutoff")
ggplot(dfm, aes(x = cutoff, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("#F8766D", "#00BFC4", "#619CFF")) +  # Custom line colors
  labs(x = "Cops genomic cutoff", y = "Value", color = "Metric") +
  theme_minimal()
