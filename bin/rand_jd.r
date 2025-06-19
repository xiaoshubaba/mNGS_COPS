#!/usr/bin/Rscript
args=commandArgs(T)
##
# path.file: must have a sample head name
# random number: randomly chose sample number, 2-10
# overlap ratio: jaccard distance, we use 0.8 here
# output.prefix: will out put randomly chose samples and seed number
if (length(args) !=4){
stop("1) path.file 2) random number 3)overlap ratio 4)output.prefix")
}

if (file.exists(args[1])){
dl <- read.table(args[1],sep="\t",header=TRUE)
}else if (!file.exists(args[1])){
stop(paste0(args[1]," path file not eixsts"))
}
dl_CASE = dl[which(dl$cohort == "Case"),]
num = dim(dl_CASE)[1]

ratio = as.numeric(args[3])
num_r = as.numeric(args[2])
i = 1
rand_matrix = c()
rand_samples = c()
seeds = c()
i_seed = i
while (i <= 100){
set.seed(i_seed)
dl_rand = sort(sample(1:num,num_r))
        if (i == 1){
		rand_samples = as.character(dl_CASE[dl_rand,]$sample)
		rand_matrix = dl_rand
		seeds[i] = 1
        }
        else{
		#check
		while(1){
		FLAG = "TRUE";
			j = 1
			matrixNum = dim(rand_matrix)[1]
			if (is.null(matrixNum)){
			matrixNum = 1
			}
			while (j<=matrixNum){
				if (matrixNum == 1){
				mV = rand_matrix
				}
				else{
				mV = rand_matrix[j,]
				}
			reV = c(dl_rand,mV)
			reR = ( length(reV) - length(unique((reV))) ) / num_r
				if (reR >= ratio){
				FLAG = "FALSE"
				print(i_seed)
				print(dl_rand)
				print(mV)
				print(reR)
				print("has mixed ratio")
				break;
				}
			j = j + 1
			}

			if (FLAG == "TRUE"){
			rand_samples = rbind(rand_samples,as.character(dl_CASE[dl_rand,]$sample))
			rand_matrix = rbind(rand_matrix,dl_rand)
			seeds[i] = i_seed
			break;
			}
			else{
			i_seed = i_seed + 1
			set.seed(i_seed)
			dl_rand = sort(sample(1:num,num_r))
			}
		}
        }
i = i + 1
i_seed = i_seed + 1
}
rownames(rand_matrix) = seeds
rownames(rand_samples) = seeds
write.table(rand_matrix,paste0(args[4],".rand.numbers.txt"),sep="\t",quote=FALSE,col.names=FALSE)
write.table(rand_samples,paste0(args[4],".rand.samples.txt"),sep="\t",quote=FALSE,col.names=FALSE)
