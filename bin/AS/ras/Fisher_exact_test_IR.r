
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]
outdir <- args[3]

setwd(outdir)
Data <- read.delim(file=infile,header=T)
Data <- Data[,1:15]
Pvalue <-c()
for (i in 1:nrow(Data)) {
	compare <- matrix(c(Data[i,8],Data[i,9],Data[i,12],Data[i,13]),nr=2)
	compare <- round(compare)
	result <- fisher.test(compare)		#,alternative = "greater"
	# result <- chisq.test(compare,simulate.p.value=F,B=40000)
	result$p.value
	Pvalue <- c(Pvalue,result$p.value)
}

Data <- cbind(Data,Pvalue)
write.table(Data,file=outfile,row.names=F,append=F,quote=F,sep='\t')