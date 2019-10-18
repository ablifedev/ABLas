
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]
outdir <- args[3]

setwd(outdir)
Data <- read.delim(file=infile,header=T)
Pvalue <-c()
for (i in 1:nrow(Data)) {
	compare <- matrix(c(Data[i,10],Data[i,11],Data[i,15],Data[i,16]),nr=2)
	result <- fisher.test(compare)		#,alternative = "greater"
	# result <- chisq.test(compare,simulate.p.value=F,B=40000)
	Pvalue <- c(Pvalue,result$p.value)
}

Data <- cbind(Data,Pvalue)
write.table(Data,file=outfile,row.names=F,append=F,quote=F,sep='\t')