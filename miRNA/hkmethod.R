

library("limma")
library(edgeR)
setwd("C:/study/lulab_research/毕设/THCA")                 #???ù???Ŀ¼
inputFile="data/filtered_hFTC_miRNA_RPM.txt"                                             
fdrFilter=0.1                                                  
logFCfilter=0.5                                                     
conNum=9                                             
treatNum=9                                                      


outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))


rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
#rownames(rt)=rt[,1]
#rt = rt[,1:41]
exp=rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data) 
data=data[rowMeans(data)>0,]

samples = colnames(rt)
y <- DGEList(counts=rt, samples=samples, group=grade)
y <- calcNormFactors(y, method='TMM')
design <- model.matrix(~grade)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
test <- glmLRT(fit, coef=2)
res <- topTags(test, n=nrow(rt), sort.by='none')
res <- cbind(res$table, baseMean=2^(res$table$logCPM))
# rename columns
mapped_names <- colnames(res)
for(i in 1:ncol(res)){
  if(colnames(res)[i] == 'logFC'){
    mapped_names[i] <- 'log2FoldChange'
  }else if(colnames(res)[i] == 'PValue'){
    mapped_names[i] <- 'pvalue'
  }else if(colnames(res)[i] == 'FDR') {
    mapped_names[i] <- 'padj'
  }else{
    mapped_names[i] <- colnames(res)[i]
  }
}
colnames(res) <- mapped_names
write.table(res, 'TMM.txt', sep='\t', quote=FALSE, row.names=TRUE)




