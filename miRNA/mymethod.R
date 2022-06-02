

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")

###normalization#####
library("limma")
library(edgeR)
setwd("C:/study/lulab_research/毕设/THCA_1_10/")                 
                                  
fdrFilter=0.05                                            
logFCfilter=1
EVNum=41  #EV                                          
csNum=37 #cf                                        


outTab=data.frame()
grade=c(rep(1,EVNum),rep(2,csNum))

rt=read.table("data/filtered_all_miRNA_RPM.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
exp=rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data) 
data=data[rowMeans(data)>0,]

# rawcount = read.table("data/filtered_all_miRNA_count_matrix.txt",sep="\t",header=T,check.names=F)
# rownames(rawcount) <- rawcount[,1]
# rawcount <- rawcount[,-1]
# rawcount = as.matrix(rawcount)
# rawdata=matrix(as.numeric(rawcount),nrow=nrow(rawcount),dimnames=list(rownames(rawcount),colnames(rawcount)))
# 


tmp1 <- data[2,1:EVNum]
tmp2 <- data[2,(EVNum+1):ncol(data)]
for(i in 1:nrow(data)){
  data[i,1:EVNum] <- data[i,1:EVNum]/tmp1
  data[i,(EVNum+1):ncol(data)] <-  data[i,(EVNum+1):ncol(data)]/tmp2
}
write.table(data,file="result/norm.txt",sep="\t",row.names=T,quote=F)

remove(tmp1,tmp2)
for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  #geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(data[i,1:EVNum])
  treatGeneMeans=mean(data[i,(EVNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:EVNum])
  treatMed=median(data[i,(EVNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)
write.table(outTab,file="result/all.xls",sep="\t",row.names=F,quote=F)

outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
#outDiff=outTab[( (as.numeric(as.vector(outTab$logFC))>0.5 | as.numeric(as.vector(outTab$logFC)) < -0.5) & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="result/diff.xls",sep="\t",row.names=F,quote=F)


heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="result/diffExp.txt",sep="\t",col.names=F,quote=F)


pdf(file="result/vol_hERV.pdf",height=5,width=5)
outTab[outTab==Inf]<-NA
outTab[is.na(outTab)] <- 0
xMax=max(abs(as.numeric(as.vector(outTab$logFC))))
yMax=max(-log10(outTab$fdr))+1
plot(as.numeric(as.vector(outTab$logFC)), -log10(outTab$fdr), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano hERV", ylim=c(0,yMax),xlim=c(-10,10),yaxs="i",pch=20, cex=0.8)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="red",cex=0.8)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))<(logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="green",cex=0.8)
abline(v=0,lty=2,lwd=3)
dev.off()


library(pheatmap)
hmExp2=data[as.vector(outDiff[,1]),]
for(i in 1:nrow(hmExp2)){
  hmExp2[i,]=scale(hmExp2[i,])
}


# hmExp=data[as.vector(outDiff[,1]),]
# hmExp=log2(hmExp+0.01)
Type=c(rep("EV",EVNum),rep("cf",csNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
#z-score ??ͼ
pdf(file="result/heatmap_hERV_zscore2.pdf",height=6,width=10)
out <- pheatmap(hmExp2, 
         annotation=Type, 
         #scale = "row",
         color = colorRampPalette(c("blue","white","red"))(500),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = T,
         fontsize = 12,
         fontsize_row=1,
         fontsize_col=1)
dev.off()


# diff <- read.table("data/FTCvsFTA_up_down.txt",sep="\t",header=T,check.names=F)
# diff <- rownames(diff)
# intersect(diff,common_miRNA)

common_miRNA <- row.names(hmExp2[out$tree_row[["order"]],])
rownames(outDiff) <- outDiff[,1]
outDiff <- outDiff[,-1]
diff <- outDiff[common_miRNA,]


enrich_EV <- rownames(diff[as.numeric(outDiff[common_miRNA,3]) < 0,])
enrich_cf <- rownames(diff[as.numeric(outDiff[common_miRNA,3]) > 0,])
write.table(HeatMap_final,file="result/HeatMap_final.txt",sep="\t",col.names=T,quote=F)



pf_FTA <- read.table("data/pf_FTA.txt",sep="\t",check.names=F)
pf_FTC <- read.table("data/pf_FTC.txt",sep="\t",check.names=F)
rownames(pf_FTA) <- pf_FTA[,1]
pf_FTA <- pf_FTA[,-1]
rownames(pf_FTC) <- pf_FTC[,1]
pf_FTC <- pf_FTC[,-1]
colnames(pf_FTA) <- colnames(up)
colnames(pf_FTC) <- colnames(up)
pf_FTA_EV <- pf_FTA[pf_FTA[,1]>0,]
pf_FTA_cf <- pf_FTA[pf_FTA[,1]<0,]
pf_FTC_EV <- pf_FTC[pf_FTC[,1]>0,]
pf_FTC_cf <- pf_FTC[pf_FTC[,1]<0,]
intersect(rownames(pf_FTA_EV),rownames(up))
intersect(rownames(pf_FTA_EV),rownames(down))
intersect(rownames(pf_FTA_cf),rownames(up))
intersect(rownames(pf_FTA_cf),rownames(down))

intersect(rownames(pf_FTC_EV),rownames(up))
intersect(rownames(pf_FTC_EV),rownames(down))
intersect(rownames(pf_FTC_cf),rownames(up))
intersect(rownames(pf_FTC_cf),rownames(down))

EV_RPM <- read.table("data/filtered_miRNA_RPM.txt",sep="\t",check.names=F)
cf_RPM <- read.table("data/cf.txt",sep="\t",check.names=F)

pf_non <- setdiff(c(rownames(up),rownames(down)),c(rownames(pf_FTA) ,rownames(pf_FTC) ))
uni_EV_up <- intersect(intersect(rownames(pf_FTA_EV),rownames(up)),intersect(rownames(pf_FTC_EV),rownames(up)))
uni_cf_down <- intersect(intersect(rownames(pf_FTA_cf),rownames(down)),intersect(rownames(pf_FTC_cf),rownames(down)))
FTA_EV_up <- setdiff(intersect(rownames(pf_FTA_EV),rownames(up)),uni_EV_up)
FTC_EV_up <- setdiff(intersect(rownames(pf_FTC_EV),rownames(up)),uni_EV_up)
FTA_EV_down <- intersect(rownames(pf_FTA_EV),rownames(down))
FTA_cf_up <- intersect(rownames(pf_FTA_cf),rownames(up))
FTC_cf_down <- setdiff(intersect(rownames(pf_FTC_cf),rownames(down)),uni_cf_down)
HeatMap_final <- rbind(cbind(EV_RPM[pf_non,],cf_RPM[pf_non,]),
                       cbind(EV_RPM[FTA_EV_up,],cf_RPM[FTA_EV_up ,]),
                       cbind(EV_RPM[FTA_EV_down,],cf_RPM[FTA_EV_down,]),
                       cbind(EV_RPM[FTA_cf_up,],cf_RPM[FTA_cf_up,]),
                       cbind(EV_RPM[FTC_EV_up,],cf_RPM[FTC_EV_up,]),
                       cbind(EV_RPM[FTC_cf_down,],cf_RPM[FTC_cf_down ,]),
                       cbind(EV_RPM[uni_EV_up,],cf_RPM[uni_EV_up,]),
                       cbind(EV_RPM[uni_cf_down,],cf_RPM[uni_cf_down,]))
HeatMap_final <- as.matrix(HeatMap_final)
HeatMap_final <- log2(HeatMap_final + 1)
for(i in 1:nrow(HeatMap_final)){
  HeatMap_final[i,]=scale(HeatMap_final[i,])
}
write.table(HeatMap_final,file="result/HeatMap_final.txt",sep="\t",col.names=T,quote=F)



up <- read.table("data/FTCvsFTA_up_down.txt",sep="\t",check.names=F)
pf <- read.table("data/pf2.txt",sep="\t",check.names=F)
rownames(pf) <- pf[,1]
pf <- pf[,-1]
pf_EV <- pf[pf[,1]>0,]
pf_cf <- pf[pf[,1]<0,]

pf_EV_up <- intersect(rownames(pf_EV),rownames(up))
intersect(rownames(pf_EV),rownames(down))
intersect(rownames(pf_cf),rownames(up))
pf_cf_down <- intersect(rownames(pf_cf),rownames(down))
pf_non <- setdiff(c(rownames(up),rownames(down)),c(pf_EV_up,pf_cf_down))
HeatMap_pf <- rbind(cbind(EV_RPM[pf_non,],cf_RPM[pf_non,]),cbind(EV_RPM[pf_EV_up,],cf_RPM[pf_EV_up,]),cbind(EV_RPM[pf_cf_down,],cf_RPM[pf_cf_down,]))
HeatMap_pf <- as.matrix(HeatMap_pf)
for(i in 1:nrow(HeatMap_pf)){
  HeatMap_pf[i,]=scale(HeatMap_pf[i,])
}
write.table(HeatMap_pf,file="result/HeatMap_pf.txt",sep="\t",col.names=T,quote=F)



write.table(rawdata[pf_EV_up,],file="result/pf_EV_up.txt",sep="\t",col.names=T,quote=F)
rt <- rawdata[pf_EV_up,]
for(i in 1:nrow(rt)){
  rt[i,]=scale(rt[i,])
}
write.table(rt,file="result/pf_EV_up_scaled.txt",sep="\t",col.names=T,quote=F)

down <- up[up[,1] < 0,]
up <- up[up[,1] > 0,]

up_EV <- intersect(rownames(up),enrich_EV)
down_EV <-intersect(rownames(down),enrich_EV)
up_cf <- intersect(rownames(up),enrich_cf)
down_cf <- intersect(rownames(down),enrich_cf)

up_non <- setdiff(rownames(up),c(up_EV,up_cf))
down_non <- setdiff(rownames(down),c(down_cf))

HeatMap <- rbind(data[up_non,],data[down_non,],data[up_EV,],data[up_cf,],data[down_cf,])
for(i in 1:nrow(HeatMap)){
  HeatMap[i,]=scale(HeatMap[i,])
}
write.table(HeatMap,file="result/HeatMap.txt",sep="\t",col.names=T,quote=F)

pdf(file="result/vol.pdf",height=5,width=5)
# outTab[outTab==Inf]<-NA
# outTab[is.na(outTab)] <- 0
volcanol <- up <- read.table("data/FTCvsFTA.edgeR.glmlrt.txt",sep="\t",check.names=F)

xMax=max(abs(as.numeric(as.vector(volcanol$log2FoldChange))))
yMax=max(-log10(volcanol$padj))+1
plot(as.numeric(as.vector(volcanol$log2FoldChange)), -log10(volcanol$padj), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano of FTC/FTA differential analysis", ylim=c(0,yMax),xlim=c(-3,3),yaxs="i",pch=20, cex=0.8)
diffSub=subset(volcanol, padj<fdrFilter & as.numeric(as.vector(log2FoldChange))>logFCfilter)
points(as.numeric(as.vector(diffSub$log2FoldChange)), -log10(diffSub$padj), pch=20, col="red",cex=0.8)
diffSub=subset(volcanol, padj<fdrFilter & as.numeric(as.vector(log2FoldChange))< -logFCfilter)
points(as.numeric(as.vector(diffSub$log2FoldChange)), -log10(diffSub$padj), pch=20, col="green",cex=0.8)

abline(v=0,lty=2,lwd=2)
abline(h=1,lty=2,lwd=1)

dev.off()








###### special #####
hmExp=data[as.vector(outTab[,1]),]
for(i in 1:nrow(hmExp)){
  hmExp[i,]=scale(hmExp[i,])
}

# hmExp_spec_FTC <- cbind(hmExp_FTC[FTC_spec,],hmExp[FTC_spec,]) # FTC special
# hmExp_spec_FTA <- cbind(hmExp[FTA_spec,],hmExp_FTA[FTA_spec,]) # FTA special

final_hm <- rbind(rbind(hmExp_spec_FTC,hmExp_spec_FTA),cbind(hmExp_rank_FTC[common_miRNA,],hmExp_rank_FTA[common_miRNA,]))
write.table(final_hm,file="final_heatmap.txt",sep="\t",col.names=T,quote=F)




#######longRNA#############
cf_diff <-  read.table("data/cf_FTCvsFTA_up_down_0.05.txt",sep="\t",check.names=F)
EV_diff <-  read.table("data/EV_FTCvsFTA_up_down_0.05.txt",sep="\t",check.names=F)
intersect(rownames(cf_diff),rownames(EV_diff))

pf_cf_diff <-  read.table("data/pf_cf_FTCvsFTA_up_down_0.05.txt",sep="\t",check.names=F)
row.names(pf_cf_diff) <- pf_cf_diff[,1]
pf_cf_diff <- pf_cf_diff[,-1]
pf_EV_diff <-  read.table("data/pF_EV_FTCvsFTA_up_down_0.05.txt",sep="\t",check.names=F)
row.names(pf_EV_diff) <- pf_EV_diff[,1]
pf_EV_diff <- pf_EV_diff[,-1]
intersect(rownames(pf_cf_diff),rownames(pf_EV_diff))

common_EV_diff<- EV_diff[intersect(rownames(pf_EV_diff),rownames(EV_diff)),]
common_cf_diff<- cf_diff[intersect(rownames(pf_cf_diff),rownames(cf_diff)),]

FTA_enrich <- read.table("data/FTA_EVvscf_up_down_0.05.txt",sep="\t",check.names=F)
FTC_enrich <- read.table("data/FTC_EVvscf_up_down_0.05.txt",sep="\t",check.names=F)
FTA_enriched <- intersect(rownames(FTA_enrich),rownames(common_EV_diff))
FTC_enriched <- intersect(rownames(FTC_enrich),rownames(common_EV_diff))


noenriched <- setdiff(rownames(common_EV_diff),union(rownames(FTA_enrich),rownames(FTC_enrich)))
FTA_EV <- FTA_enriched[FTA_enrich[FTA_enriched,1] > 0]
FTA_cf <- FTA_enriched[FTA_enrich[FTA_enriched,1] < 0]
FTC_EV <- FTC_enriched[FTC_enrich[FTC_enriched,1] > 0]
FTC_cf <- FTC_enriched[FTC_enrich[FTC_enriched,1] < 0]

FTA_spec_EVenriched <- setdiff(FTA_EV,universal_EVenriched)
FTA_spec_cfenriched <- setdiff(FTA_cf,universal_cfenriched)
FTC_spec_EVenriched <- setdiff(FTC_EV,universal_EVenriched)
FTC_spec_cfenriched <- setdiff(FTC_cf,universal_cfenriched) 

universal_EVenriched <- intersect(FTA_EV,FTC_EV)
universal_cfenriched <- intersect(FTA_cf,FTC_cf)


TPM_EV <- read.table("data/EV_filtered_RNA_TPM.txt",sep="\t")
TPM_cf <- read.table("data/cf_filtered_RNA_TPM.txt",sep="\t")

heatmap_draw <- rbind(
cbind(TPM_EV[noenriched,21:41],TPM_cf[noenriched,9:18],TPM_EV[noenriched,1:20],TPM_cf[noenriched,1:8]),
cbind(TPM_EV[FTA_spec_EVenriched,21:41],TPM_cf[FTA_spec_EVenriched,9:18],TPM_EV[FTA_spec_EVenriched,1:20],TPM_cf[FTA_spec_EVenriched,1:8]),
cbind(TPM_EV[FTA_spec_cfenriched,21:41],TPM_cf[FTA_spec_cfenriched,9:18],TPM_EV[FTA_spec_cfenriched,1:20],TPM_cf[FTA_spec_cfenriched,1:8]),
cbind(TPM_EV[FTC_spec_EVenriched,21:41],TPM_cf[FTC_spec_EVenriched,9:18],TPM_EV[FTC_spec_EVenriched,1:20],TPM_cf[FTC_spec_EVenriched,1:8]),
cbind(TPM_EV[FTC_spec_cfenriched,21:41],TPM_cf[FTC_spec_cfenriched,9:18],TPM_EV[FTC_spec_cfenriched,1:20],TPM_cf[FTC_spec_cfenriched,1:8]),
cbind(TPM_EV[universal_EVenriched,21:41],TPM_cf[universal_EVenriched,9:18],TPM_EV[universal_EVenriched,1:20],TPM_cf[universal_EVenriched,1:8]),
cbind(TPM_EV[universal_cfenriched,21:41],TPM_cf[universal_cfenriched,9:18],TPM_EV[universal_cfenriched,1:20],TPM_cf[universal_cfenriched,1:8]))



heatmap_draw <- as.matrix(heatmap_draw)
heatmap_draw <- log2(heatmap_draw+1)
for(i in 1:nrow(heatmap_draw)){
  heatmap_draw[i,]=scale(heatmap_draw[i,])
}
write.table(heatmap_draw,file="result/heatmap.txt",sep="\t",col.names=T,quote=F)


###########volcanol####################
pdf(file="result/vol.pdf",height=5,width=5)
# outTab[outTab==Inf]<-NA
# outTab[is.na(outTab)] <- 0
volcanol  <- read.table("data/EV_FTCvsFTA.edgeR.glmlrt.txt",sep="\t",check.names=F)

xMax=max(abs(as.numeric(as.vector(volcanol$log2FoldChange))))
yMax=max(-log10(volcanol$padj))+1
plot(as.numeric(as.vector(volcanol$log2FoldChange)), -log10(volcanol$padj), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano of FTC/FTA differential analysis", ylim=c(0,yMax),xlim=c(-3,3),yaxs="i",pch=20, cex=0.8)
diffSub=subset(volcanol, padj<fdrFilter & as.numeric(as.vector(log2FoldChange))>logFCfilter)
points(as.numeric(as.vector(diffSub$log2FoldChange)), -log10(diffSub$padj), pch=20, col="red",cex=0.8)
diffSub=subset(volcanol, padj<fdrFilter & as.numeric(as.vector(log2FoldChange))< -logFCfilter)
points(as.numeric(as.vector(diffSub$log2FoldChange)), -log10(diffSub$padj), pch=20, col="green",cex=0.8)

abline(v=0,lty=2,lwd=2)
abline(h=1.3,lty=2,lwd=1)

dev.off()

pdf(file="result/FTC.pdf",height=5,width=5)
# outTab[outTab==Inf]<-NA
# outTab[is.na(outTab)] <- 0
volcanol  <- read.table("data/FTC_EVvscf.edgeR.glmlrt.txt",sep="\t",check.names=F)

xMax=max(abs(as.numeric(as.vector(volcanol$log2FoldChange))))
yMax=max(-log10(volcanol$padj))+1
plot(as.numeric(as.vector(volcanol$log2FoldChange)), -log10(volcanol$padj), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano of FTC EV/cf differential analysis", ylim=c(0,40),xlim=c(-8,8),yaxs="i",pch=20, cex=0.8)
diffSub=subset(volcanol, padj<fdrFilter & as.numeric(as.vector(log2FoldChange))>logFCfilter)
points(as.numeric(as.vector(diffSub$log2FoldChange)), -log10(diffSub$padj), pch=20, col="red",cex=0.8)
diffSub=subset(volcanol, padj<fdrFilter & as.numeric(as.vector(log2FoldChange))< -logFCfilter)
points(as.numeric(as.vector(diffSub$log2FoldChange)), -log10(diffSub$padj), pch=20, col="green",cex=0.8)

abline(v=0,lty=2,lwd=2)
abline(h=1.3,lty=2,lwd=1)

dev.off()

pdf(file="result/FTA.pdf",height=5,width=5)
# outTab[outTab==Inf]<-NA
# outTab[is.na(outTab)] <- 0
volcanol  <- read.table("data/FTA_EVvscf.edgeR.glmlrt.txt",sep="\t",check.names=F)

xMax=max(abs(as.numeric(as.vector(volcanol$log2FoldChange))))
yMax=max(-log10(volcanol$padj))+1
plot(as.numeric(as.vector(volcanol$log2FoldChange)), -log10(volcanol$padj), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano of FTA EV/cf differential analysis", ylim=c(0,40),xlim=c(-8,8),yaxs="i",pch=20, cex=0.8)
diffSub=subset(volcanol, padj<fdrFilter & as.numeric(as.vector(log2FoldChange))>logFCfilter)
points(as.numeric(as.vector(diffSub$log2FoldChange)), -log10(diffSub$padj), pch=20, col="red",cex=0.8)
diffSub=subset(volcanol, padj<fdrFilter & as.numeric(as.vector(log2FoldChange))< -logFCfilter)
points(as.numeric(as.vector(diffSub$log2FoldChange)), -log10(diffSub$padj), pch=20, col="green",cex=0.8)

abline(v=0,lty=2,lwd=2)
abline(h=1.3,lty=2,lwd=1)

dev.off()

######waterfall######
i <- as.numeric(outTab$logFC)
i <- i[order(i)]
i[which(i == '-Inf')] = min(i[which(i != '-Inf')])
i[which(i == 'Inf')] = max(i[which(i != 'Inf')] )
barplot(-i,ylim = c(-6,4),main ='FTC',ylab = 'log2(FC)(EV/cf)')
abline(h = c(1,-2), col = "gray60",lwd = 2, lty = 2)
text(300,3, "FC > 2",cex = 1.5)
text(300,-3, "FC < -4",cex = 1.5)


######maketable######
######Venn######
FTC_EV <- outDiff[outDiff$logFC < 0,]$gene
FTC_cf <- outDiff[outDiff$logFC > 0,]$gene

FTA_EV <- outDiff[outDiff$logFC < 0,]$gene
FTA_cf <- outDiff[outDiff$logFC > 0,]$gene

length(intersect(FTC_EV,FTA_EV))
length(intersect(FTC_cf,FTA_cf))

list(FTC_EV,FTC_cf,FTA_EV,FTA_cf)

write.table(FTC_EV,file="FTC_EV.txt",sep="\t",col.names=F,quote=F,row.names = F)
write.table(FTC_cf,file="FTC_cf.txt",sep="\t",col.names=F,quote=F,row.names = F)
write.table(FTA_EV,file="FTA_EV.txt",sep="\t",col.names=F,quote=F,row.names = F)
write.table(FTA_cf,file="FTA_cf.txt",sep="\t",col.names=F,quote=F,row.names = F)

















#####rerank FTA_spec####
rownames(outDiff)=outDiff[,1]  #取出第一???
outDiff=outDiff[,-1]
tmppos <- c() #cf
tmpneg <- c() #EV
for(i in rownames(hmExp_spec_FTA)){
  if(outDiff[i,3]>0){
    tmppos <- append(tmppos,i)
  }
  else{
    tmpneg <-append(tmpneg,i)
  }
}

tmp <- append(tmppos,tmpneg)

hmExp_spec_FTA_rank <- hmExp_spec_FTA[tmp,]

#####rerank FTC_spec####
rownames(outDiff)=outDiff[,1]  #取出第一???
outDiff=outDiff[,-1]
tmppos2 <- c() #cf
tmpneg2 <- c() #EV
for(i in rownames(hmExp_spec_FTC)){
  if(outDiff[i,3]>0){
    tmppos2 <- append(tmppos2,i)
  }
  else{
    tmpneg2 <-append(tmpneg2,i)
  }
}

tmp2 <- append(tmppos2,tmpneg2)

hmExp_spec_FTC_rank <- hmExp_spec_FTC[tmp2,]

saigo_hm <- rbind(rbind(hmExp_spec_FTC_rank,hmExp_spec_FTA_rank),cbind(hmExp_rank_FTC[common_miRNA,],hmExp_rank_FTA[common_miRNA,]))
write.table(saigo_hm,file="saigo_heatmap.txt",sep="\t",col.names=T,quote=F)


########FTA/FTC diff gene analysis###########
diff=read.table('data/FTCvsFTA_up_down.txt',sep="\t",header=T,check.names=F)

rownames(outTab)=outTab[,1]  #取出第一???
outTab=outTab[,-1]
outTab[rownames(diff),]
diff_order <- c() #cf
for(i in rownames(diff)){
  
  
  if(i %in% common_miRNA){
    #print(hmExp_rank_FTC[i,])
    print('common')
    print(i)
  }
  if(i %in% FTC_EV){
    #print(hmExp_rank_FTC[i,])
    print('FTC_EV')
    print(i)
  }
  if(i %in% FTC_cf){
    #print(hmExp_rank_FTC[i,])
    print('FTC_cf')
    print(i)
  }
  
  if(i %in% FTA_EV){
    #print(hmExp_rank_FTC[i,])
    print('FTA_EV')
    print(i)
    diff_order <- append(diff_order,i)
  }
  
  if(i %in% FTA_cf){
    #print(hmExp_rank_FTA[i,])
    print('FTA_cf')
    print(i)
  }
}
diff_order <- append(diff_order,"hsa-miR-134-5p")
diff_order <- append(diff_order,"hsa-miR-146a-5p")
diff_order <- append(diff_order,"hsa-miR-223-5p")
diff_order <- append(diff_order,"hsa-miR-4306")
diff_order <- append(diff_order,"hsa-miR-432-5p")

diff_FTC <- data[diff_order,]
for(i in 1:nrow(diff_FTC)){
  diff_FTC[i,]=scale(diff_FTC[i,])
}


diff_FTA <- data[diff_order,]


for(i in 1:nrow(diff_FTA)){
  diff_FTA[i,]=scale(diff_FTA[i,])
}

diff_total <- cbind(diff_FTC,diff_FTA)
write.table(diff_total,file="diff_total.txt",sep="\t",col.names=T,quote=F)

data[rownames(diff),]


write.table(row.names(pf_FTA_EV),file="FTA_EV.txt",sep="\t",col.names=F,row.names = F,quote=F)
write.table(row.names(pf_FTA_cf),file="FTA_cf.txt",sep="\t",col.names=F,row.names = F,quote=F)
write.table(row.names(pf_FTC_EV),file="FTC_EV.txt",sep="\t",col.names=F,row.names = F,quote=F)
write.table(row.names(pf_FTC_cf),file="FTC_cf.txt",sep="\t",col.names=F,row.names = F,quote=F)

