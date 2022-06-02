# test RBP & motif
# last 220118 by bpf
# b.p.f@qq.com

setwd("C:/study/lulab_research/毕设/zq/miRNA")
options(scipen = 4,stringsAsFactors = F,digits = 5)
library('mice')


# def func. ---------------------------------------------------------------
pMiss <- function(x){round(sum(is.na(x))/length(x),3)}
## tpm
tpm <- function(count,gene.len){
  mat <- as.matrix(count)
  matrix_tpm <- 1000*mat / gene.len
  matrix_tpm <- t(t(matrix_tpm) * 1e6 / colSums(matrix_tpm))
  return(matrix_tpm)
}

cpm <- function(count){
  mat <- as.matrix(count)
  matrix_cpm <- mat
  matrix_cpm <- t(t(matrix_cpm) * 1e6 / colSums(matrix_cpm))
  return(matrix_cpm)
}

## def limma diff func.
#with weight
limma.trend <- function(logcpm,group){
  #logcpm = mat.norm.tmp
  #group = sample.table.tmp$source 
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(edgeR))
  #y <- DGEList(counts=mat, samples=samples, group=group)
  #y <- calcNormFactors(y, method=norm_method)
  model <- model.matrix(~group)
  #y <- voom(y, model, plot=FALSE)
  #fit <- lmFit(y, model)
  
  ## add weight
  we <- limma::arrayWeights(logcpm, design = model, method = "genebygene")
  fit <- lmFit(logcpm, model, weights = we)
  fit <- eBayes(fit, robust=TRUE, trend=TRUE) # , trend=TRUE (limma.trend), voom比limma-trend更适用于样本库大小不一的情况
  #fit2 <- contrasts.ft(fit)
  #fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
  #top_table <- topTable(fit2, sort.by='none', n=Inf)
  top_table <- topTable(fit, coef=2, sort.by='none', n=Inf)
  # rename columns
  mapped_names <- colnames(top_table)
  for(i in 1:ncol(top_table)){
    if(colnames(top_table)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(top_table)[i] == 'P.Value'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(top_table)[i] == 'adj.P.Val') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(top_table)[i]
    }
  }
  colnames(top_table) <- mapped_names
  res <- top_table
  return(res)
}

## limma.trend.paired deprecated !!!
limma.trend.paired <- function(logcpm,group,patient){
  #logcpm <- mat.norm.tmp
  #group <- sample.table.tmp$source
  #patient <- sample.table.tmp$cell_id
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(edgeR))
  #y <- DGEList(counts=mat, samples=samples, group=group)
  #y <- calcNormFactors(y, method=norm_method)
  model <- model.matrix(~patient+group)
  #y <- voom(y, model, plot=FALSE)
  #fit <- lmFit(y, model)
  fit <- lmFit(logcpm, model)
  fit <- eBayes(fit, robust=TRUE, trend=TRUE) # , trend=TRUE (limma.trend), voom比limma-trend更适用于样本库大小不一的情况
  #fit2 <- contrasts.ft(fit)
  #fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
  #top_table <- topTable(fit2, sort.by='none', n=Inf)
  top_table <- topTable(fit, coef=ncol(model), sort.by='none', n=Inf)
  # rename columns
  mapped_names <- colnames(top_table)
  for(i in 1:ncol(top_table)){
    if(colnames(top_table)[i] == 'logFC'){
      mapped_names[i] <- 'log2FoldChange'
    }else if(colnames(top_table)[i] == 'P.Value'){
      mapped_names[i] <- 'pvalue'
    }else if(colnames(top_table)[i] == 'adj.P.Val') {
      mapped_names[i] <- 'padj'
    }else{
      mapped_names[i] <- colnames(top_table)[i]
    }
  }
  colnames(top_table) <- mapped_names
  res <- top_table
  return(res)
}

# limma.voom <- function(logcpm,group){
#   suppressPackageStartupMessages(library(limma))
#   suppressPackageStartupMessages(library(edgeR))
#   #y <- DGEList(counts=mat, samples=samples, group=group)
#   #y <- calcNormFactors(y, method=norm_method)
#   model <- model.matrix(~group)
#   #y <- voom(y, model, plot=FALSE)
#   #fit <- lmFit(y, model)
#   fit <- lmFit(logcpm, model)
#   fit <- eBayes(fit, robust=TRUE) # , trend=TRUE (limma.trend), voom比limma-trend更适用于样本库大小不一的情况
#   #fit2 <- contrasts.ft(fit)
#   #fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
#   #top_table <- topTable(fit2, sort.by='none', n=Inf)
#   top_table <- topTable(fit, coef=2, sort.by='none', n=Inf)
#   # rename columns
#   mapped_names <- colnames(top_table)
#   for(i in 1:ncol(top_table)){
#     if(colnames(top_table)[i] == 'logFC'){
#       mapped_names[i] <- 'log2FoldChange'
#     }else if(colnames(top_table)[i] == 'P.Value'){
#       mapped_names[i] <- 'pvalue'
#     }else if(colnames(top_table)[i] == 'adj.P.Val') {
#       mapped_names[i] <- 'padj'
#     }else{
#       mapped_names[i] <- colnames(top_table)[i]
#     }
#   }
#   colnames(top_table) <- mapped_names
#   res <- top_table
#   return(res)
# }


## define deseq2 diff func.
#paired patient
diff <- function(mat,samples,group,patient,method,norm_method, filterType="small"){
  suppressPackageStartupMessages(library(edgeR))
  mat <- mat[,samples]
  print(dim(mat))
  y <- DGEList(counts=mat, samples=samples, group=group)
  
  ## filter low expr
  message(filterType)
  counts <- edgeR::getCounts(y)
  if (filterType=="small"){
    keep <- rowSums(counts>=1) >= 0.5*length(samples)
  }else if(filterType=="long"){
    tpm_mat <- tpm(count = counts, 
                   gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[2]))))
    keep <- rowSums(tpm_mat >= 1) >= 0.5*length(samples)
  }
  #keep <- filterByExpr(y,min.count = 1, min.total.count = 5, min.prop = 0.1) 
  #min.count: Minimum count required for at least some samples.
  #min.total.count: Minimum total count required.
  #min.prop: Minimum proportion of samples in the smallest group that express the gene.
  y <- y[keep,,keep.lib.sizes=FALSE]
  print(dim(y))
  
  y <- calcNormFactors(y, method=norm_method)
  
  design <- model.matrix(~ factor(patient) + factor(group))  # partial paired mode
  #design <- model.matrix(~  group + patient)  # partial paired mode, design order has great impact
  #design <- model.matrix(~ group) 
  y <- estimateDisp(y, design) 
  if(method == 'edger_glmqlf'){
    fit <- glmQLFit(y, design)
    test <- glmQLFTest(fit, coef=ncol(design))   
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_glmlrt'){
    fit <- glmFit(y, design)
    test <- glmLRT(fit, coef=ncol(design))  # coef should be the pos col
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_exact'){
    test <- exactTest(y)
    res <- topTags(test, n=nrow(mat), sort.by='none')
  }
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
  return(res)
}

#not paired
diff.v2 <- function(mat,samples,group,method,norm_method, filterType="small"){
  suppressPackageStartupMessages(library(edgeR))
  mat <- mat[,samples]
  print(dim(mat))
  y <- DGEList(counts=mat, samples=samples, group=group)
  
  ## filter low expr
  #keep <- filterByExpr(y,min.count = 1, min.total.count = 5, min.prop = 0.1)
  message(filterType)
  counts <- edgeR::getCounts(y)
  if (filterType=="small"){
    keep <- rowSums(counts>=1) >= 0.5*length(samples)
  }else if(filterType=="long"){
    tpm_mat <- tpm(count = counts, 
                   gene.len = as.numeric(unlist(lapply(strsplit(rownames(counts), "|", fixed = T), function(x) x[2]))))
    keep <- rowSums(tpm_mat >= 1) >= 0.5*length(samples)
  }
  #min.count: Minimum count required for at least some samples.
  #min.total.count: Minimum total count required.
  #min.prop: Minimum proportion of samples in the smallest group that express the gene.
  y <- y[keep,,keep.lib.sizes=FALSE]
  print(dim(y))
  
  y <- calcNormFactors(y, method=norm_method)
  
  #design <- model.matrix(~ patient + group)  # partial paired mode
  #design <- model.matrix(~  group + patient)  # partial paired mode, design order has great impact
  design <- model.matrix(~ factor(group))
  y <- estimateDisp(y, design)
  if(method == 'edger_glmqlf'){
    fit <- glmQLFit(y, design)
    test <- glmQLFTest(fit, coef=)
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_glmlrt'){
    fit <- glmFit(y, design)
    test <- glmLRT(fit, coef=2) # coef should be the pos col, ncol(design)=2
    res <- topTags(test, n=nrow(mat), sort.by='none')
  } else if(method == 'edger_exact'){
    test <- exactTest(y)
    res <- topTags(test, n=nrow(mat), sort.by='none')
  }
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
  return(res)
}


## diff wilcox paired 
wilcox.paired <- function(matrix_cpm,group,paired=TRUE){
  # matrix_cpm = mat.norm.tmp
  # group = sample.table.tmp$source
  # paired = FALSE
  suppressPackageStartupMessages(library(edgeR))
  test_func <- function(x){
    wilcox.test(as.numeric(x[group == levels(group)[1]]), as.numeric(x[group == levels(group)[2]]), alternative='two.sided', paired = paired)$p.value
  }
  #matrix_cpm <- cpm(mat)
  #table(group == levels(group)[2])
  #wilcox.test(as.numeric(matrix_cpm[1,group == levels(group)[1]]), as.numeric(matrix_cpm[1,group == levels(group)[2]]), alternative='two.sided', paired = paired)$p.value
  #matrix_cpm
  pvalues <- apply(matrix_cpm, 1, test_func)
  treatMeans <- apply(matrix_cpm[,which(group == levels(group)[2])], 1, mean)
  ctrlMeans <- apply(matrix_cpm[,which(group == levels(group)[1])], 1, mean)
  logFC <- treatMeans - ctrlMeans
  res <- data.frame(log2FoldChange=logFC,
                    pvalue=pvalues,
                    padj=p.adjust(pvalues, method='BH'),
                    baseMean=apply(matrix_cpm, 1, mean),
                    treatMean=treatMeans,
                    ctrlMean=ctrlMeans)
  return(res)
}



# test public data: 2021,Nature -------------------------------------------
# smRNA-seq CPM
# mat <- read.table("./data/2022-Nature-EVmotif/41586_2021_4234_MOESM14_ESM-Supplementary-Table-10.txt",header = T,sep = "\t",row.names = 1,check.names = T)
# colnames(mat) <- gsub(".","_",fixed = T,colnames(mat))
# colnames(mat) <- sub("_",".",fixed = T,colnames(mat))
# colnames(mat)[c(1,5,9,13)] <- paste0(colnames(mat)[c(1,5,9,13)],"_0")
# colnames(mat) <- gsub("Cell","Lysate",fixed = T,colnames(mat))
# ## log transform (needed for limma, not needed for wilcox, but both doesnot hurt)
# mat <- log2(mat+1)
# mat.norm[1:3,1:3]

## normlize
#mat.norm <- apply(mat,2,function(x) x/exp(mean(log(x)))) #  geometric mean CPM normalized
#mat.norm <- as.data.frame(t(mat.norm))


# qPCR Reltive Ct
mat <- read.table("../paper_data/rawdata.txt",header = T,sep = "\t",row.names = 1,check.names = F)

## normlize
mat.norm <- apply(mat,1,function(x) x-as.numeric(mat["mmu-miR-501-5p",])) # miR-138-5p or miR-501-5p
mat.norm <- as.data.frame(t(mat.norm))
write.table(mat.norm, '../paper_data/normdata.txt',sep = "\t")
#fill NA ?
#mat.norm <- na.omit(mat.norm)

## get sample table
sample.table <- as.data.frame(cbind(data_id=colnames(mat),
                                    cell_id=unlist(lapply(strsplit(colnames(mat),".",fixed=T),function(x) x[2])),
                                    source=unlist(lapply(strsplit(colnames(mat),".",fixed=T),function(x) x[1]))
                                    #source=unlist(lapply(strsplit(colnames(mat),"_",fixed=T),function(x) x[1]))
                                    ))
rownames(sample.table) <- sample.table$data_id
sample.table$cell_line <- unlist(lapply(strsplit(sample.table$cell_id,"_",fixed=T),function(x) x[1]))
sample.table$source <- factor(sample.table$source,levels = c("Lysate","Exosome"))
sample.table$cell_id <- factor(sample.table$cell_id)

## filter unpaired samples (only for qPCR)
#sample.table <- sample.table[!(sample.table$cell_id %in% c("BAT_4","3T3L1_4","AML12_4","SVEC_4")),]

## run limma diff
#mat has >600 genes, trend=T, robust=T might be better (https://support.bioconductor.org/p/p134303/)
#diff <- limma.voom(logcpm = mat.norm, group = sample.table$source )
#diff <- limma.trend(logcpm = mat.norm, group = sample.table$source )

res <- list()
mymaxit <- 5
for (i in unique(sample.table$cell_line)){
  #i <- "C2C12"
  print(i)
  sample.table.tmp <- sample.table[sample.table$cell_line==i,]
  #mat.norm.tmp <- mat.norm[,sample.table.tmp$data_id]
  mat.norm.tmp <- mat[,sample.table.tmp$data_id] # do not use normalized mat
  print(nrow(mat.norm.tmp))
  
  # omit NA rows
  # mat.norm.tmp <- na.omit(mat.norm.tmp)
  # print(nrow(mat.norm.tmp))
  mat.norm.tmp[is.na(mat.norm.tmp)] <- 0
  
  # mat.norm.tmp <- mat.norm.tmp[apply(mat.norm.tmp, 1, pMiss) < 0.15,]
  # if(i == 'C2C12'){
  #   tempData1 <- mice(mat.norm.tmp[,1:4],m=5,maxit=mymaxit,meth='pmm',seed=500)
  #   tempData2 <- mice(mat.norm.tmp[,5:8],m=5,maxit=mymaxit,meth='pmm',seed=500)
  #   mat.norm.tmp <- cbind(complete(tempData1,action = 5),complete(tempData2,action = 5))
  # }
  # else{
  #   tempData1 <- mice(mat.norm.tmp[,1:3],m=5,maxit=mymaxit,meth='pmm',seed=500)
  #   tempData2 <- mice(mat.norm.tmp[,4:7],m=5,maxit=mymaxit,meth='pmm',seed=500)
  #   mat.norm.tmp <- cbind(complete(tempData1,action = 5),complete(tempData2,action = 5))
  # }

  diff.tmp <- limma.trend(logcpm = mat.norm.tmp, group = sample.table.tmp$source )
  #diff.tmp <- limma.voom(logcpm = mat.norm.tmp, group = sample.table.tmp$source )
  #diff.tmp <- limma.trend.paired(logcpm = mat.norm.tmp, group = factor(sample.table.tmp$source), patient = factor(sample.table.tmp$cell_id) )
  
  #diff.tmp <- wilcox.paired(matrix_cpm = mat.norm.tmp,group = sample.table.tmp$source,paired = FALSE) # need omit na rows for wilcox test
  res[[i]] <- diff.tmp
}

for (i in unique(sample.table$cell_line)){
print(i)

res.tmp <- res[[i]]
# write(rownames(res.tmp[res.tmp$padj<0.1 & res.tmp$log2FoldChange>0,]),paste('result/',i,'_EV.txt',sep=''))
# write(rownames(res.tmp[res.tmp$padj<0.1 & res.tmp$log2FoldChange<0,]),paste('result/',i,'_cf.txt',sep=''))

print(table(res.tmp$padj<0.1 & res.tmp$log2FoldChange>=0))
print(table(res.tmp$padj<0.1 & res.tmp$log2FoldChange<=0))
}
#195 TRUE, trend
#188 TRUE, voom
#195 TRUE, trend
#365 TRUE, 194 cell, 171 EV, voom
#366 TRUE, 200 cell, 166 EV, trend

#suppl tab 6
#205 cell,173 EV

res.tmp <- res[["C2C12"]]
write.table(na.omit(rownames(res.tmp)[res.tmp$padj<0.1 & res.tmp$log2FoldChange>0]),"./tmp.txt",quote = F,sep = "\t",row.names = F,col.names = F)



# get ev enriched mir
res <- rio::import("./data/2022-Nature-EVmotif/41586_2021_4234_MOESM10_ESM-Supplementary-Table-6.xlsx",)
rownames(res) <- res$...1
res <- res[,-1]
cellline <- c("BAT","3T3L1","C2C12","SVEC","AML12")
cellline.col <- paste0("sEV.vs.Lysate.in.",cellline,".sig")
colnames(res)

#for (i in 1:length(cellline)){
#  print(i)
#  tmp.ev <- rownames(res)[res[[cellline[i]]]==1]
#}

ev5 <-  rownames(res)[rowSums(res[,1:5]==1)==5]
ev3.4 <-  rownames(res)[rowSums(res[,1:5]==1)==3 | rowSums(res[,1:5]==1)==4]
cl5 <-  rownames(res)[rowSums(res[,1:5]==-1)==5]
cl3.4 <-  rownames(res)[rowSums(res[,1:5]==-1)==3 | rowSums(res[,1:5]==-1)==4]
notenrich <- rownames(res)[rowSums(res[,1:5]==0)==5]
#those sorted into the sEV in all cells are shown in red (n=13); 
#those sorted into sEV in 3 or 4 of the five cell types are shown in green (n=90); 
#those not enriched in either sEV or cells are shown in black (n=109); 
#those retained in 3 or 4 cell types are shown in pink (n=97);  fig1h:2+14+34+4*4+3+2+10+14+2+3=100 
#those retained in the cell bodies of all cell types are shown in blue (n=43). 
l <- list(ev5,ev3.4,cl5,cl3.4,notenrich)
names(l) <- c("EV5","EV3-4","CL5","CL3-4","notenrich")
for (i in names(l)){
  write.table(l[[i]],paste0("./output/Nat2022/",i,".txt"),row.names = F,col.names = F,sep = "\t",quote = F)
}

  
  

# GC content
library(seqinr)
convGC <- function(seqs){
  #print(seqs)
  #as.character(mySequences$`mmu-miR-15b-5p`)
  #seqs <- "cuauacaaucuacugucuuucc"
  seqs <- gsub("U|u","t",seqs)  # important, seqinr::GC seems neglect U|u base
  tmp <- seqinr::GC(seqinr::s2c(seqs),forceToLower = TRUE, exact = FALSE, NA.GC = NA, 
             oldGC = T # not default
             #alphabet = s2c("acgtswmkryvhdb")
             )
  #cuauacaaucuacugucuuucc	22	0.363636364
  return(tmp)
}

getGC <- function(x){
  #x <- "CL3-4"
  print(x)
  mySequences <- Biostrings::readRNAStringSet(paste0("./myNat/",x,".fa"))
  #miR <- names(mySequences)
  res <- lapply(as.character(mySequences),convGC)
  #names(res) <- names(mySequences)
  res <- do.call("rbind",res)
  res2 <- as.data.frame(cbind(GC=res[,1],mir=names(mySequences),type=x))
  return(res2)
}

dat <- lapply(c("EV5","EV3-4","CL5","CL3-4","notenrich"),getGC)
dat <- do.call("rbind",dat)
dat$type <- factor(dat$type,levels = c("EV5","EV3-4","notenrich","CL3-4","CL5"))
dat$GC <- as.numeric(dat$GC)
table(dat$type)

b <- runif(nrow(dat), -0.1, 0.1)
ggplot(dat,aes(x=type,y=GC,group=type,fill=type))+
  #geom_bar(stat = "identity",position = "dodge", color="grey30") + # 
  geom_boxplot(fill="white",outlier.alpha = 0)+
  geom_point(aes(x=as.numeric(type)+b,y=GC,color=type,fill=type,shape=type),size=2)+
  ylim(c(0,1)) +
  #ggsci::scale_fill_jama() +
  #scale_color_brewer(palette = "Set3") +
  #scale_x_discrete(labels=RNA_label)+
  #scale_fill_manual(values=RNA_color,labels=RNA_label)+
  #geom_text(aes(x=type,y=roc+0.017,label=format(roc,digits=3)),size=6,angle=0,color="grey30")+
  labs(x="",y="miRNA %CG",title = "") + 
  scale_x_discrete(label=c("sEV-enriched 5 cell types","sEV-enriched 3-4 cell types","all unchanged","Cell-enriched 3-4 cell types","Cell-enriched 5 cell types"))+
  #ggsci::scale_fill_aaas()+
  scale_fill_manual(values = c("firebrick","seagreen","black","steelblue","purple3"))+
  scale_color_manual(values = c("firebrick","seagreen","black","steelblue","purple3"))+
  #scale_shape_manual(values = c(1,2,6,5,0))+
  scale_shape_manual(values = c(16,17,25,18,15))+
  ggpubr::stat_compare_means(ref.group = "notenrich",
                             #comparisons = my_comparisons,
                             label.x.npc = 0.5,hide.ns=F,size =8,#paired = TRUE,  
                             #aes(group = Group,label = p.format), # p.signif,p.format
                             label = "p.signif",
                             method = "wilcox.test",
                             #symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns"))
  ) +
  theme_classic() +
  theme(plot.title = element_text(size=32,color="black",hjust = 0.5),
        axis.title = element_text(size=24,color="black"), 
        axis.text = element_text(size=16,color="black"), # ,face="bold"
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 0.9, vjust = 0.5),
        panel.grid=element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold"),
        strip.background = element_blank() ) +
  #panel.spacing = unit(1.5, "lines"),
  #panel.border = element_rect(colour = "black", fill=NA, size=2)
  guides(color=guide_legend("Type"))
ggsave("./result/myGC.pdf",width = 12,height = 10)




# MFE:RNAfold
#library(LncFinder)
#mySequences <- seqinr::read.fasta("./output/NC2013/cell.fa")
#test <- LncFinder::run_RNAfold(mySequences, RNAfold.path = "/apps/bin/RNAfold -p -d2 --noLP ", parallel.cores = 6)
#length(mySequences)
#res has no MFE
#tmp <- seqinr::read.fasta("./output/Nat2022/CL3-4.txt.preMIR.RNAfold")
#tmp <- ips::read.fas("./output/Nat2022/CL3-4.txt.preMIR.RNAfold")
#dat <- list()
conv <- function(x){
  #x <- "EV3-4"
  print(x)
  #rnafold <- read.delim2(paste0("./output/Nat2022/",x,".txt.preMIR.RNAfold"),header = F,sep = " ") # CL3-4, CL5, EV5, EV3-4, notenrich
  rnafold <- read.delim2(paste0("./myNat/",x,".RNAfold"),header = F,sep = "") # CL3-4, CL5, EV5, EV3-4, notenrich
  #rnafold <- read.delim2(paste0("./Nat2022/",x,".txt.RNAfold"),header = F,sep = "") # CL3-4, CL5, EV5, EV3-4, notenrich
  
  mir <- rnafold$V1[grepl(">",rnafold$V1)]
  mir <- gsub(">","",mir)
  mfe <- rnafold$V3[rnafold$V3!=""]
#  mfe <- gsub("(","",mfe,fixed = T)
#  mfe <- gsub(")","",mfe,fixed = T)
  mfe <- gsub("(","",mfe,fixed = T)
  mfe <- gsub(")","",mfe,fixed = T)
  mfe <- as.numeric(mfe)
  rnafold.out <- as.data.frame(cbind(mir=mir,mfe=mfe))
  rnafold.out$mfe <- as.numeric(rnafold.out$mfe)
  #colnames(rnafold.out)[2] <- x
  rnafold.out$type <- x
  #max(rnafold.out$mfe)
  #rnafold.out <- rnafold.out$mfe-max(rnafold.out$mfe)
  #hist(rnafold.out$mfe) # ,xlim = c(-40,0)
  #dat[[x]] <- rnafold.out
  return(rnafold.out)
}

dat <- lapply(c("EV5","EV3-4","CL5","CL3-4","notenrich"),conv)
dat <- do.call("rbind",dat)
dat$type <- factor(dat$type,levels = c("EV5","EV3-4","notenrich","CL3-4","CL5"))

b <- runif(nrow(dat), -0.1, 0.1)
ggplot(dat,aes(x=type,y=mfe,group=type,fill=type))+
  #geom_bar(stat = "identity",position = "dodge", color="grey30") + # 
  geom_boxplot(fill="white",outlier.alpha = 0)+
  geom_point(aes(x=as.numeric(type)+b,y=mfe,color=type,fill=type,shape=type),size=2)+
  #ylim(c(0,1)) +
  #ggsci::scale_fill_jama() +
  #scale_color_brewer(palette = "Set3") +
  #scale_x_discrete(labels=RNA_label)+
  #scale_fill_manual(values=RNA_color,labels=RNA_label)+
  #geom_text(aes(x=type,y=roc+0.017,label=format(roc,digits=3)),size=6,angle=0,color="grey30")+
  labs(x="",y="delta G (kcal/mol)",title = "") + 
  scale_x_discrete(label=c("sEV-enriched 5 cell types","sEV-enriched 3-4 cell types","all unchanged","Cell-enriched 3-4 cell types","Cell-enriched 5 cell types"))+
  #ggsci::scale_fill_aaas()+
  scale_fill_manual(values = c("firebrick","seagreen","black","purple3","steelblue"))+
  scale_color_manual(values = c("firebrick","seagreen","black","purple3","steelblue"))+
  scale_shape_manual(values = c(16,17,25,18,15))+
  ggpubr::stat_compare_means(ref.group = "notenrich",
                              #comparisons = my_comparisons,
                             label.x.npc = 0.5,hide.ns=F,size =8,#paired = TRUE,  
                             #aes(group = Group,label = p.format), # p.signif
                             label = "p.signif",
                             method = "wilcox.test",
                             #symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns"))
  ) +
  theme_classic() +
  theme(plot.title = element_text(size=32,color="black",hjust = 0.5),
        axis.title = element_text(size=24,color="black"), 
        axis.text = element_text(size=16,color="black"), # ,face="bold"
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 0.9, vjust = 0.5),
        panel.grid=element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold"),
        strip.background = element_blank() ) +
  #panel.spacing = unit(1.5, "lines"),
  #panel.border = element_rect(colour = "black", fill=NA, size=2)
  guides(color=guide_legend("Type"))
ggsave("result/MFE_mine.pdf",width = 12,height = 10)






  
  





## reply qpcr array
options(stringsAsFactors = FALSE)

library(limma)
library(ezlimma)
library(ezlimmaplot)
tab1 <- read.csv("./data/2022-Nature-EVmotif/41586_2021_4234_MOESM5_ESM-Supplementary-Table-1.csv", row.names = 1)

grp <- gsub(pattern = "_.", replacement = "", x=colnames(tab1))
grp <- factor(grp, levels=unique(grp))
design <- model.matrix(~0+grp)
colnames(design) <- gsub("grp", "", colnames(design))
aw <- limma::arrayWeights(tab1, design = design, method = "genebygene")
names(aw) <- colnames(tab1)

contr.v <- c(Exosome.vs.Lysate.in.3T3L1 = "Exosome.3T3L1 - Lysate.3T3L1", 
             Exosome.vs.Lysate.in.C2C12 = "Exosome.C2C12 - Lysate.C2C12",
             Exosome.vs.Lysate.in.SVEC = "Exosome.SVEC - Lysate.SVEC", 
             Exosome.vs.Lysate.in.AML12 = "Exosome.AML12 - Lysate.AML12",  
             Exosome.vs.Lysate.in.BAT = "Exosome.BAT - Lysate.BAT")
mtt <- ezlimma::limma_contrasts(tab1, grp = grp, contrast.v = contr.v, # weights = aw, 
                                cols=c("t", "P.Value", "adj.P.Val", "logFC"))

enrich.tab <- ezlimmaplot::ezvenn(tab=mtt, fdr.cutoff = 0.1, plot=FALSE)
fig1g <- apply(enrich.tab, 2, FUN=function(xx) summary(as.factor(xx)))
rownames(fig1g) <- c("Cell", "Not enriched", "sEV")
colnames(fig1g) <- gsub("\\.sig", "",  gsub("Exosome.vs.Lysate.in.", "", colnames(fig1g), fixed = TRUE))
write.csv(t(fig1g[, c("BAT", "C2C12", "3T3L1", "AML12", "SVEC")]), "fig_1g.csv")
fig1g <- t(fig1g)
fig1g


## reply (test smRNA-seq)
mat <- read.table("./data/2022-Nature-EVmotif/41586_2021_4234_MOESM14_ESM-Supplementary-Table-10.txt",header = T,sep = "\t",row.names = 1,check.names = T)
colnames(mat) <- gsub(".","_",fixed = T,colnames(mat))
colnames(mat) <- sub("_",".",fixed = T,colnames(mat))
colnames(mat)[c(1,5,9,13)] <- paste0(colnames(mat)[c(1,5,9,13)],"_0")
colnames(mat) <- gsub("Cell","Lysate",fixed = T,colnames(mat))
## log transform (needed for limma, not needed for wilcox, but both doesnot hurt)
mat <- log2(mat+1)

options(stringsAsFactors = FALSE)
library(limma)
library(ezlimma)
library(ezlimmaplot)
#tab1 <- read.csv("./data/2022-Nature-EVmotif/41586_2021_4234_MOESM5_ESM-Supplementary-Table-1.csv", row.names = 1)
tab1 <- mat
colnames(mat)
grp <- gsub(pattern = "_.", replacement = "", x=colnames(tab1))
grp <- factor(grp, levels=unique(grp))
design <- model.matrix(~0+grp)
colnames(design) <- gsub("grp", "", colnames(design))
aw <- limma::arrayWeights(tab1, design = design, method = "genebygene")
names(aw) <- colnames(tab1)

contr.v <- c(Exosome.vs.Lysate.in.3T3L1 = "Exosome.3T3L1 - Lysate.3T3L1", 
             #Exosome.vs.Lysate.in.C2C12 = "Exosome.C2C12 - Lysate.C2C12",
             #Exosome.vs.Lysate.in.SVEC = "Exosome.SVEC - Lysate.SVEC", 
             Exosome.vs.Lysate.in.AML12 = "Exosome.AML12 - Lysate.AML12"
             #Exosome.vs.Lysate.in.BAT = "Exosome.BAT - Lysate.BAT"
             )
mtt <- ezlimma::limma_contrasts(tab1, grp = grp, contrast.v = contr.v, weights = aw, trend = T,
                                cols=c("t", "P.Value", "adj.P.Val", "logFC"))

enrich.tab <- ezlimmaplot::ezvenn(tab=mtt, fdr.cutoff = 0.1, plot=FALSE)
fig1g <- apply(enrich.tab, 2, FUN=function(xx) summary(as.factor(xx)))
rownames(fig1g) <- c("Cell", "Not enriched", "sEV")
colnames(fig1g) <- gsub("\\.sig", "",  gsub("Exosome.vs.Lysate.in.", "", colnames(fig1g), fixed = TRUE))
#write.csv(t(fig1g[, c("BAT", "C2C12", "3T3L1", "AML12", "SVEC")]), "fig_1g.csv")
fig1g <- t(fig1g)
fig1g
tmp <- limma.trend(logcpm = tab1, group = grp)
table(tmp$padj<0.1 & tmp$log2FoldChange>0)
table(tmp$padj<0.1 & tmp$log2FoldChange<0)
(149+128)/459



# QC WHK-WHCSU sRNA&lRNA FTC&FTA  ----------------------------------------------
## small RNA
qc.cf.genome.map <- read.table("./output/lulab/FTC_cf_small/multiqc/map_genome/map_genome_multiqc_data/multiqc_bowtie2.txt",header = T,sep = "\t")
#uniq_reads: csFTC-21,csFTA-4
#map_ratio: csFTC-20,csFTA-6
qc.EV.genome.map <- read.table("./output/lulab/FTC_EV_small/multiqc/map_genome/map_genome_multiqc_data/multiqc_bowtie2.txt",header = T,sep = "\t")
#uniq_reads: FTA-3,FTA-20
#map_ratio<0.2: FTA-14

qc.cf <- data.frame(row.names = qc.cf.genome.map$Sample, mapping_ratio=qc.cf.genome.map$overall_alignment_rate, uniq_reads=qc.cf.genome.map$total_reads)
qc.cf$source <- "cf-small"
qc.EV <- data.frame(row.names = qc.EV.genome.map$Sample, mapping_ratio=qc.EV.genome.map$overall_alignment_rate, uniq_reads=qc.EV.genome.map$total_reads)
qc.EV$source <- "EV-small"

qc <- rbind(qc.cf,qc.EV)
rownames(qc) <- gsub("_map_genome","",rownames(qc))
qc$patient <- gsub("cs","",rownames(qc))
qc$uniq_reads <- qc$uniq_reads/1000000

## plot
dat <- qc
dat$group <- factor(dat$source)
#qc$patient <- factor(qc$patient)
dat$label <- unlist(lapply(strsplit(rownames(dat),"_"),function(x) x[1]))
b <- runif(nrow(dat), -0.1, 0.1)

library(ggplot2)
ggplot(dat,aes(x=group,y=mapping_ratio,group=group,fill=group))+
  #geom_bar(stat = "identity",position = "dodge", color="grey30") + # 
  geom_violin(alpha=0.4)+
  #geom_boxplot(aes(x=group,y=mapping_ratio,fill=group,color=group),fill="white",outlier.alpha = 0, alpha=0.8,width=0.8)+
  geom_point(aes(x=as.numeric(group)+b, y=mapping_ratio, group=group,fill=group,color=group))+
  geom_line(aes(x=as.numeric(group)+b, y=mapping_ratio, group=patient), linetype="dashed", alpha=0.6)+
  # ggrepel::geom_text_repel(aes(x=as.numeric(group)+b,y=mapping_ratio,label=label), #nudge_x=.5,nudge_y=.5,
  #                         #box.padding = 0.5,box.padding = .1,point.padding = .1,
  #                         size = 4,color="black",max.overlaps = Inf ) + 
  geom_hline(yintercept = 20,linetype="dashed",color="firebrick")+
  labs(x="",y="Genome Mapping Ratio %",title = "small RNA QC") + 
  #scale_x_discrete(label=c("sEV-enriched 5 cell types","sEV-enriched 3-4 cell types","all unchanged","Cell-enriched 3-4 cell types","Cell-enriched 5 cell types"))+
  ggsci::scale_fill_aaas()+
  ggsci::scale_color_aaas()+
  ggpubr::stat_compare_means(#ref.group = "IGF2BP1",
                             #comparisons = list(c("bg_cell","cell"),c("bg_exosome","exosome")),
                             label.x.npc = 0.5,hide.ns=F,size =8,#paired = TRUE,  
                             #aes(group = Group,label = p.format), # p.signif,p.format
                             label = "p.signif",
                             method = "wilcox.test",
                             #symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns"))
  ) +
  theme_classic() +
  theme(plot.title = element_text(size=24,color="black",hjust = 0.5),
        axis.title = element_text(size=24,color="black"), 
        axis.text = element_text(size=24,color="black"), # ,face="bold"
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 0.9, vjust = 0.5),
        panel.grid=element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold"),
        strip.background = element_blank() ) +
  #panel.spacing = unit(1.5, "lines"),
  #panel.border = element_rect(colour = "black", fill=NA, size=2)
  guides(color=guide_legend("Source"))
ggsave("./output/lulab/FTC_small/QC1.pdf")

ggplot(dat,aes(x=group,y=uniq_reads,group=group,fill=group))+
  #geom_bar(stat = "identity",position = "dodge", color="grey30") + # 
  geom_violin(alpha=0.4)+
  #geom_boxplot(aes(x=group,y=uniq_reads,fill=group,color=group),fill="white",outlier.alpha = 0.8, alpha=0.8,width=0.8)+
  geom_point(aes(x=as.numeric(group)+b, y=uniq_reads, group=group,fill=group,color=group))+
  geom_line(aes(x=as.numeric(group)+b, y=uniq_reads, group=patient), linetype="dashed", alpha=0.4)+
  # ggrepel::geom_text_repel(aes(x=as.numeric(group)+b,y=uniq_reads,label=label), #nudge_x=.5,nudge_y=.5,
  #                          #box.padding = 0.5,box.padding = .1,point.padding = .1,
  #                          size = 4,color="black",max.overlaps = Inf ) + 
  labs(x="",y="Genome Total Reads M",title = "small RNA QC") + 
  #geom_hline(yintercept = 20,linetype="dashed",color="firebrick")+
  #ylim(c(0,5))+
  #scale_x_discrete(label=c("sEV-enriched 5 cell types","sEV-enriched 3-4 cell types","all unchanged","Cell-enriched 3-4 cell types","Cell-enriched 5 cell types"))+
  ggsci::scale_fill_aaas()+
  ggsci::scale_color_aaas()+
  ggpubr::stat_compare_means(#ref.group = "IGF2BP1",
    #comparisons = list(c("bg_cell","cell"),c("bg_exosome","exosome")),
    label.x.npc = 0.5,hide.ns=F,size =8,#paired = TRUE,  
    #aes(group = Group,label = p.format), # p.signif,p.format
    label = "p.signif",
    method = "wilcox.test",
    #symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns"))
  ) +
  theme_classic() +
  theme(plot.title = element_text(size=24,color="black",hjust = 0.5),
        axis.title = element_text(size=24,color="black"), 
        axis.text = element_text(size=24,color="black"), # ,face="bold"
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 0.9, vjust = 0.5),
        panel.grid=element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold"),
        strip.background = element_blank() ) +
  #panel.spacing = unit(1.5, "lines"),
  #panel.border = element_rect(colour = "black", fill=NA, size=2)
  guides(color=guide_legend("Source"))
ggsave("./output/lulab/FTC_small/QC2.pdf")




## long RNA
setwd("/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/")
options(scipen = 4,stringsAsFactors = F,digits = 5)
cf <- read.table("output/cf-long/log/QC.txt",header = T,sep = "\t")
cf <- t(cf)
colnames(cf) <- cf[1,]
cf <- cf[-1,]
EV <- read.table("output/EV-long/log/QC.txt",header = T,sep = "\t")
EV <- t(EV)
colnames(EV) <- EV[1,]
EV <- EV[-1,]

all.table <- as.data.frame(rbind(cf,EV),colClasses=rep("numeric",35))
for(i in 1:ncol(all.table)){
  all.table[,i] <- as.numeric(all.table[,i])
}
# QC:
# Raw Filter	10000000.00 
# Clean Filter	5000000.00 
# Spike in 	0.50 
# Genome Aligned Reads	500000.00 
# rRNA	0.50 
# Long RNA ratio	0.20 
# Unclassified	0.30 
# Intron Spanning	100000.00 
#colnames(all.table)
failQC <- rownames(all.table)[all.table$raw<10000000.00 
                              | all.table$trimGC<5000000.00 
                              | all.table$spikein_long/all.table$trimGC>0.5 
                              | all.table$hg38_long<500000.00 
                              | all.table$rRNA/all.table$trimGC>0.5 
                              | (all.table$mRNA+all.table$lncRNA)/all.table$hg38_long_dedup<0.2 
                              | all.table$unclassified/all.table$hg38_long_dedup>0.3]
str(all.table )
dim(all.table)

hist(all.table$unclassified/all.table$hg38_long,breaks = 100)
table(all.table$rRNA/all.table$trimGC>0.5 )

all.table$unclassified.hg38_long.ratio <- all.table$unclassified/all.table$hg38_long_dedup
rownames(all.table)
all.table$group <- substr(rownames(all.table),1,2)
all.table$group <- ifelse(all.table$group=="cl","cf-long","EV-long")
all.table$group <- factor(all.table$group)
all.table$patient <- unlist(sapply(strsplit(rownames(all.table),"_",fixed = T),"[",1))
all.table$patient <- gsub("cl|el","",all.table$patient)
all.table$patient <- factor(all.table$patient)

#elFTA.13_L4
#elFTA.4_L4
#elFTA.6_L3
colnames(all.table)
dat <- all.table[,c("unclassified.hg38_long.ratio","group","patient")]
b <- runif(nrow(dat), -0.1, 0.1)

library(ggplot2)
ggplot(dat,aes(x=group,y=unclassified.hg38_long.ratio,group=group,fill=group))+
  geom_violin(alpha=0.4)+
  geom_point(aes(x=as.numeric(group)+b, y=unclassified.hg38_long.ratio, group=group,fill=group,color=group))+
  geom_line(aes(x=as.numeric(group)+b, y=unclassified.hg38_long.ratio, group=patient), linetype="dashed", alpha=0.6)+
  geom_hline(yintercept = 3,linetype="dashed",color="firebrick")+
  geom_hline(yintercept = 0.5,linetype="dashed",color="firebrick")+
  labs(x="",y="Unclassified/Genome Mapped",title = "long RNA QC") + 
  #scale_x_discrete(label=c("sEV-enriched 5 cell types","sEV-enriched 3-4 cell types","all unchanged","Cell-enriched 3-4 cell types","Cell-enriched 5 cell types"))+
  ggsci::scale_fill_aaas()+
  ggsci::scale_color_aaas()+
  ggpubr::stat_compare_means(#ref.group = "IGF2BP1",
    #comparisons = list(c("bg_cell","cell"),c("bg_exosome","exosome")),
    label.x.npc = 0.5,hide.ns=F,size =8,#paired = TRUE,  
    #aes(group = Group,label = p.format), # p.signif,p.format
    label = "p.signif",
    method = "wilcox.test",
    #symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns"))
  ) +
  theme_classic() +
  theme(plot.title = element_text(size=24,color="black",hjust = 0.5),
        axis.title = element_text(size=24,color="black"), 
        axis.text = element_text(size=24,color="black"), # ,face="bold"
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 0.9, vjust = 0.5),
        panel.grid=element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold"),
        strip.background = element_blank() ) +
  #panel.spacing = unit(1.5, "lines"),
  #panel.border = element_rect(colour = "black", fill=NA, size=2)
  guides(color=guide_legend("Source"))
ggsave("./output/FTC_long/QC1.pdf")
table(dat$unclassified.hg38_long.ratio>3)
dat[dat$unclassified.hg38_long.ratio>13,]
#elFTA.13_L4
#elFTA.17_L3                  
#elFTA.4_L4                  
#elFTA.6_L3





# diff WHK-WHCSU sRNA EVvsCF --------------------------------------------------------------------
# FTC (pair and unpair)
## read sample table
sample.table <- read.table("data/FTC_sample_table.txt", header = TRUE, check.names=FALSE, sep='\t') # row.names=1, 
sample.table$group <- ifelse(grepl("FTC",sample.table$patient_id),"FTC","FTA")
sample.table <- sample.table[,c("ev_small","cf_small","group","patient_id")]
sample.table <- reshape2::melt(sample.table,id.var=c("patient_id","group"))
colnames(sample.table)[3:4] <- c("source","data_id")
rownames(sample.table) <- sample.table$data_id

## read mat
#matrix <- "./data/lulab/FTC_EV_small/count_matrix/miRNA_count_matrix.txt"
#matrix <- "./data/lulab/FTC_EV_long/count_matrix/gencode.txt"
mat.cf.small <- read.table("data/cf_small/miRNA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
mat.ev.small <- read.table("data/EV_small/miRNA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
#colnames(mat.ev.small) <- paste0("es",colnames(mat.ev.small) )
table(rownames(mat.ev.small) == rownames(mat.cf.small))
mat.ev.small <- mat.ev.small[rownames(mat.cf.small),]
mat <- cbind(mat.cf.small,mat.ev.small)
s <- intersect(sample.table$data_id,colnames(mat))
mat <- mat[,s]
sample.table <- sample.table[s,]

## unpair mode 
positive_samples <- sample.table[sample.table$group=="FTC" & sample.table$source=="ev_small","data_id"]
negative_samples <- sample.table[sample.table$group=="FTC" & sample.table$source=="cf_small","data_id"]
samples <- c(positive_samples, negative_samples)
sample.table <- sample.table[samples,]
mat <- mat[,samples]
group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
patient <- sample.table$patient_id[sample.table$group=="FTC" & (sample.table$source=="ev_small" | sample.table$source=="cf_small")]
method <- "edger_glmlrt"
norm_method <- "TMM"
#res1 <- res
#table(patient)
res.FTC.unpair <- diff.v2(mat = mat, samples = samples, group = group, method = method, norm_method = norm_method)  #
write.table(res.FTC.unpair,"./output/result/unpaired/miR_diff_FTC_EVvsCF_unpair.txt",quote = F,sep = "\t",row.names = T,col.names = T)

res.FTC.unpair <- res.FTC.unpair[order(res.FTC.unpair$pvalue,decreasing = F),]
res.FTC.unpair.sig <- res.FTC.unpair[res.FTC.unpair$padj<=0.1,] # abs(res.FTC.unpair$log2FoldChange)>=2 & 
res.FTC.unpair.sig <- res.FTC.unpair.sig[order(-res.FTC.unpair.sig$pvalue,abs(res.FTC.unpair.sig$log2FoldChange),decreasing = T),]
table(res.FTC.unpair.sig$log2FoldChange>0)


## pair mode 
sample.table.pair <- sample.table[!(sample.table$patient_id %in% c("FTA-14","FTA-19")),]
positive_samples <- sample.table.pair[sample.table.pair$group=="FTC" & sample.table.pair$source=="ev_small","data_id"]
negative_samples <- sample.table.pair[sample.table.pair$group=="FTC" & sample.table.pair$source=="cf_small","data_id"]
samples <- c(positive_samples, negative_samples)
sample.table.pair <- sample.table.pair[samples,]
mat <- mat[,samples]
group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
patient <- sample.table.pair$patient_id[sample.table.pair$group=="FTC" & (sample.table.pair$source=="ev_small" | sample.table.pair$source=="cf_small")]
table(patient)

res.FTC <- diff(mat = mat, samples = samples, group = group, patient = patient, method = method, norm_method = norm_method)
write.table(res.FTC,"./output/result/paired/miR_diff_FTC_EVvsCF_pair.txt",quote = F,sep = "\t",row.names = T,col.names = T)

res.FTC <- res.FTC[order(res.FTC$pvalue,decreasing = F),]
res.FTC.sig <- res.FTC[res.FTC$padj<=0.1,] # abs(res.FTC$log2FoldChange)>=2 & 
res.FTC.sig <- res.FTC.sig[order(-res.FTC.sig$pvalue,abs(res.FTC.sig$log2FoldChange),decreasing = T),]
table(res.FTC.sig$log2FoldChange>0)

# tmp <- read.table("./output/lulab/diff_FTC.txt",sep = "\t",header = T,row.names = 1)
# tmp <- tmp[order(tmp$pvalue,decreasing = F),]
# tmp.sig <- tmp[tmp$padj<=0.1,]

VennDiagram::venn.diagram(x = list(rownames(res.FTC.unpair)[1:100],rownames(res.FTC)[1:100]),category.names = c("FTC_unpair","FTC_pair"), 
                          #fill=RColorBrewer::brewer.pal(4, "Set3"),
                          force.unique = T,print.mode = "raw",
                          filename = "./venn_FTC_pair_unpair.png",imagetype="png",height = 4000, width = 4000, 
                          cex = 2.5,
                          cat.cex = 2
                          #cat.pos = c(-10,+10,0),cat.dist = c(0.09,0.02)
)

#EV domain
{
  ## read sample table
  sample.table <- read.table("data/FTC_sample_table.txt", header = TRUE, check.names=FALSE, sep='\t') # row.names=1, 
  sample.table$group <- ifelse(grepl("FTC",sample.table$patient_id),"FTC","FTA")
  sample.table <- sample.table[,c("ev_small","cf_small","group","patient_id")]
  sample.table <- reshape2::melt(sample.table,id.var=c("patient_id","group"))
  colnames(sample.table)[3:4] <- c("source","data_id")
  rownames(sample.table) <- sample.table$data_id
  

  
  
  
  ## read mat
  mat.cf.domain <- read.table("data/cf_small/small_domain.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
  mat.ev.domain <- read.table("data/EV_small/small_domain.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
  #colnames(mat.ev.small) <- paste0("es",colnames(mat.ev.small) )
  table(rownames(mat.ev.domain) == rownames(mat.cf.domain))
  mat <- cbind(mat.cf.domain,mat.ev.domain)
  colnames(mat) <- unlist(lapply(strsplit(colnames(mat),"_",fixed=T),function(x) x[1]))
  
  s <- intersect(sample.table$data_id,colnames(mat))
  mat <- mat[,s]
  sample.table <- sample.table[s,]
  table(duplicated(sample.table$patient_id))
  
  
  sample.table.pair <- sample.table[!(sample.table$patient_id %in% c("FTA-14","FTA-19")),]
  positive_samples <- sample.table.pair[sample.table.pair$group=="FTC" & sample.table.pair$source=="ev_small","data_id"]
  negative_samples <- sample.table.pair[sample.table.pair$group=="FTC" & sample.table.pair$source=="cf_small","data_id"]
  samples <- c(positive_samples, negative_samples)
  sample.table.pair <- sample.table.pair[samples,]
  mat <- mat[,samples]
  group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
  patient <- sample.table.pair$patient_id[sample.table.pair$group=="FTC" & (sample.table.pair$source=="ev_small" | sample.table.pair$source=="cf_small")]
  table(patient)
  
  res.FTC <- diff(mat = mat, samples = samples, group = group, patient = patient, method = method, norm_method = norm_method)
  write.table(res.FTC,"./output/result/paired/smallRNA_diff_FTC_EVvsCF_pair.txt",quote = F,sep = "\t",row.names = T,col.names = T)
  
  res.FTC <- res.FTC[order(res.FTC$pvalue,decreasing = F),]
  res.FTC.sig <- res.FTC[res.FTC$padj<=0.1,] # abs(res.FTC$log2FoldChange)>=2 & 
  res.FTC.sig <- res.FTC.sig[order(-res.FTC.sig$pvalue,abs(res.FTC.sig$log2FoldChange),decreasing = T),]
  table(res.FTC.sig$log2FoldChange>0)
  
}

#cf domain
{
  ## read sample table
  sample.table <- read.table("data/FTC_sample_table.txt", header = TRUE, check.names=FALSE, sep='\t') # row.names=1, 
  sample.table$group <- ifelse(grepl("FTC",sample.table$patient_id),"FTC","FTA")
  sample.table <- sample.table[,c("ev_small","cf_small","group","patient_id")]
  sample.table <- reshape2::melt(sample.table,id.var=c("patient_id","group"))
  colnames(sample.table)[3:4] <- c("source","data_id")
  rownames(sample.table) <- sample.table$data_id
  ## read mat
  mat.cf.domain <- read.table("data/cf_small/small_domain.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
  mat.ev.domain <- read.table("data/EV_small/small_domain.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
  #colnames(mat.ev.small) <- paste0("es",colnames(mat.ev.small) )
  table(rownames(mat.ev.domain) == rownames(mat.cf.domain))
  mat <- cbind(mat.cf.domain,mat.ev.domain)
  colnames(mat) <- unlist(lapply(strsplit(colnames(mat),"_",fixed=T),function(x) x[1]))
  
  s <- intersect(sample.table$data_id,colnames(mat))
  mat <- mat[,s]
  sample.table <- sample.table[s,]
  table(duplicated(sample.table$patient_id))
  
  ## pair mode 
  sample.table.pair <- sample.table[!(sample.table$patient_id %in% c("FTA-14","FTA-19")),]
  positive_samples <- sample.table.pair[sample.table.pair$group=="FTA" & sample.table.pair$source=="ev_small","data_id"]
  negative_samples <- sample.table.pair[sample.table.pair$group=="FTA" & sample.table.pair$source=="cf_small","data_id"]
  samples <- c(positive_samples, negative_samples)
  sample.table.pair <- sample.table.pair[samples,]
  mat <- mat[,samples]
  group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
  patient <- sample.table.pair$patient_id[sample.table.pair$group=="FTA" & (sample.table.pair$source=="ev_small" | sample.table.pair$source=="cf_small")]
  table(patient)
  
  res.FTA <- diff(mat = mat, samples = samples, group = group, patient = patient, method = method, norm_method = norm_method)
  write.table(res.FTA,"./output/result/paired/miR_diff_FTA_EVvsCF_pair.txt",quote = F,sep = "\t",row.names = T,col.names = T)
  
  res.FTA <- res.FTA[order(res.FTA$pvalue,decreasing = F),]
  res.FTA.sig <- res.FTA[res.FTA$padj<=0.1,] # abs(res.FTA$log2FoldChange)>=2 & 
  res.FTA.sig <- res.FTA.sig[order(-res.FTA.sig$pvalue,abs(res.FTA.sig$log2FoldChange),decreasing = T),]
  table(res.FTA.sig$log2FoldChange>0)
  
}

# FTA (pair and unpair)
## read sample table
sample.table <- read.table("data/FTC_sample_table.txt", header = TRUE, check.names=FALSE, sep='\t') # row.names=1, 
sample.table$group <- ifelse(grepl("FTC",sample.table$patient_id),"FTC","FTA")
sample.table <- sample.table[,c("ev_small","cf_small","group","patient_id")]
sample.table <- reshape2::melt(sample.table,id.var=c("patient_id","group"))
colnames(sample.table)[3:4] <- c("source","data_id")
rownames(sample.table) <- sample.table$data_id

## read mat
#matrix <- "./data/lulab/FTC_EV_small/count_matrix/miRNA_count_matrix.txt"
#matrix <- "./data/lulab/FTC_EV_long/count_matrix/gencode.txt"
mat.cf.small <- read.table("data/cf_small/miRNA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
mat.ev.small <- read.table("data/EV_small/miRNA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
#colnames(mat.ev.small) <- paste0("es",colnames(mat.ev.small) )
table(rownames(mat.ev.small) == rownames(mat.cf.small))
mat.ev.small <- mat.ev.small[rownames(mat.cf.small),]
mat <- cbind(mat.cf.small,mat.ev.small)
s <- intersect(sample.table$data_id,colnames(mat))
mat <- mat[,s]
sample.table <- sample.table[s,]

## unpair mode 
positive_samples <- sample.table[sample.table$group=="FTA" & sample.table$source=="ev_small","data_id"]
negative_samples <- sample.table[sample.table$group=="FTA" & sample.table$source=="cf_small","data_id"]
samples <- c(positive_samples, negative_samples)
sample.table <- sample.table[samples,]
mat <- mat[,samples]
group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
patient <- sample.table$patient_id[sample.table$group=="FTA" & (sample.table$source=="ev_small" | sample.table$source=="cf_small")]
method <- "edger_glmlrt"
norm_method <- "TMM"
#res1 <- res
#table(patient)
res.FTA.unpair <- diff.v2(mat = mat, samples = samples, group = group, method = method, norm_method = norm_method)  #
write.table(res.FTA.unpair,"./output/result/unpaired/miR_diff_FTA_EVvsCF_unpair.txt",quote = F,sep = "\t",row.names = T,col.names = T)

res.FTA.unpair <- res.FTA.unpair[order(res.FTA.unpair$pvalue,decreasing = F),]
res.FTA.unpair.sig <- res.FTA.unpair[res.FTA.unpair$padj<=0.1,] # abs(res.FTA.unpair$log2FoldChange)>=2 & 
res.FTA.unpair.sig <- res.FTA.unpair.sig[order(-res.FTA.unpair.sig$pvalue,abs(res.FTA.unpair.sig$log2FoldChange),decreasing = T),]
table(res.FTA.unpair.sig$log2FoldChange>0)


## pair mode 
sample.table.pair <- sample.table[!(sample.table$patient_id %in% c("FTA-14","FTA-19")),]
positive_samples <- sample.table.pair[sample.table.pair$group=="FTA" & sample.table.pair$source=="ev_small","data_id"]
negative_samples <- sample.table.pair[sample.table.pair$group=="FTA" & sample.table.pair$source=="cf_small","data_id"]
samples <- c(positive_samples, negative_samples)
sample.table.pair <- sample.table.pair[samples,]
mat <- mat[,samples]
group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
patient <- sample.table.pair$patient_id[sample.table.pair$group=="FTA" & (sample.table.pair$source=="ev_small" | sample.table.pair$source=="cf_small")]
table(patient)

res.FTA <- diff(mat = mat, samples = samples, group = group, patient = patient, method = method, norm_method = norm_method)
write.table(res.FTA,"./output/result/paired/miR_diff_FTA_EVvsCF_pair.txt",quote = F,sep = "\t",row.names = T,col.names = T)

res.FTA <- res.FTA[order(res.FTA$pvalue,decreasing = F),]
res.FTA.sig <- res.FTA[res.FTA$padj<=0.1,] # abs(res.FTA$log2FoldChange)>=2 & 
res.FTA.sig <- res.FTA.sig[order(-res.FTA.sig$pvalue,abs(res.FTA.sig$log2FoldChange),decreasing = T),]
table(res.FTA.sig$log2FoldChange>0)


# FTC + FTA (combine mode)
#建议在模型中不要多组比对，即不要建立4个group因子，cf比EV的变动要大，可能会降低EV比较时的敏感性
#详见https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#multi-factor-designs
#（If I have multiple groups, should I run all together or split into pairs of groups?）
#unpaired
#FTC-12
#FTA-5
#FTA-14
#FTA-19
## read sample table
# sample.table <- read.table("meta/lulab/FTC_sample_table.txt", header = TRUE, check.names=FALSE, sep='\t') # row.names=1,
# sample.table$group <- ifelse(grepl("FTC",sample.table$patient_id),"FTC","FTA")   # cannot change
# sample.table <- sample.table[,c("ev_small","cf_small","group","patient_id")]
# sample.table <- reshape2::melt(sample.table,id.var=c("patient_id","group"))
# colnames(sample.table)[3:4] <- c("source","data_id")
# rownames(sample.table) <- sample.table$data_id
# sample.table$group <-  paste0(sample.table$group,"_",sample.table$source) ## multi comparison
# 
# ## read mat
# #matrix <- "./data/lulab/FTC_EV_small/count_matrix/miRNA_count_matrix.txt"
# #matrix <- "./data/lulab/FTC_EV_long/count_matrix/gencode.txt"
# mat.cf.small <- read.table("./data/lulab/FTC_cf_small/count_matrix/miRNA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
# mat.ev.small <- read.table("./data/lulab/FTC_EV_small/count_matrix/miRNA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
# #colnames(mat.ev.small) <- paste0("es",colnames(mat.ev.small) )
# table(rownames(mat.ev.small) == rownames(mat.cf.small))
# mat.ev.small <- mat.ev.small[rownames(mat.cf.small),]
# mat <- cbind(mat.cf.small,mat.ev.small)
# s <- intersect(sample.table$data_id,colnames(mat))
# mat <- mat[,s]
# sample.table <- sample.table[s,]
# table(sample.table$group)
# 
# ## filter
# #positive_samples <- sample.table[,"data_id"] # sample.table$group=="FTA_ev_small"
# #negative_samples <- sample.table[,"data_id"]  # sample.table$group=="FTA_cf_small"
# samples <- sample.table$data_id
# mat <- mat[,samples]
# 
# group <- sample.table$group
# patient <- sample.table$patient_id
# method <- "edger_glmlrt"
# norm_method <- "TMM"
# #res1 <- res
# 
# res <- diff(mat = mat, samples = samples, group = group, patient = patient, method = method, norm_method = norm_method)
# res <- diff.v2(mat = mat, samples = samples, group = group, method = method, norm_method = norm_method)
# #err: patientFTC-8, not full rank
# res.sig <- res[res$padj<=0.1,] # abs(res.FTC$log2FoldChange)>=2 & 
# res.sig <- res.sig[order(-res.sig$pvalue,abs(res.sig$log2FoldChange),decreasing = T),]
# table(res.sig$log2FoldChange>0)


## mds/PCA plot
sample.table <- read.table("meta/lulab/FTC_sample_table.txt", header = TRUE, check.names=FALSE, sep='\t') # row.names=1, 
sample.table$group <- ifelse(grepl("FTC",sample.table$patient_id),"FTC","FTA")   # cannot change
sample.table <- sample.table[,c("ev_small","cf_small","group","patient_id")]
sample.table <- reshape2::melt(sample.table,id.var=c("patient_id","group"))
colnames(sample.table)[3:4] <- c("source","data_id")
rownames(sample.table) <- sample.table$data_id

mat.cf.small <- read.table("./data/lulab/FTC_cf_small/count_matrix/miRNA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
mat.ev.small <- read.table("./data/lulab/FTC_EV_small/count_matrix/miRNA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
table(rownames(mat.ev.small) == rownames(mat.cf.small))
g <- intersect(rownames(mat.ev.small),rownames(mat.cf.small))
mat.ev.small <- mat.ev.small[g,]
mat.cf.small <- mat.cf.small[g,]
mat <- cbind(mat.cf.small,mat.ev.small)
#mat.cf.small <- read.table("./output/lulab/FTC_cf_small/count_matrix/filtered_miRNA_RPM.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
#mat.ev.small <- read.table("./output/lulab/FTC_EV_small/count_matrix/filtered_miRNA_RPM.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')

s <- intersect(sample.table$data_id,colnames(mat))
sample.table <- sample.table[s,]
mat <- mat[,sample.table$data_id]
#hist(as.matrix(mat),xlim = c(1,10),breaks = 200000)
y <- DGEList(counts=mat)
keep <- rowSums(edgeR::getCounts(y)>=1) >= 0.5*length(samples)
y <- y[keep,,keep.lib.sizes=FALSE]
print(dim(y))
y <- calcNormFactors(y, method="TMM")
logcpm <- edgeR::cpm(y, log=TRUE, normalized.lib.sizes = T,prior.count = 1)

# opar <- par(no.readonly = T)
# par(opar)
# par(mfrow=c(2,2))
# dim.plot = c(1,2)
# plotMDS(x = logcpm, #top = 1000, 
#         dim.plot = dim.plot, ndim = max(dim.plot), gene.selection = "pairwise", 
#         col=as.numeric(as.factor(sample.table$group)), pch=as.numeric(as.factor(sample.table$source)),
#         plot = TRUE)
#dev.off()

# PCAtools
# p <- PCAtools::pca(logcpm, metadata = sample.table, removeVar = 0.1)
# PCAtools::screeplot(p, axisLabSize = 18, titleLabSize = 22)
# PCAtools::biplot(p, colby = "group", shape = "source", legendPosition = "right", lab = "",pointSize = 5, legendIconSize = 12, legendLabSize = 30,labSize = 0) # , lab = NULL


pca <- stats::prcomp(t(logcpm), scale=TRUE)  ###prcomp函数的横行必须是样本，所以倒置一下
#choose top2 PC
pca.var <- pca$sdev^2  ## sdev是标准偏差，十个样本，就有十个标准偏差，平方是避免负数的干扰
pca.var.per <- round(pca.var/sum(pca.var)*ncol(logcpm), 1)  ##求每个样本的variation
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")  ##用柱状图可视化
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])

library(ggplot2)
pca.plot <- 
  ggplot(data=pca.data,aes(x=X,y=Y,color=sample.table$group,shape=sample.table$source))+
  geom_point(size=6) + 
  xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep="")) + 
  ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep="")) + 
  theme_bw()+
  ggsci::scale_color_d3(name="Disease")+
  #scale_color_discrete()+
  scale_shape_discrete(name="Source")+
  theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
      axis.title = element_text(size = 28,color ="black",face="bold"), 
      axis.text = element_text(size= 24,color = "black"), #,face="bold
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
      axis.text.y = element_text( hjust = 0 ), # angle = 45,
      panel.grid=element_blank(),
      legend.position = c(0.5,0.3),
      legend.text = element_text(size= 26,color = "black",face="bold"),
      legend.title= element_text(size= 26,color = "black",face="bold"))

ggsave(filename = "./pca.pdf",plot = pca.plot,)




## venn
FTA_EV <- rownames(res.FTA.sig)[res.FTA.sig$log2FoldChange>0]
FTA_cf <- rownames(res.FTA.sig)[res.FTA.sig$log2FoldChange<0]
FTC_EV <- rownames(res.FTC.sig)[res.FTC.sig$log2FoldChange>0]
FTC_cf <- rownames(res.FTC.sig)[res.FTC.sig$log2FoldChange<0]
#FTC:410
#FTC_EV:117
#FTC_cf: 69
#FTA:428
#FTA_EV:82
#FTA_cf: 62

VennDiagram::venn.diagram(x = list(FTA_EV,FTA_cf,FTC_EV,FTC_cf),category.names = c("FTA_EV","FTA_cf","FTC_EV","FTC_cf"), 
                          fill=RColorBrewer::brewer.pal(4, "Set3"),
                          force.unique = T,print.mode = "raw",
                          filename = "./venn_pair.png",imagetype="png",height = 4000, width = 4000, 
                          cex = 2.5,
                          cat.cex = 2
                          #cat.pos = c(-10,+10,0),cat.dist = c(0.09,0.02)
                          )




## heatmap
suppressPackageStartupMessages(library(pheatmap))
table(res.FTA.sig$log2FoldChange>=2 & res.FTA.sig$padj<=0.1)  # bug, cannot emit
# FTA_EV <- rownames(res.FTA.sig)[res.FTA.sig$log2FoldChange>=2 & res.FTA.sig$padj<=0.1]
# FTA_cf <- rownames(res.FTA.sig)[res.FTA.sig$log2FoldChange<=-2 & res.FTA.sig$padj<=0.1]
# table(res.FTC.sig$log2FoldChange>=2 & res.FTC.sig$padj<=0.1)  # bug, cannot emit
# FTC_EV <- rownames(res.FTC.sig)[res.FTC.sig$log2FoldChange>=2 & res.FTC.sig$padj<=0.1]
# FTC_cf <- rownames(res.FTC.sig)[res.FTC.sig$log2FoldChange<=-2 & res.FTC.sig$padj<=0.1]
# bg1 <- rownames(res.FTA)[res.FTA$padj>0.1]
# bg2 <- rownames(res.FTC)[res.FTC$padj>0.1]
FTA_EV <- rownames(res.FTA.sig)[res.FTA.sig$log2FoldChange>0 & res.FTA.sig$padj<=0.1]
FTA_cf <- rownames(res.FTA.sig)[res.FTA.sig$log2FoldChange<0 & res.FTA.sig$padj<=0.1]
FTC_EV <- rownames(res.FTC.sig)[res.FTC.sig$log2FoldChange>0 & res.FTC.sig$padj<=0.1]
FTC_cf <- rownames(res.FTC.sig)[res.FTC.sig$log2FoldChange<0 & res.FTC.sig$padj<=0.1]
bg1 <- rownames(res.FTA)[res.FTA$padj>0.1]
bg2 <- rownames(res.FTC)[res.FTC$padj>0.1]
u <- unique(rownames(res.FTA),rownames(res.FTC))
bg11 <- u[!(u %in% unique(c(FTA_EV,FTC_EV)))]  # union might be better than intersect
bg22 <- u[!(u %in% unique(c(FTA_cf,FTC_cf)))]

write.table(FTA_EV,"./output/lulab/EV_FTA.txt",sep = "\t",row.names = F,col.names = F,quote = F)
write.table(FTC_EV,"./output/lulab/EV_FTC.txt",sep = "\t",row.names = F,col.names = F,quote = F)
write.table(intersect(FTA_EV,FTC_EV),"./output/lulab/EV_intersect.txt",sep = "\t",row.names = F,col.names = F,quote = F)
write.table(FTA_cf,"./output/lulab/cf_FTA.txt",sep = "\t",row.names = F,col.names = F,quote = F)
write.table(FTC_cf,"./output/lulab/cf_FTC.txt",sep = "\t",row.names = F,col.names = F,quote = F)
write.table(intersect(FTA_cf,FTC_cf),"./output/lulab/cf_intersect.txt",sep = "\t",row.names = F,col.names = F,quote = F)
write.table(intersect(bg1,bg2),"./output/lulab/bg_intersect.txt",sep = "\t",row.names = F,col.names = F,quote = F) # not enriched
write.table(bg11,"./output/lulab/bg_EV.txt",sep = "\t",row.names = F,col.names = F,quote = F)
write.table(bg22,"./output/lulab/bg_cf.txt",sep = "\t",row.names = F,col.names = F,quote = F)

uni.ev <- intersect(FTA_EV,FTC_EV)
uni.cf <- intersect(FTA_cf,FTC_cf)
FTAonly.ev <- FTA_EV[!(FTA_EV %in% intersect(FTA_EV,FTC_EV))]
FTConly.ev <- FTC_EV[!(FTC_EV %in% intersect(FTA_EV,FTC_EV))]
FTAonly.cf <- FTA_cf[!(FTA_cf %in% intersect(FTA_cf,FTC_cf))]
FTConly.cf <- FTC_cf[!(FTC_cf %in% intersect(FTA_cf,FTC_cf))]
write.table(FTAonly.ev,"./output/lulab/EV_FTA_only.txt",sep = "\t",row.names = F,col.names = F,quote = F)
write.table(FTConly.ev,"./output/lulab/EV_FTC_only.txt",sep = "\t",row.names = F,col.names = F,quote = F)
write.table(FTAonly.cf,"./output/lulab/cf_FTA_only.txt",sep = "\t",row.names = F,col.names = F,quote = F)
write.table(FTConly.cf,"./output/lulab/cf_FTC_only.txt",sep = "\t",row.names = F,col.names = F,quote = F)


FTA.ev <- sample.table$data_id[sample.table$group=="FTA" & sample.table$source=="ev_small"]
FTA.cf <- sample.table$data_id[sample.table$group=="FTA" & sample.table$source=="cf_small"]
FTC.ev <- sample.table$data_id[sample.table$group=="FTC" & sample.table$source=="ev_small"]
FTC.cf <- sample.table$data_id[sample.table$group=="FTC" & sample.table$source=="cf_small"]

vsd.filter <- logcpm[c(uni.ev,FTAonly.ev,FTConly.ev,FTAonly.cf,FTConly.cf,uni.cf),c(FTA.ev,FTA.cf,FTC.ev,FTC.cf)]
dim(vsd.filter)
#all(rownames(res) == rownames(vsd))
#is.na(rld.filter) %>% table()
vsd.filter[1:4,1:4]

annotation_col <- sample.table[,!(colnames(sample.table) %in% c("patient_id","data_id"))]
annotation_row <- data.frame( type= c(rep("uni.ev",length(uni.ev)),
                                           rep("FTAonly.ev",length(FTAonly.ev)),
                                           rep("FTConly.ev",length(FTConly.ev)),
                                           rep("FTAonly.cf",length(FTAonly.cf)),
                                           rep("FTConly.cf",length(FTConly.cf)),
                                           rep("uni.cf",length(uni.cf))
                                           ),
                              row.names=c(uni.ev,
                                    FTAonly.ev,
                                    FTConly.ev,
                                    FTAonly.cf,
                                    FTConly.cf,
                                    uni.cf
                                    ))
annotation_row$type <- factor(annotation_row$type,levels = c("uni.ev","FTAonly.ev","FTConly.ev","FTAonly.cf","FTConly.cf","uni.cf"))
# annotation_col <- as.data.frame(colnames(vsd.filter))
# annotation_col$sampletype <- "NC"
# annotation_col$sampletype[grep(pattern = "CRC",x = annotation_col[,1])] <- "CRC"
# rownames(annotation_col) <- gsub("-me|-wgs","",annotation_col[,1])
# annotation_col <- annotation_col[-1]
# ann_colors = list(
#   sampletype = c(CRC="orange2",NC="steelblue")
# )

pheatmap(mat = vsd.filter,
         annotation_col = annotation_col, #data frame that specifies the annotations shown on left side of the heatmap.
         annotation_row = annotation_row,
         #annotation_colors = c(RColorBrewer::brewer.pal(6, "Set3")),
         scale = "row",
         #labels_col = 3, labels_row = 6,
         cluster_cols = F,cluster_rows = F,
         #cutree_cols = 2,cutree_rows = 3,
         show_colnames=T, show_rownames=F, 
         #fontsize = 8, 
         height = 9,width = 10,
         colorRampPalette(c("blue","white","red"))(100),  # steelblue/dodgerblue4-firebrick2
         #fontsize_row = 5,
         filename ="./heatmap_paired.pdf")   #filename = paste0("./output/lulab/",input.type,"_heatmap.pdf")




# GC content
library(seqinr)
convGC <- function(seqs){
  #print(seqs)
  #as.character(mySequences$`mmu-miR-15b-5p`)
  #seqs <- "cuauacaaucuacugucuuucc"
  seqs <- gsub("U|u","t",seqs)  # important, seqinr::GC seems neglect U|u base
  tmp <- seqinr::GC(seqinr::s2c(seqs),forceToLower = TRUE, exact = FALSE, NA.GC = NA, 
                    oldGC = T # not default
                    #alphabet = s2c("acgtswmkryvhdb")
  )
  #cuauacaaucuacugucuuucc	22	0.363636364
  return(tmp)
}

getGC <- function(x){
  #x <- "CL3-4"
  print(x)
  mySequences <- Biostrings::readRNAStringSet(paste0("./myNat/",x,".fa"))
  #miR <- names(mySequences)
  res <- lapply(as.character(mySequences),convGC)
  #names(res) <- names(mySequences)
  res <- do.call("rbind",res)
  res2 <- as.data.frame(cbind(GC=res[,1],mir=names(mySequences),type=x))
  return(res2)
}

## plot FTC GC box
#dat <- lapply(c("EV_intersect","cf_intersect","EV_FTC_only","cf_FTC_only","bg_intersect"),getGC) 
dat <- lapply(c("EV_FTC_only","cf_FTC_only","bg_intersect"),getGC) 
#bg_cf,bg_EV,bg_intersect,cf_FTA_only,cf_FTC_only,EV_FTA_only,EV_FTC_only,EV_FTA,EV_FTC,cf_FTA,cf_FTC,cf_intersect,EV_intersect
dat <- do.call("rbind",dat)
#dat$type <- factor(dat$type,levels = c("EV_intersect","EV_FTC_only","bg_intersect","cf_FTC_only","cf_intersect"))
dat$type <- factor(dat$type,levels = c("EV_FTC_only","bg_intersect","cf_FTC_only"))
dat$GC <- as.numeric(dat$GC)
table(dat$type)
table(is.na(dat$GC))

b <- runif(nrow(dat), -0.1, 0.1)
ggplot(dat,aes(x=type,y=GC,group=type,fill=type))+
  #geom_bar(stat = "identity",position = "dodge", color="grey30") + # 
  geom_boxplot(fill="white",outlier.alpha = 0)+
  geom_point(aes(x=as.numeric(type)+b,y=GC,color=type,fill=type),size=3)+ #,shape=type
  ylim(c(0,1.67)) +
  #ggsci::scale_fill_jama() +
  #scale_color_brewer(palette = "Set3") +
  #scale_x_discrete(labels=RNA_label)+
  #scale_fill_manual(values=RNA_color,labels=RNA_label)+
  #geom_text(aes(x=type,y=roc+0.017,label=format(roc,digits=3)),size=6,angle=0,color="grey30")+
  labs(x="",y="miRNA %CG",title = "") + 
  #scale_x_discrete(label=c("sEV-enriched 5 cell types","sEV-enriched 3-4 cell types","all unchanged","Cell-enriched 3-4 cell types","Cell-enriched 5 cell types"))+
  ggsci::scale_fill_aaas()+
  ggsci::scale_color_aaas()+
  
  #scale_fill_manual(values = c("firebrick","seagreen","black","purple3","steelblue"))+
  #scale_color_manual(values = c("firebrick","seagreen","black","purple3","steelblue"))+
  #scale_shape_manual(values = c(1,2,6,5,0))+
  ggpubr::stat_compare_means(#ref.group = "notenrich",
    comparisons = list(c("cf_FTC_only","EV_FTC_only"),
                       #c("EV_intersect","bg_intersect"),
                       c("EV_FTC_only","bg_intersect"),
                       #c("cf_intersect","bg_intersect"),
                       c("cf_FTC_only","bg_intersect")
                       ),
    label.x.npc = 0.2,hide.ns=F,size =10,step.increase=0.2,
#paired = TRUE,  
    #aes(group = Group,label = p.format), # p.signif,p.format
    label = "p.signif",
    method = "wilcox.test",
    #symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns"))
  ) +
  theme_classic() +
  theme(plot.title = element_text(size=32,color="black",hjust = 0.5),
        axis.title = element_text(size=24,color="black"), 
        axis.text = element_text(size=24,color="black"), # ,face="bold"
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 0.9, vjust = 0.5),
        panel.grid=element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold"),
        strip.background = element_blank() ) +
  #panel.spacing = unit(1.5, "lines"),
  #panel.border = element_rect(colour = "black", fill=NA, size=2)
  guides(color=guide_legend("Type"))
ggsave("./result/GC_FTC_mine.pdf",width = 12,height = 10)


## plot FTA GC box
dat <- lapply(c("EV_FTA_only","cf_FTA_only","bg_intersect"),getGC) 
#bg_cf,bg_EV,bg_intersect,cf_FTA_only,cf_FTC_only,EV_FTA_only,EV_FTC_only,EV_FTA,EV_FTC,cf_FTA,cf_FTC,cf_intersect,EV_intersect
dat <- do.call("rbind",dat)
dat$type <- factor(dat$type,levels = c("EV_FTA_only","bg_intersect","cf_FTA_only"))
dat$GC <- as.numeric(dat$GC)

b <- runif(nrow(dat), -0.1, 0.1)
ggplot(dat,aes(x=type,y=GC,group=type,fill=type))+
  #geom_bar(stat = "identity",position = "dodge", color="grey30") + # 
  geom_boxplot(fill="white",outlier.alpha = 0)+
  geom_point(aes(x=as.numeric(type)+b,y=GC,color=type,fill=type),size=3)+ #,shape=type
  ylim(c(0,1.67)) +
  #ggsci::scale_fill_jama() +
  #scale_color_brewer(palette = "Set3") +
  #scale_x_discrete(labels=RNA_label)+
  #scale_fill_manual(values=RNA_color,labels=RNA_label)+
  #geom_text(aes(x=type,y=roc+0.017,label=format(roc,digits=3)),size=6,angle=0,color="grey30")+
  labs(x="",y="miRNA %CG",title = "") + 
  #scale_x_discrete(label=c("sEV-enriched 5 cell types","sEV-enriched 3-4 cell types","all unchanged","Cell-enriched 3-4 cell types","Cell-enriched 5 cell types"))+
  ggsci::scale_fill_aaas()+
  ggsci::scale_color_aaas()+
  
  #scale_fill_manual(values = c("firebrick","seagreen","black","purple3","steelblue"))+
  #scale_color_manual(values = c("firebrick","seagreen","black","purple3","steelblue"))+
  #scale_shape_manual(values = c(1,2,6,5,0))+
  ggpubr::stat_compare_means(#ref.group = "notenrich",
    comparisons = list(
                       c("cf_FTA_only","EV_FTA_only"),
                       #c("EV_intersect","bg_intersect"),
                       c("EV_FTA_only","bg_intersect"),
                       #c("cf_intersect","bg_intersect"),
                       c("cf_FTA_only","bg_intersect")
                       
    ),
    label.x.npc = 0.2,hide.ns=F,size =10,step.increase=0.2,
    #paired = TRUE,  
    #aes(group = Group,label = p.format), # p.signif,p.format
    label = "p.signif",
    method = "wilcox.test",
    #symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns"))
  ) +
  theme_classic() +
  theme(plot.title = element_text(size=32,color="black",hjust = 0.5),
        axis.title = element_text(size=24,color="black"), 
        axis.text = element_text(size=24,color="black"), # ,face="bold"
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 0.9, vjust = 0.5),
        panel.grid=element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold"),
        strip.background = element_blank() ) +
  #panel.spacing = unit(1.5, "lines"),
  #panel.border = element_rect(colour = "black", fill=NA, size=2)
  guides(color=guide_legend("Type"))
ggsave("result/GC_FTA_mine.pdf",width = 12,height = 10)


# MFE:RNAfold
conv <- function(x){
  #x <- "EV3-4"
  print(x)
  #rnafold <- read.delim2(paste0("./output/Nat2022/",x,".txt.preMIR.RNAfold"),header = F,sep = " ") # CL3-4, CL5, EV5, EV3-4, notenrich
  rnafold <- read.delim2(paste0("./myNat/",x,".RNAfold"),header = F,sep = "") # CL3-4, CL5, EV5, EV3-4, notenrich
  mir <- rnafold$V1[grepl(">",rnafold$V1)]
  mir <- gsub(">","",mir)
  mfe <- rnafold$V3[rnafold$V3!=""]
  #  mfe <- gsub("(","",mfe,fixed = T)
  #  mfe <- gsub(")","",mfe,fixed = T)
  mfe <- gsub("(","",mfe,fixed = T)
  mfe <- gsub(")","",mfe,fixed = T)
  mfe <- as.numeric(mfe)
  rnafold.out <- as.data.frame(cbind(mir=mir,mfe=mfe))
  rnafold.out$mfe <- as.numeric(rnafold.out$mfe)
  #colnames(rnafold.out)[2] <- x
  rnafold.out$type <- x
  #max(rnafold.out$mfe)
  #rnafold.out <- rnafold.out$mfe-max(rnafold.out$mfe)
  #hist(rnafold.out$mfe) # ,xlim = c(-40,0)
  #dat[[x]] <- rnafold.out
  return(rnafold.out)
}


## FTC
#dat <- lapply(c("EV_intersect","cf_intersect","EV_FTC_only","cf_FTC_only","bg_intersect"),conv) 
dat <- lapply(c("EV_FTC_only","cf_FTC_only","bg_intersect"),conv) 

#bg_cf,bg_EV,bg_intersect,cf_FTC_only,cf_FTC_only,EV_FTC_only,EV_FTC_only,EV_FTC,EV_FTC,cf_FTC,cf_FTC,cf_intersect,EV_intersect
dat <- do.call("rbind",dat)
dat$type <- factor(dat$type,levels = c("EV_FTC_only","bg_intersect","cf_FTC_only"))
dat$mfe <- as.numeric(dat$mfe)
unique(dat$type)

b <- runif(nrow(dat), -0.1, 0.1)
ggplot(dat,aes(x=type,y=mfe,group=type,fill=type))+
  #geom_bar(stat = "identity",position = "dodge", color="grey30") + # 
  geom_boxplot(fill="white",outlier.alpha = 0)+
  geom_point(aes(x=as.numeric(type)+b,y=mfe,color=type,fill=type),size=3)+ #,shape=type
  ylim(c(-10,6)) +
  #ggsci::scale_fill_jama() +
  #scale_color_brewer(palette = "Set3") +
  #scale_x_discrete(labels=RNA_label)+
  #scale_fill_manual(values=RNA_color,labels=RNA_label)+
  #geom_text(aes(x=type,y=roc+0.017,label=format(roc,digits=3)),size=6,angle=0,color="grey30")+
  labs(x="",y="delta G (kcal/mol)",title = "") + 
  #scale_x_discrete(label=c("sEV-enriched 5 cell types","sEV-enriched 3-4 cell types","all unchanged","Cell-enriched 3-4 cell types","Cell-enriched 5 cell types"))+
  ggsci::scale_fill_aaas()+
  ggsci::scale_color_aaas()+
  
  #scale_fill_manual(values = c("firebrick","seagreen","black","purple3","steelblue"))+
  #scale_color_manual(values = c("firebrick","seagreen","black","purple3","steelblue"))+
  #scale_shape_manual(values = c(1,2,6,5,0))+
  ggpubr::stat_compare_means(#ref.group = "notenrich",
    comparisons = list(
      c("cf_FTC_only","EV_FTC_only"),
      #c("EV_intersect","bg_intersect"),
      c("EV_FTC_only","bg_intersect"),
      #c("cf_intersect","bg_intersect"),
      c("cf_FTC_only","bg_intersect")
      
    ),
    label.x.npc = 0.2,hide.ns=F,size =10,step.increase=0.2,
    #paired = TRUE,  
    #aes(group = Group,label = p.format), # p.signif,p.format
    label = "p.signif",
    method = "wilcox.test",
    #symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns"))
  ) +
  theme_classic() +
  theme(plot.title = element_text(size=32,color="black",hjust = 0.5),
        axis.title = element_text(size=24,color="black"), 
        axis.text = element_text(size=24,color="black"), # ,face="bold"
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 0.9, vjust = 0.5),
        panel.grid=element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold"),
        strip.background = element_blank() ) +
  #panel.spacing = unit(1.5, "lines"),
  #panel.border = element_rect(colour = "black", fill=NA, size=2)
  guides(color=guide_legend("Type"))
ggsave("./result/MFE_FTC_mine.pdf",width = 12,height = 10)

## FTA
dat <- lapply(c("EV_FTA_only","cf_FTA_only","bg_intersect"),conv) 
#bg_cf,bg_EV,bg_intersect,cf_FTA_only,cf_FTC_only,EV_FTA_only,EV_FTC_only,EV_FTA,EV_FTC,cf_FTA,cf_FTC,cf_intersect,EV_intersect
dat <- do.call("rbind",dat)
dat$type <- factor(dat$type,levels = c("EV_FTA_only","bg_intersect","cf_FTA_only"))
dat$mfe <- as.numeric(dat$mfe)
unique(dat$type)

b <- runif(nrow(dat), -0.1, 0.1)
ggplot(dat,aes(x=type,y=mfe,group=type,fill=type))+
  #geom_bar(stat = "identity",position = "dodge", color="grey30") + # 
  geom_boxplot(fill="white",outlier.alpha = 0)+
  geom_point(aes(x=as.numeric(type)+b,y=mfe,color=type,fill=type),size=3)+ #,shape=type
  ylim(c(-10,6)) +
  #ggsci::scale_fill_jama() +
  #scale_color_brewer(palette = "Set3") +
  #scale_x_discrete(labels=RNA_label)+
  #scale_fill_manual(values=RNA_color,labels=RNA_label)+
  #geom_text(aes(x=type,y=roc+0.017,label=format(roc,digits=3)),size=6,angle=0,color="grey30")+
  labs(x="",y="delta G (kcal/mol)",title = "") + 
  #scale_x_discrete(label=c("sEV-enriched 5 cell types","sEV-enriched 3-4 cell types","all unchanged","Cell-enriched 3-4 cell types","Cell-enriched 5 cell types"))+
  ggsci::scale_fill_aaas()+
  ggsci::scale_color_aaas()+
  
  #scale_fill_manual(values = c("firebrick","seagreen","black","purple3","steelblue"))+
  #scale_color_manual(values = c("firebrick","seagreen","black","purple3","steelblue"))+
  #scale_shape_manual(values = c(1,2,6,5,0))+
  ggpubr::stat_compare_means(#ref.group = "notenrich",
    comparisons = list(
      c("cf_FTA_only","EV_FTA_only"),
      #c("EV_intersect","bg_intersect"),
      c("EV_FTA_only","bg_intersect"),
      #c("cf_intersect","bg_intersect"),
      c("cf_FTA_only","bg_intersect")
      
    ),
    label.x.npc = 0.2,hide.ns=F,size =10,step.increase=0.2,
    #paired = TRUE,  
    #aes(group = Group,label = p.format), # p.signif,p.format
    label = "p.signif",
    method = "wilcox.test",
    #symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01,1),symbols = c( "***", "**", "*", "ns"))
  ) +
  theme_classic() +
  theme(plot.title = element_text(size=32,color="black",hjust = 0.5),
        axis.title = element_text(size=24,color="black"), 
        axis.text = element_text(size=24,color="black"), # ,face="bold"
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 0.9, vjust = 0.5),
        panel.grid=element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size= 24),
        legend.title= element_text(size= 24,face="bold"),
        strip.text = element_text(size = 20,face="bold"),
        strip.background = element_blank() ) +
  #panel.spacing = unit(1.5, "lines"),
  #panel.border = element_rect(colour = "black", fill=NA, size=2)
  guides(color=guide_legend("Type"))
ggsave("./result/MFE_FTA_mine.pdf",width = 12,height = 10)





# diff WHK-WHCSU sRNA FTCvsFTA --------------------------------------------------------------------
#EV domain
{
## read sample table
sample.table <- read.table("data/FTC_sample_table.txt", header = TRUE, check.names=FALSE, sep='\t') # row.names=1, 
sample.table$group <- ifelse(grepl("FTC",sample.table$patient_id),"FTC","FTA")
sample.table <- sample.table[,c("ev_small","cf_small","group","patient_id")]
sample.table <- reshape2::melt(sample.table,id.var=c("patient_id","group"))
colnames(sample.table)[3:4] <- c("source","data_id")
rownames(sample.table) <- sample.table$data_id
## read mat
mat.cf.domain <- read.table("data/cf_small/small_domain.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
mat.ev.domain <- read.table("data/EV_small/small_domain.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
#colnames(mat.ev.small) <- paste0("es",colnames(mat.ev.small) )
table(rownames(mat.ev.domain) == rownames(mat.cf.domain))
mat <- cbind(mat.cf.domain,mat.ev.domain)
colnames(mat) <- unlist(lapply(strsplit(colnames(mat),"_",fixed=T),function(x) x[1]))

s <- intersect(sample.table$data_id,colnames(mat))
mat <- mat[,s]
sample.table <- sample.table[s,]
table(duplicated(sample.table$patient_id))

## filter
positive_samples <- sample.table[sample.table$group=="FTC" & sample.table$source=="ev_small","data_id"]
negative_samples <- sample.table[sample.table$group=="FTA" & sample.table$source=="ev_small","data_id"]
samples <- c(positive_samples, negative_samples)
mat <- mat[,samples]
group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
method <- "edger_glmlrt"
norm_method <- "TMM"

res.ev.domain <- diff.v2(mat = mat, samples = samples, group = group,  method = method, norm_method = norm_method) # patient = patient,
write.table(res.ev.domain,"output/result/unpaired/miR_diff_EV_FTCvsFTA_domain.txt",quote = F,sep = "\t",row.names = T,col.names = T)

## summary
# res.ev.sig <- res.ev[res.ev$padj<=0.1,] 
res.ev.domain.sig <- res.ev.domain[res.ev.domain$pvalue < 0.01 & abs(res.ev.domain$log2FoldChange) > 1,]
res.ev.domain.sig <- res.ev.domain.sig[order(-res.ev.domain.sig$pvalue,abs(res.ev.domain.sig$log2FoldChange),decreasing = T),]
table(res.ev.domain.sig$log2FoldChange>0)


}
#cf domain
{
  ## read sample table
  sample.table <- read.table("data/FTC_sample_table.txt", header = TRUE, check.names=FALSE, sep='\t') # row.names=1, 
  sample.table$group <- ifelse(grepl("FTC",sample.table$patient_id),"FTC","FTA")
  sample.table <- sample.table[,c("ev_small","cf_small","group","patient_id")]
  sample.table <- reshape2::melt(sample.table,id.var=c("patient_id","group"))
  colnames(sample.table)[3:4] <- c("source","data_id")
  rownames(sample.table) <- sample.table$data_id
  ## read mat
  mat.cf.domain <- read.table("data/cf_small/small_domain.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
  mat.ev.domain <- read.table("data/EV_small/small_domain.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
  #colnames(mat.ev.small) <- paste0("es",colnames(mat.ev.small) )
  table(rownames(mat.ev.domain) == rownames(mat.cf.domain))
  mat <- cbind(mat.cf.domain,mat.ev.domain)
  colnames(mat) <- unlist(lapply(strsplit(colnames(mat),"_",fixed=T),function(x) x[1]))
  
  s <- intersect(sample.table$data_id,colnames(mat))
  mat <- mat[,s]
  sample.table <- sample.table[s,]
  table(duplicated(sample.table$patient_id))
  
  ## filter
  positive_samples <- sample.table[sample.table$group=="FTC" & sample.table$source=="cf_small","data_id"]
  negative_samples <- sample.table[sample.table$group=="FTA" & sample.table$source=="cf_small","data_id"]
  samples <- c(positive_samples, negative_samples)
  mat <- mat[,samples]
  group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
  method <- "edger_glmlrt"
  norm_method <- "TMM"
  
  res.cf.domain <- diff.v2(mat = mat, samples = samples, group = group,  method = method, norm_method = norm_method) # patient = patient,
  write.table(res.cf.domain,"output/result/unpaired/miR_diff_cf_FTCvsFTA_domain.txt",quote = F,sep = "\t",row.names = T,col.names = T)
  
  ## summary
  # res.ev.sig <- res.ev[res.ev$padj<=0.1,] 
  res.cf.domain.sig <- res.cf.domain[res.cf.domain$pvalue < 0.01 & abs(res.cf.domain$log2FoldChange) > 1,]
  res.cf.domain.sig <- res.cf.domain.sig[order(-res.cf.domain.sig$pvalue,abs(res.cf.domain.sig$log2FoldChange),decreasing = T),]
  table(res.cf.domain.sig$log2FoldChange>0)
}

# EV (not paired mode)
## read sample table
sample.table <- read.table("data/FTC_sample_table.txt", header = TRUE, check.names=FALSE, sep='\t') # row.names=1, 
sample.table$group <- ifelse(grepl("FTC",sample.table$patient_id),"FTC","FTA")
sample.table <- sample.table[,c("ev_small","cf_small","group","patient_id")]
sample.table <- reshape2::melt(sample.table,id.var=c("patient_id","group"))
colnames(sample.table)[3:4] <- c("source","data_id")
rownames(sample.table) <- sample.table$data_id

## read mat
mat.cf.small <- read.table("data/cf_small/miRNA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
mat.ev.small <- read.table("data/EV_small/miRNA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
#colnames(mat.ev.small) <- paste0("es",colnames(mat.ev.small) )
table(rownames(mat.ev.small) == rownames(mat.cf.small))
mat.ev.small <- mat.ev.small[rownames(mat.cf.small),]
mat <- cbind(mat.cf.small,mat.ev.small)
s <- intersect(sample.table$data_id,colnames(mat))
mat <- mat[,s]
sample.table <- sample.table[s,]
table(duplicated(sample.table$patient_id))

## filter
positive_samples <- sample.table[sample.table$group=="FTC" & sample.table$source=="ev_small","data_id"]
negative_samples <- sample.table[sample.table$group=="FTA" & sample.table$source=="ev_small","data_id"]
samples <- c(positive_samples, negative_samples)
mat <- mat[,samples]
group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
#patient <- sample.table$patient_id[sample.table$source=="ev_small" & (sample.table$group=="FTA" | sample.table$group=="FTC")]  # !!! need match()
method <- "edger_glmlrt"
norm_method <- "TMM"
#res1 <- res

res.ev <- diff.v2(mat = mat, samples = samples, group = group,  method = method, norm_method = norm_method) # patient = patient,
write.table(res.ev,"output/result/unpaired/miR_diff_EV_FTCvsFTA_unpair.txt",quote = F,sep = "\t",row.names = T,col.names = T)

## summary
# res.ev.sig <- res.ev[res.ev$padj<=0.1,] 
res.ev.sig <- res.ev[res.ev$pvalue < 0.01 & abs(res.ev$log2FoldChange) > 1,]
res.ev.sig <- res.ev.sig[order(-res.ev.sig$pvalue,abs(res.ev.sig$log2FoldChange),decreasing = T),]
table(res.ev.sig$log2FoldChange>0)


# cf (not paired mode)
## read sample table
sample.table <- read.table("data/FTC_sample_table.txt", header = TRUE, check.names=FALSE, sep='\t') # row.names=1, 
sample.table$group <- ifelse(grepl("FTC",sample.table$patient_id),"FTC","FTA")
sample.table <- sample.table[,c("ev_small","cf_small","group","patient_id")]
sample.table <- reshape2::melt(sample.table,id.var=c("patient_id","group"))
colnames(sample.table)[3:4] <- c("source","data_id")
rownames(sample.table) <- sample.table$data_id

## read mat
#matrix <- "./data/lulab/FTC_EV_small/count_matrix/miRNA_count_matrix.txt"
#matrix <- "./data/lulab/FTC_EV_long/count_matrix/gencode.txt"
mat.cf.small <- read.table("data/cf_small/miRNA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
mat.ev.small <- read.table("data/EV_small/miRNA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
#colnames(mat.ev.small) <- paste0("es",colnames(mat.ev.small) )
table(rownames(mat.ev.small) == rownames(mat.cf.small))
mat.ev.small <- mat.ev.small[rownames(mat.cf.small),]
mat <- cbind(mat.cf.small,mat.ev.small)
s <- intersect(sample.table$data_id,colnames(mat))
mat <- mat[,s]
sample.table <- sample.table[s,]
table(duplicated(sample.table$patient_id))

## filter
positive_samples <- sample.table[sample.table$group=="FTC" & sample.table$source=="cf_small","data_id"]
negative_samples <- sample.table[sample.table$group=="FTA" & sample.table$source=="cf_small","data_id"]
samples <- c(positive_samples, negative_samples)
mat <- mat[,samples]
group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
#patient <- sample.table$patient_id[sample.table$source=="ev_small" & (sample.table$group=="FTA" | sample.table$group=="FTC")]  # !!! need match()
method <- "edger_glmlrt"
norm_method <- "TMM"
#res1 <- res

res.cf <- diff.v2(mat = mat, samples = samples, group = group,  method = method, norm_method = norm_method) # patient = patient,
write.table(res.cf,"output/result/unpaired/miR_diff_cf_FTCvsFTA_unpair.txt",quote = F,sep = "\t",row.names = T,col.names = T)

## summary
# res.cf.sig <- res.cf[res.cf$padj<=0.1,] # abs(res.cf$log2FoldChange)>=2 & 
# res.cf.sig <- res.cf.sig[order(-res.cf.sig$pvalue,abs(res.cf.sig$log2FoldChange),decreasing = T),]
# table(res.cf.sig$log2FoldChange>0)

res.cf.sig <- res.cf[res.cf$pvalue < 0.01 & abs(res.cf$log2FoldChange) > 1,]
res.cf.sig <- res.cf.sig[order(-res.cf.sig$pvalue,abs(res.cf.sig$log2FoldChange),decreasing = T),]
table(res.cf.sig$log2FoldChange>0)
#[1] 348  17
#  TRUE 
#   2

### volcano plot
{
  res.cf.plt <- res.cf
  res.cf.plt <- res.cf.plt[order(-res.cf.plt$pvalue,abs(res.cf.plt$log2FoldChange),decreasing = T),]
  res.cf.plt$threshold <- factor(ifelse(res.cf.plt$pvalue < 0.01 & abs(res.cf.plt$log2FoldChange) >=1, ifelse(res.cf.plt$log2FoldChange > 1 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
  
  res.ev.plt <- res.ev
  res.ev.plt <- res.ev.plt[order(-res.ev.plt$pvalue,abs(res.ev.plt$log2FoldChange),decreasing = T),]
  res.ev.plt$threshold <- factor(ifelse(res.ev.plt$pvalue < 0.01 & abs(res.ev.plt$log2FoldChange) >=1, ifelse(res.ev.plt$log2FoldChange > 1 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
  
  res.cf.domain.plt <- res.cf.domain
  res.cf.domain.plt <- res.cf.domain.plt[order(-res.cf.domain.plt$pvalue,abs(res.cf.domain.plt$log2FoldChange),decreasing = T),]
  res.cf.domain.plt$threshold <- factor(ifelse(res.cf.domain.plt$pvalue < 0.01 & abs(res.cf.domain.plt$log2FoldChange) >=1, ifelse(res.cf.domain.plt$log2FoldChange > 1 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
  
  res.ev.domain.plt <- res.ev.domain
  res.ev.domain.plt <- res.ev.domain.plt[order(-res.ev.domain.plt$pvalue,abs(res.ev.domain.plt$log2FoldChange),decreasing = T),]
  res.ev.domain.plt$threshold <- factor(ifelse(res.ev.domain.plt$pvalue < 0.01 & abs(res.ev.domain.plt$log2FoldChange) >=1, ifelse(res.ev.domain.plt$log2FoldChange > 1 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
  
  res.FTC.unpair.plt <- res.FTC.unpair 
  res.FTC.unpair.plt <- res.FTC.unpair.plt[order(-res.FTC.unpair.plt$pvalue,abs(res.FTC.unpair.plt$log2FoldChange),decreasing = T),]
  res.FTC.unpair.plt$threshold <- factor(ifelse(res.FTC.unpair.plt$padj < 0.1, ifelse(res.FTC.unpair.plt$log2FoldChange > 0 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
  
  res.FTA.unpair.plt <- res.FTA.unpair 
  res.FTA.unpair.plt <- res.FTA.unpair.plt[order(-res.FTA.unpair.plt$pvalue,abs(res.FTA.unpair.plt$log2FoldChange),decreasing = T),]
  res.FTA.unpair.plt$threshold <- factor(ifelse(res.FTA.unpair.plt$padj < 0.1, ifelse(res.FTA.unpair.plt$log2FoldChange > 0 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
  
  
  res.FTC.plt <- res.FTC
  res.FTC.plt <- res.FTC.plt[order(-res.FTC.plt$pvalue,abs(res.FTC.plt$log2FoldChange),decreasing = T),]
  res.FTC.plt$threshold <- factor(ifelse(res.FTC.plt$padj < 0.1, ifelse(res.FTC.plt$log2FoldChange > 0 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
  
  res.FTA.plt <- res.FTA
  res.FTA.plt <- res.FTA.plt[order(-res.FTA.plt$pvalue,abs(res.FTA.plt$log2FoldChange),decreasing = T),]
  res.FTA.plt$threshold <- factor(ifelse(res.FTA.plt$padj < 0.1, ifelse(res.FTA.plt$log2FoldChange > 0 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
  
  col = brewer.pal(3,"Set1")
  {
    p=ggplot(res.FTA.plt, aes(x=log2FoldChange, y =-log10(padj), color=threshold)) +  #y =-log10(pvalue)
      scale_color_manual(values=c(col[2],col[1],"grey"))+
      geom_point(size=2, alpha=0.8, shape=16, stroke = 0)+
      #scale_y_continuous(limits=c(0,2.8))+
      theme_bw(base_size = 12) +
      geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+
      geom_hline(yintercept = -log10(0.1),lty=4,col="grey",lwd=0.6)+
      theme(legend.position="right",
            panel.grid=element_blank(),
            panel.border=element_rect(size=1, fill="transparent"),
            legend.title = element_blank(),
            legend.text= element_text(face="bold", color="black", size=15),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(face="bold", color="black", size=15),
            axis.text.y = element_text(face="bold",  color="black", size=15),
            axis.title.x = element_text(face="bold", color="black", size=15),
            axis.title.y = element_text(face="bold",color="black", size=15))+
      labs(x="log2FoldChange",y="-log10 (FDR)",title="FTA", face="bold")
    p
  }
  ggsave("output/result/paired/diffexp_FTA_pair_volcano.pdf",p,width = 5,height = 3)
}

### venn
{
  color = brewer.pal(6, "Set3")
  venn <- venn.diagram(
    x = list(      
      # EV=rownames(res.ev.plt[res.ev.plt$threshold == "Up",]),
      # cf=rownames(res.cf.plt[res.cf.plt$threshold == "Up",])
      # EV=rownames(res.ev.plt[res.ev.plt$threshold == "Down",]),
      # cf=rownames(res.cf.plt[res.cf.plt$threshold == "Down",])
      # FTC=rownames(res.FTC.unpair.plt[res.FTC.unpair.plt$threshold == "Up",]),
      # FTA=rownames(res.FTA.unpair.plt[res.FTA.unpair.plt$threshold == "Up",])
      # FTC=rownames(res.FTC.unpair.plt[res.FTC.unpair.plt$threshold == "Down",]),
      # FTA=rownames(res.FTA.unpair.plt[res.FTA.unpair.plt$threshold == "Down",])
      # 
      FTC=rownames(res.FTC.plt[res.FTC.plt$threshold == "Up",]),
      FTA=rownames(res.FTA.plt[res.FTA.plt$threshold == "Up",])
      # FTC=rownames(res.FTC.plt[res.FTC.plt$threshold == "Down",]),
      # FTA=rownames(res.FTA.plt[res.FTA.plt$threshold == "Down",])

    ),
    filename = NULL,
    fill = color[1:2],
    alpha = 0.7,
    lwd = 2,
    lty = "blank",
    col = "transparent", #线条色
    cex = 2,
    fontface = "bold", #加粗
    fontfamily = "sans",
    cat.col = "black",
    cat.default.pos = "text", #位置, outer内 text外
    cat.cex = 1.5,
    cat.fontfamily = "sans",
    cat.fontface = "bold",     
    cat.dist = c(0.15, 0.15), 
    cat.pos = c(0,0),
    main = "EV-enriched",
    main.col = "black",
    main.cex = 1.5,
    main.fontfamily = "sans",
    main.fontface = "bold"
  );
  grid.draw(venn);
  ggsave("output/result/paired/EV_venn.pdf",venn,width = 5,height = 5)
}



### heatmap####
sample.table <- read.table("data/FTC_sample_table.txt", header = TRUE, check.names=FALSE, sep='\t') # row.names=1, 
sample.table$group <- ifelse(grepl("FTC",sample.table$patient_id),"FTC","FTA")
sample.table <- sample.table[,c("ev_small","cf_small","group","patient_id")]
sample.table <- reshape2::melt(sample.table,id.var=c("patient_id","group"))
colnames(sample.table)[3:4] <- c("source","data_id")
rownames(sample.table) <- sample.table$data_id

mat.cf.small <- read.table("data/cf_small/miRNA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
mat.ev.small <- read.table("data/EV_small/miRNA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')

mat.cf.domain <- read.table("data/cf_small/small_domain.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
mat.ev.domain <- read.table("data/EV_small/small_domain.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
colnames(mat.ev.domain) <- unlist(lapply(strsplit(colnames(mat.ev.domain),"_1",fixed=TRUE), function(x) x[1]))
mat.ev.domain <- mat.ev.domain[,colnames(mat.ev.small)]
#colnames(mat.ev.small) <- paste0("es",colnames(mat.ev.small) )
table(rownames(mat.ev.small) == rownames(mat.cf.small))


# EV
positive_samples <- sample.table[sample.table$group=="FTC" & sample.table$source=="ev_small","data_id"]
negative_samples <- sample.table[sample.table$group=="FTA" & sample.table$source=="ev_small","data_id"]

# cf
positive_samples <- sample.table[sample.table$group=="FTC" & sample.table$source=="cf_small","data_id"]
negative_samples <- sample.table[sample.table$group=="FTA" & sample.table$source=="cf_small","data_id"]
negative_samples <- setdiff(negative_samples,c('csFTA-14','csFTA-19'))



# FTC 
positive_samples <- sample.table[sample.table$group=="FTC" & sample.table$source=="ev_small","data_id"]
negative_samples <- sample.table[sample.table$group=="FTC" & sample.table$source=="cf_small","data_id"]

# FTA unpaired
positive_samples <- sample.table[sample.table$group=="FTA" & sample.table$source=="ev_small","data_id"]
negative_samples <- sample.table[sample.table$group=="FTA" & sample.table$source=="cf_small","data_id"]
negative_samples <- setdiff(negative_samples,c('csFTA-14','csFTA-19'))

# FTA paired
positive_samples <- sample.table[sample.table$group=="FTA" & sample.table$source=="ev_small","data_id"]
negative_samples <- sample.table[sample.table$group=="FTA" & sample.table$source=="cf_small","data_id"]
negative_samples <- setdiff(negative_samples,c('csFTA-14','csFTA-19'))
positive_samples <- setdiff(positive_samples,c('FTA-14','FTA-19'))

samples <- c(positive_samples, negative_samples)
sample.table <- sample.table[match(samples,sample.table$data_id),]

mat.ev.small <- mat.ev.small[rownames(mat.cf.small),]
tmp <- cbind(mat.ev.small,mat.cf.small)
colnames(tmp) <- unlist(lapply(strsplit(colnames(tmp),"_",fixed=T),function(x) x[1]))
mat <- cbind(mat.ev.small,mat.cf.small)
colnames(mat) <- unlist(lapply(strsplit(colnames(mat),"_",fixed=T),function(x) x[1]))
mat <- mat[,colnames(tmp)]

logRPM <- log2(cpm(rbind(mat,tmp)) + 1)
s <- intersect(sample.table$data_id,colnames(mat))
mat <- mat[,s]
sample.table <- sample.table[s,]

{
  #FTCvsFTA domain
  difgenes <- c(rownames(res.ev.domain.plt[res.ev.domain.plt$threshold == "Up",]),rownames(res.ev.domain.plt[res.ev.domain.plt$threshold == "Down",]))
  difgenes <- c(rownames(res.cf.domain.plt[res.cf.domain.plt$threshold == "Up",]),rownames(res.cf.domain.plt[res.cf.domain.plt$threshold == "Down",]))
  
  
  #FTCvsFTA
  difgenes <- c(rownames(res.ev.plt[res.ev.plt$threshold == "Up",]),rownames(res.ev.plt[res.ev.plt$threshold == "Down",]))
  difgenes <- c(rownames(res.cf.plt[res.cf.plt$threshold == "Up",]),rownames(res.cf.plt[res.cf.plt$threshold == "Down",]))
  
  #EVvscf unpair
  difgenes <- c(rownames(res.FTC.unpair.plt[res.FTC.unpair.plt$threshold == "Up",]),rownames(res.FTC.unpair.plt[res.FTC.unpair.plt$threshold == "Down",]))
  difgenes <- c(rownames(res.FTA.unpair.plt[res.FTA.unpair.plt$threshold == "Up",]),rownames(res.FTA.unpair.plt[res.FTA.unpair.plt$threshold == "Down",]))

  #EVvscf pair
  difgenes <- c(rownames(res.FTC.plt[res.FTC.plt$threshold == "Up",]),rownames(res.FTC.plt[res.FTC.plt$threshold == "Down",]))
  difgenes <- c(rownames(res.FTA.plt[res.FTA.plt$threshold == "Up",]),rownames(res.FTA.plt[res.FTA.plt$threshold == "Down",]))
  

  logRPM <- log2(cpm(mat) + 0.1)
  logRPM <- logRPM[geneset$V1[1:10],]
  rownames(logRPM) <- gene.
  
  tmp <- res.ev.sig[geneset$V1[1:10],]
  rownames(tmp) <- gene.
  tmp['hsa-miR-4306',] <- res.ev.sig['hsa-miR-4306',]
  tmp['hsa-miR-6721-5p',] <- res.ev.sig['hsa-miR-6721-5p',]

  logRPM.scale <- scale(t(logRPM[order(tmp$log2FoldChange),]), center = T, scale = T)
  logRPM.scale <- t(logRPM.scale)
  tmp[order(tmp$log2FoldChange),]
  
  write.table(logRPM.scale,"output/result/unpaired/heatmap_EV_miRNA.txt",quote = F,sep = "\t",row.names = T,col.names = T)
  write.table(tmp[order(tmp$log2FoldChange),],"output/result/unpaired/top10_ev_smallRNA.txt",quote = F,sep = "\t",row.names = T,col.names = T)
  
  
  tbs.order <- read.table("data/tbs_order.txt", header = TRUE, check.names=F,sep='\t')
  row.names(tbs.order) <- tbs.order[,1]
  logRPM.scale[,as.factor(tbs.order['miRNA',2:22])]
  logRPM.scale <- cbind(logRPM.scale[,as.factor(tbs.order['miRNA',2:22])],logRPM.scale[,paste('cs',as.factor(tbs.order['miRNA',2:22]),sep='')])
  
  class. <- sample.table$source[order(sample.table$source,sample.table$group)]
  cancer. <- sample.table$group[order(sample.table$source,sample.table$group)]
  tbs. <- as.factor(c(rep("",21),colnames(tbs.order)[2:22]))
  
  class. <- sample.table$source[order(sample.table$source,sample.table$group)]
  cancer. <- sample.table$group[order(sample.table$source,sample.table$group)]
  # 
  #tbs.order <- tbs.order[,-1]
  # 
  # logRPM.scale <- scale(t(logRPM), center = T, scale = T)
  # logRPM.scale <- t(logRPM.scale)[,order(sample.table$group,sample.table$source)]
  # class. <- sample.table$source[order(sample.table$group,sample.table$source)]
  # cancer. <- sample.table$group[order(sample.table$group,sample.table$source)]
  
  ## pheatmap::pheatmap 
  {
    library(pheatmap)
    ann_col = data.frame(#cancer. = as.character(cancer.),
                         class. = as.character(class.),
                         tbs. = as.character(tbs.))
    #ann_col = data.frame(class. = as.character(class.))
    rownames(ann_col) = colnames(logRPM.scale)
    #ann_row = data.frame(GeneClass = factor(rep(c("CRC.cf","LUAD.cf","NC.EV", "CRC.EV", "NC.CRC.EV","LUAD.EV"), c(5,148,8,0,1,5))))
    #rownames(ann_row) = difgenes
    #set colors of each group
    ann_colors = list( class. = brewer.pal(2,"Set1")[1:2], tbs. = c(brewer.pal(5,"YlOrRd")[1:5],'white'))
    #ann_colors = list(class. = brewer.pal(3,"Set1")[1:2]) 
    names(ann_colors$class.) = c("cf_small","ev_small")
    #names(ann_colors$cancer.) = c("FTC","FTA")  cancer. = brewer.pal(3,"Set1")[1:2],
    names(ann_colors$tbs.) = c("primary","M-liver","M-bone","M-bone-lung","M-lung","")
    
    col = brewer.pal(3,"Set1")
    pdf("output/result/unpaired/pheatmap_FTC_TBS.pdf", width = 6, height = 7)
    par(mar=c(2,2,2,2))
    pheatmap(logRPM.scale, 
             #color = colorRampPalette(c(col[2],"white",col[1]))(1000),
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
             cutree_col = 2, 
             #cutree_row = 3, #break up the heatmap by clusters you define
             cluster_rows = F, 
             cluster_cols = F, #by default, pheatmap clusters by both row and col
             show_rownames = F,
             show_colnames = F,
             fontsize_col = 12,
             angle_col = 0,
             annotation_col = ann_col,
             border_color = "NA",
             #annotation_row = ann_row,
             annotation_colors = ann_colors,
             annotation_names_row = F)
    dev.off()
  }
  ## pheatmap::pheatmap 
  {
    library(pheatmap)
    # ann_col = data.frame(cancer. = as.character(cancer.),
    #                      class. = as.character(class.))
    ann_col = data.frame(cancer = as.character(cancer.))
    rownames(ann_col) = colnames(logRPM.scale)
    #ann_row = data.frame(GeneClass = factor(rep(c("CRC.cf","LUAD.cf","NC.EV", "CRC.EV", "NC.CRC.EV","LUAD.EV"), c(5,148,8,0,1,5))))
    #rownames(ann_row) = difgenes
    #set colors of each group
    #ann_colors = list(cancer. = brewer.pal(3,"Set1")[1:2], class. = brewer.pal(5,"Set1")[4:5])
    ann_colors = list(cancer. = brewer.pal(3,"Set1")[1:2])
    #names(ann_colors$class.) = c("cf_small","ev_small")
    names(ann_colors$cancer.) = c("FTA","FTC")
    
    col = brewer.pal(3,"Set1")
    pdf("output/result/unpaired/pheatmap_ev_small_top10.pdf", width = 7, height = 2.4)
    par(mar=c(2,2,2,2))
    pheatmap(logRPM.scale, 
             #color = colorRampPalette(c(col[2],"white",col[1]))(1000),
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
             cutree_col = 2, 
             #cutree_row = 3, #break up the heatmap by clusters you define
             cluster_rows = F, 
             cluster_cols = F, #by default, pheatmap clusters by both row and col
             show_rownames = T,
             show_colnames = F,
             fontsize_col = 12,
             angle_col = 0,
             annotation_col = ann_col,
             border_color = "NA",
             #annotation_row = ann_row,
             main ='top 10 miRNA',
             annotation_colors = ann_colors,
             annotation_names_row = F)
    dev.off()
  }
}


{
  library(pheatmap)
  library(grid)
  library(gtable)
  
  # Modified pheatmap:::heatmap_motor
  heatmap_motor <- function (matrix, border_color, cellwidth, cellheight, tree_col, 
                             tree_row, treeheight_col, treeheight_row, filename, width, 
                             height, breaks, color, legend, annotation_row, annotation_col, 
                             annotation_colors, annotation_legend, annotation_names_row, 
                             annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
                             hjust_col, vjust_col, angle_col, fmat, fontsize_number, number_color, 
                             gaps_col, gaps_row, labels_row, labels_col, ...) 
  {
    lo = pheatmap:::lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), 
                       ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, 
                       treeheight_col = treeheight_col, treeheight_row = treeheight_row, 
                       legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, 
                       annotation_colors = annotation_colors, annotation_legend = annotation_legend, 
                       annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, 
                       main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
                       fontsize_col = fontsize_col, angle_col = angle_col, gaps_row = gaps_row, 
                       gaps_col = gaps_col, ...)
    res = lo$gt
    mindim = lo$mindim
    if (!is.na(filename)) {
      if (is.na(height)) {
        height = convertHeight(gtable_height(res), "inches", valueOnly = T)
      }
      if (is.na(width)) {
        width = convertWidth(gtable_width(res), "inches", valueOnly = T)
      }
      r = regexpr("\\.[a-zA-Z]*$", filename)
      if (r == -1) 
        stop("Improper filename")
      ending = substr(filename, r + 1, r + attr(r, "match.length"))
      f = switch(ending, pdf = function(x, ...) pdf(x, ...), 
                 png = function(x, ...) png(x, units = "in", res = 300, 
                                            ...), jpeg = function(x, ...) jpeg(x, units = "in", 
                                                                               res = 300, ...), jpg = function(x, ...) jpeg(x, 
                                                                                                                            units = "in", res = 300, ...), tiff = function(x, 
                                                                                                                                                                           ...) tiff(x, units = "in", res = 300, compression = "lzw", 
                                                                                                                                                                                     ...), bmp = function(x, ...) bmp(x, units = "in", 
                                                                                                                                                                                                                      res = 300, ...), stop("File type should be: pdf, png, bmp, jpg, tiff"))
      f(filename, height = height, width = width)
      gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, 
                         border_color = border_color, tree_col = tree_col, 
                         tree_row = tree_row, treeheight_col = treeheight_col, 
                         treeheight_row = treeheight_row, breaks = breaks, 
                         color = color, legend = legend, annotation_col = annotation_col, 
                         annotation_row = annotation_row, annotation_colors = annotation_colors, 
                         annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, 
                         annotation_names_col = annotation_names_col, filename = NA, 
                         main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
                         fontsize_col = fontsize_col, hjust_col = hjust_col, 
                         vjust_col = vjust_col, angle_col = angle_col, fmat = fmat, 
                         fontsize_number = fontsize_number, number_color = number_color, 
                         labels_row = labels_row, labels_col = labels_col, 
                         gaps_col = gaps_col, gaps_row = gaps_row, ...)
      grid.draw(gt)
      dev.off()
      return(gt)
    }
    if (mindim < 3) 
      border_color = NA
    if (!is.na(main)) {
      elem = pheatmap:::draw_main(main, fontsize = 1.3 * fontsize, ...)
      res = gtable_add_grob(res, elem, t = 1, l = 3, name = "main", 
                            clip = "off")
    }
    if (!pheatmap:::is.na2(tree_col) & treeheight_col != 0) {
      elem = pheatmap:::draw_dendrogram(tree_col, gaps_col, horizontal = T)
      res = gtable_add_grob(res, elem, t = 2, l = 3, name = "col_tree")
    }
    if (!pheatmap:::is.na2(tree_row) & treeheight_row != 0) {
      elem = pheatmap:::draw_dendrogram(tree_row, gaps_row, horizontal = F)
      res = gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
    }
    elem = pheatmap:::draw_matrix(matrix, border_color, gaps_row, gaps_col, 
                                  fmat, fontsize_number, number_color)
    res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
                          name = "matrix")
    if (length(labels_col) != 0) {
      pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, 
                  hjust_col = hjust_col, vjust_col = vjust_col, angle_col = angle_col, 
                  ...)
      elem = do.call(pheatmap:::draw_colnames, pars)
      res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", 
                            name = "col_names")
    }
    if (length(labels_row) != 0) {
      pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row, 
                  ...)
      elem = do.call(pheatmap:::draw_rownames, pars)
      res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
                            name = "row_names")
    }
    if (!pheatmap:::is.na2(annotation_col)) {
      converted_annotation = pheatmap:::convert_annotations(annotation_col, 
                                                            annotation_colors)
      elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
                                         gaps_col, fontsize, horizontal = T)
      res = gtable_add_grob(res, elem, t = 3, l = 3, clip = "off", 
                            name = "col_annotation")
      if (annotation_names_col) {
        elem = pheatmap:::draw_annotation_names(annotation_col, fontsize, 
                                                horizontal = T)
        res = gtable_add_grob(res, elem, t = 3, l = 4, clip = "off", 
                              name = "col_annotation_names")
      }
    }
    if (!pheatmap:::is.na2(annotation_row)) {
      converted_annotation = pheatmap:::convert_annotations(annotation_row, 
                                                            annotation_colors)
      elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
                                         gaps_row, fontsize, horizontal = F)
      res = gtable_add_grob(res, elem, t = 4, l = 2, clip = "off", 
                            name = "row_annotation")
      if (annotation_names_row) {
        elem = pheatmap:::draw_annotation_names(annotation_row, fontsize, 
                                                horizontal = F, hjust_col = hjust_col, vjust_col = vjust_col, 
                                                angle_col = angle_col)
        res = gtable_add_grob(res, elem, t = 5, l = 2, clip = "off", 
                              name = "row_annotation_names")
      }
    }
    annotation = c(annotation_col[length(annotation_col):1], 
                   annotation_row[length(annotation_row):1])
    annotation = annotation[unlist(lapply(annotation, function(x) !pheatmap:::is.na2(x)))]
    if (length(annotation) > 0 & annotation_legend) {
      elem = pheatmap:::draw_annotation_legend(annotation, annotation_colors, 
                                               border_color, fontsize = fontsize, ...)
      t = ifelse(is.null(labels_row), 4, 3)
      res = gtable_add_grob(res, elem, t = t, l = 6, b = 5, 
                            clip = "off", name = "annotation_legend")
    }
    if (!pheatmap:::is.na2(legend)) {
      elem = pheatmap:::draw_legend(color, breaks, legend, fontsize = fontsize, 
                                    ...)
      t = ifelse(is.null(labels_row), 4, 3)
      res = gtable_add_grob(res, elem, t = t, l = 5, b = 5, 
                            clip = "off", name = "legend")
    }
    return(res)
  }
  
  # Modified pheatmap:::lo    
  lo <- function (rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, 
                  treeheight_col, treeheight_row, legend, annotation_row, annotation_col, 
                  annotation_colors, annotation_legend, annotation_names_row, 
                  annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
                  angle_col, gaps_row, gaps_col, ...) 
  {
    if (!is.null(coln[1]) | (!pheatmap:::is.na2(annotation_row) & annotation_names_row)) {
      if (!is.null(coln[1])) {
        t = coln
      }
      else {
        t = ""
      }
      tw = strwidth(t, units = "in", cex = fontsize_col/fontsize)
      if (annotation_names_row) {
        t = c(t, colnames(annotation_row))
        tw = c(tw, strwidth(colnames(annotation_row), units = "in"))
      }
      longest_coln = which.max(tw)
      gp = list(fontsize = ifelse(longest_coln <= length(coln), 
                                  fontsize_col, fontsize), ...)
      coln_height = unit(1, "grobheight", textGrob(t[longest_coln], 
                                                   rot = angle_col, gp = do.call(gpar, gp))) + unit(10, 
                                                                                                    "bigpts")
    }
    else {
      coln_height = unit(5, "bigpts")
    }
    if (!is.null(rown[1])) {
      t = rown
      tw = strwidth(t, units = "in", cex = fontsize_row/fontsize)
      if (annotation_names_col) {
        t = c(t, colnames(annotation_col))
        tw = c(tw, strwidth(colnames(annotation_col), units = "in"))
      }
      longest_rown = which.max(tw)
      gp = list(fontsize = ifelse(longest_rown <= length(rown), 
                                  fontsize_row, fontsize), ...)
      rown_width = unit(1, "grobwidth", textGrob(t[longest_rown], 
                                                 rot = 0, gp = do.call(gpar, gp))) + unit(10, "bigpts")
    }
    else {
      rown_width = unit(5, "bigpts")
    }
    gp = list(fontsize = fontsize, ...)
    if (!pheatmap:::is.na2(legend)) {
      longest_break = which.max(nchar(names(legend)))
      longest_break = unit(1.1, "grobwidth", 
                           textGrob(as.character(names(legend))[longest_break], 
                                    gp = do.call(gpar, gp)))
      title_length = unit(1.1, "grobwidth", textGrob("Scale", 
                                                     gp = gpar(fontface = "bold", ...)))
      legend_width = unit(12, "bigpts") + longest_break * 1.2
      legend_width = max(title_length, legend_width)
    }
    else {
      legend_width = unit(0, "bigpts")
    }
    if (is.na(main)) {
      main_height = unit(0, "npc")
    }
    else {
      main_height = unit(1.5, "grobheight", textGrob(main, 
                                                     gp = gpar(fontsize = 1.3 * fontsize, ...)))
    }
    textheight = unit(fontsize, "bigpts")
    if (!pheatmap:::is.na2(annotation_col)) {
      annot_col_height = ncol(annotation_col) * (textheight + 
                                                   unit(2, "bigpts")) + unit(2, "bigpts")
      t = c(as.vector(as.matrix(annotation_col)), colnames(annotation_col))
      annot_col_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
                                                               gp = gpar(...))) + unit(12, "bigpts")
      if (!annotation_legend) {
        annot_col_legend_width = unit(0, "npc")
      }
    }
    else {
      annot_col_height = unit(0, "bigpts")
      annot_col_legend_width = unit(0, "bigpts")
    }
    if (!pheatmap:::is.na2(annotation_row)) {
      annot_row_width = ncol(annotation_row) * (textheight + 
                                                  unit(2, "bigpts")) + unit(2, "bigpts")
      t = c(as.vector(as.matrix(annotation_row)), colnames(annotation_row))
      annot_row_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
                                                               gp = gpar(...))) + unit(12, "bigpts")
      if (!annotation_legend) {
        annot_row_legend_width = unit(0, "npc")
      }
    }
    else {
      annot_row_width = unit(0, "bigpts")
      annot_row_legend_width = unit(0, "bigpts")
    }
    annot_legend_width = max(annot_row_legend_width, annot_col_legend_width)
    treeheight_col = unit(treeheight_col, "bigpts") + unit(5, 
                                                           "bigpts")
    treeheight_row = unit(treeheight_row, "bigpts") + unit(5, 
                                                           "bigpts")
    if (is.na(cellwidth)) {
      mat_width = unit(1, "npc") - rown_width - legend_width - 
        treeheight_row - annot_row_width - annot_legend_width
    }
    else {
      mat_width = unit(cellwidth * ncol, "bigpts") + length(gaps_col) * 
        unit(4, "bigpts")
    }
    if (is.na(cellheight)) {
      mat_height = unit(1, "npc") - main_height - coln_height - 
        treeheight_col - annot_col_height
    }
    else {
      mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * 
        unit(4, "bigpts")
    }
    gt = gtable(widths = unit.c(treeheight_row, rown_width,  
                                mat_width, treeheight_row, legend_width, annot_legend_width), 
                heights = unit.c(main_height, treeheight_col, annot_col_height, 
                                 mat_height, coln_height), vp = viewport(gp = do.call(gpar, 
                                                                                      gp)))
    cw = convertWidth(mat_width - (length(gaps_col) * unit(4, 
                                                           "bigpts")), "bigpts", valueOnly = T)/ncol
    ch = convertHeight(mat_height - (length(gaps_row) * unit(4, 
                                                             "bigpts")), "bigpts", valueOnly = T)/nrow
    mindim = min(cw, ch)
    res = list(gt = gt, mindim = mindim)
    return(res)
  }
  
  # Modified pheatmap:::draw_rownames      
  draw_rownames <- function (rown, gaps, ...) 
  {
    coord = pheatmap:::find_coordinates(length(rown), gaps)
    y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
    res = textGrob(rown, x = unit(-3, "bigpts"), y = y, vjust = 0.5, 
                   hjust = 1, gp = gpar(...))
    return(res)
  }
  
  # assignInNamespace(x="draw_rownames", value=draw_rownames, ns="pheatmap")
  # assignInNamespace(x="lo", value=lo, ns="pheatmap")
  # assignInNamespace(x="heatmap_motor", value=heatmap_motor, ns="pheatmap") 
  assignInNamespace(x="draw_rownames", value=draw_rownames, ns="pheatmap")
  assignInNamespace(x="lo", value=lo, ns="pheatmap")
  assignInNamespace(x="heatmap_motor", value=heatmap_motor, ns="pheatmap")
  
}


miRNA_ID <- read.table("data/cf_small/miRNA_count_matrix.txt", header = TRUE, row.names=1, check.names=FALSE, sep='\t')
miRNA_ID <- rownames(miRNA_ID)
miRNA_ID <- unlist(lapply(strsplit(miRNA_ID,"-3p",fixed=TRUE), function(x) x[1]))
miRNA_ID <- unlist(lapply(strsplit(miRNA_ID,"-5p",fixed=TRUE), function(x) x[1]))

write.table(miRNA_ID,'data/miRNA_ID.txt',sep='\t',quote = F,row.names = F,col.names = F)
####diff top10####

write.table(rownames(res.ev.sig[order(abs(res.ev.sig$log2FoldChange),-res.ev.sig$pvalue,decreasing = T),])[1:10],"C:/Users/86185/bioinf/lulab/diff_top10/data/miRNA_diff_top10.txt",quote = F,sep = "\t",row.names = F,col.names = F)
write.table(rownames(res.ev.domain.sig[order(abs(res.ev.domain.sig$log2FoldChange),-res.ev.domain.sig$pvalue,decreasing = T),])[1:10],"C:/Users/86185/bioinf/lulab/diff_top10/data/smallRNA_diff_top10.txt",quote = F,sep = "\t",row.names = F,col.names = F)
write.table(intersect(rownames(res.FTA.sig[res.FTA.sig$log2FoldChange > 0,]),rownames(res.FTA.sig[res.FTC.sig$log2FoldChange > 0,])),
            "C:/Users/86185/bioinf/lulab/EV_enrich_intersect/data/miRNA/EV_enrich.txt",quote = F,sep = "\t",row.names = F,col.names = F)
write.table(intersect(rownames(res.FTA.sig[res.FTA.sig$log2FoldChange > 0,]),rownames(res.FTA.sig[res.FTC.sig$log2FoldChange > 0,])),
            "C:/Users/86185/bioinf/lulab/EV_enrich_intersect/data/smallRNA/EV_enrich.txt",quote = F,sep = "\t",row.names = F,col.names = F)



# test pholegenic tree ----------------------------------------------------
library(Biostrings)
library(adegenet)
library(ape)
library(ggtree)
library(ggplot2)
library(stats)
library(ips)
library(msa)
library(seqinr)

#system.file("tex", "texshade.sty", package="msa")
#mySequenceFile <- system.file("examples", "exampleAA.fasta", package="msa")
#mySequences <- readAAStringSet(mySequenceFile)
mySequences <- Biostrings::readDNAStringSet("cf_EV_intersect.fa")
#mySequences <- readDNAStringSet("output/NC2013/cell_exosome_union_T.fa")
mySequences
myFirstAlignment <- msa::msa(mySequences, type = "dna", #order = "input", 
                        method = "ClustalW",gapOpening = 12,gapExtension = 3)   # , "ClustalOmega", "Muscle"
## use default substitution matrix
myFirstAlignment
print(myFirstAlignment, show="complete")

## plot multi-align
# msaPrettyPrint(myFirstAlignment, output="pdf", showNames="none",
#                showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
# msaPrettyPrint(myFirstAlignment, output="pdf", y=c(164, 213),
#                subset=c(1:6), showNames="none", showLogo="top",
#                logoColors="rasmol", shadingMode="similar",
#                showLegend=FALSE, askForOverwrite=FALSE)
msaPrettyPrint(myFirstAlignment, output="pdf",# y=c(164, 213),
               showNames="left", shadingMode="functional",
               shadingModeArg="structure",
               askForOverwrite=FALSE)



## prepare tree
#hemoSeq <- readAAStringSet(system.file("examples/HemoglobinAA.fasta",package="msa"))
#hemoAln <- msa(hemoSeq)
#hemoAln2 <- msaConvert(hemoAln, type="seqinr::alignment")
myFirstAlignment2 <- msaConvert(myFirstAlignment, type="seqinr::alignment")
d <- dist.alignment(myFirstAlignment2, "identity")
tmp <- as.matrix(d)  # distance matrix

#hsa-miR-148b
#dist.mat <- distance(birds.t, method = "jaccard")
#true.dist.mat <- as.dist(dist.mat)
#clust.res <- hclust(true.dist.mat, method = "complete")


#d <- colnames(d)
#library(ape)
#hemoTree <- ape::nj(d)  # Neighbor-Joining Tree Estimation, d may be an object of class “dist”
hemoTree2 <- stats::hclust(d)  # method: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"
#?hclust
# plot(hemoTree, main="Phylogenetic Tree of Hemoglobin Alpha Sequences")
# plot(hemoTree2, main="Phylogenetic Tree of Hemoglobin Alpha Sequences")
# 
# plot(as.phylo(hemoTree), type = "fan", #tip.color = as.factor(1),
#      edge.color = "steelblue", edge.width = 2, edge.lty = 2,
#      label.offset = 1, cex = 0.7)
# 
# plot(ape::as.phylo(hemoTree2), type = "fan")
#?as.phylo  # converts an object into a tree of class "phylo"


## plot with ggdendro
# library("ggplot2")
# library("ggdendro")
# ggdendrogram(as.dendrogram(hemoTree2), rotate = TRUE, theme_dendro = FALSE)+
#   ggsci::scale_color_nejm()+
#   
#   theme_classic() + 
#   theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
#         axis.title = element_text(size = 28,color ="black",face="bold"), 
#         axis.text = element_text(size= 24,color = "black"), #,face="bold
#         panel.grid.minor.y = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
#         axis.text.y = element_text( hjust = 1 ), # angle = 45,
#         panel.grid=element_blank(),
#         #legend.position = "top",
#         legend.text = element_text(size= 26,color = "black",face="bold"),
#         legend.title= element_text(size= 26,color = "black",face="bold"))


## plot with dendextend
for (m in c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")){
  #m <- "ward.D"
hemoTree2 <- stats::hclust(d,method = m)  # method: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"
  

dend <- as.dendrogram(hemoTree2,horiz=F)
nname <- hemoTree2$labels

cf <- read.fas("./cf_intersect.fa")
#cf <- read.fas("./output/NC2013/cell.fa")
color <- ifelse(nname %in% names(cf),"blue","red")
table(color)

library(dendextend)
#par(mfrow = c(1,2), mar = c(5,2,1,0))
dend <- dend %>% 
  dendextend::color_labels(labels=nname,col=color) %>%  
  dendextend::color_branches(k = 2,col = c("steelblue4","firebrick")) %>%
  #set("branches_lwd", c(2,1,2)) %>%
  #set("branches_lty", c(1,2,1)) %>% 
  set("nodes_cex", c(4,3,4)) %>% 
  set("leaves_cex", c(4,3,4))



pdf(paste0("./tree-",m,".pdf"),width = 14,height = 20)
opar <- par()
par(oma=c(3,3,3,3), mar=c(10,10,10,10))

plot(dend,las = 1,horiz = T,
     cex.lab = 2,
     cex.axis = 2,
     cex.main = 2,
     cex.sub = 2,
     main = m, 
     xlab = 'Distance', 
     ylab = 'Motif') 
legend(#x = 0.8, y = 18, 
       x = 2.5, y = 38,
       #"topleft",
       title = "Cluster",  # Title
       title.adj = 0.5,         # Horizontal adjustment of the title
       title.col = "black",      # Color of the title
       legend = c("cluster1", "cluster2"),
       lty = c(1, 1), col = c("steelblue4","firebrick"), lwd = 2, cex=1)
legend(#x = 0.3, y = 18,
       x = 1.8, y = 38,
       #"top",
       #inset = c(-0.45, 0),
       title = "Source",  # Title
       title.adj = 0.5,         # Horizontal adjustment of the title
       title.col = "black",      # Color of the title
       #title.cex=2,
       #pch = 1,
       legend = c("CFmotif", "EXOmotif"),
       lty = 2, col = c("blue", "red"), lwd = 2,cex=1)
dev.off()
#legend()
}
#dend <- color_labels(dend, k = 2)
#plot(dend)








set.seed(12)
fit <- kmeans(iris[1:4],3)
unclass(fit)
attributes(fit)
sapply(fit,class)

args(limma.trend)
formals(limma.trend)


methods(numeric)
getAnywhere(numeric)

Rprof()
print(1:1000)
Rprof(NULL)
summaryRprof()

args(mad)
mad(1:10)

debug(mad)
mad(1:10)
undebug(mad)

trace(mad,edit = T)
mad(1:10)
untrace(mad)

x <- -1
if(x<0) browser()

getAnywhere(mad)


x <- 1:4
x
get("x")
assign("y",1:4)
y

system("R CMD install npar www.statmethods.net/RiA/npar_1.0.tar.gz")

intersect(intersect(rownames(res.FTC.sig[res.FTC.sig$log2FoldChange > 0,]),rownames(res.ev.sig)),intersect(rownames(res.FTA.sig[res.FTC.sig$log2FoldChange > 0,]),rownames(res.ev.sig)))
intersect(intersect(rownames(res.FTC.sig[res.FTC.sig$log2FoldChange > 0,]),rownames(res.cf.sig)),intersect(rownames(res.FTA.sig[res.FTC.sig$log2FoldChange > 0,]),rownames(res.cf.sig)))

intersect(intersect(rownames(res.FTC.sig[res.FTC.unpair.sig$log2FoldChange > 0,]),rownames(res.cf.sig)),intersect(rownames(res.FTA.unpair.sig[res.FTC.sig$log2FoldChange > 0,]),rownames(res.cf.sig)))
intersect(intersect(rownames(res.FTC.sig[res.FTC.unpair.sig$log2FoldChange > 0,]),rownames(res.ev.sig)),intersect(rownames(res.FTA.unpair.sig[res.FTC.sig$log2FoldChange > 0,]),rownames(res.ev.sig)))
