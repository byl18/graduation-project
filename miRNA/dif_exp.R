########################################
###### differential expression for i-pico
###### zhanqing
###### 20220310

.libPaths()
library(ggplot2)
library(gg.gap)
library(dplyr)
library(ggpubr)
library(ggsci)
library(tidyverse)
library(edgeR)
library(RUVSeq)
library(RColorBrewer)
library(VennDiagram)
library(GGally)
library(gplots)
library(pheatmap)
library(scales)
library(reshape2)
display.brewer.all()

setwd("C:/study/lulab_research/毕设/zq/repeatZQ/output/")

### Annotate library ID with sample ID
annotate <- function(input){
  # input is library id, output is sample_id
  library_id <- input
  metadata <- read.table("./output/count_matrix/annotation.txt", sep = "\t", header = T)
  index=match(library_id, metadata$library_id)
  sample_id = c()
  for (i in 1:length(library_id)){
    sample_id[i] = metadata$sample_id[index[i]]
  }
  input <- sample_id
}


# gencode+circRNA
circrna <- read.delim2("count_matrix/circRNA.txt", header = T, row.names = 1, check.names = F)
circrna <- circrna[,colnames(gencode)]
data <- rbind(gencode,circrna)

########### unpaired
{
  gencode <- read.delim2("./output/count_matrix/gencode.cutoff_filter.txt", header = T, row.names = 1, check.names = F)
  dim(gencode)
  gencode <- gencode[,!colnames(gencode)%in%c("H-1","H-2","H-3","H-4","H-5")] ## HCC样本不要了
  data=gencode
  ### group
  {
    # put cell-free as reference level
    # cancer
    cancer <- unlist(lapply(strsplit(colnames(data),"-",fixed=T),function(x) x[1]))
    cancer = factor(cancer, levels = c("NC","CRC","LUAD"))
    # class
    class <- unlist(lapply(strsplit(colnames(data),"-",fixed=T),function(x) x[3]))
    class[which(is.na(class))] <- "cf"
    class = factor(class, levels = c("cf","EV"))
  }
  ### edgeR: tpm cutoff to filter genes
  {
    # lib.size to be specific as lib.size(hg38_dedup)
    # remove genes that are expressed less than 50% samples(tpm)
    {
      lib.size <- read.delim2("./output/count_matrix/libsize.txt",sep="\t",header=T,row.names = 1)
      rownames(lib.size) <- annotate(rownames(lib.size))
      lib.size <- lib.size[colnames(data),]
      y <- DGEList(counts=data, lib.size=lib.size)
      #此时对lib.size进行normalized目的是通过tpm filter genes
      RPM <- round(edgeR::cpm(y, normalzied.lib.sizes = T),3)
      geneid <- rownames(RPM)
      split <- strsplit(geneid,"|",fixed=T)
      length <- unlist(lapply(split,function(split) split[2]))
      y$genes <- data.frame(length=as.numeric(length))
      RPKM <- round(rpkm(y),3)
      TPM <- round((RPKM/colSums(RPKM))*10^6,3)
      cpm <- as.data.frame(RPM)
      tpm <- as.data.frame(TPM)
      ###################
      ### group = cancer
      {
        group = cancer
        
        nc = tpm[,group=="NC"]
        crc = tpm[,group=="CRC"]
        luad = tpm[,group=="LUAD"]
        n.nc = length(group[group=="NC"])/2
        n.crc = length(group[group=="CRC"])/2
        n.luad = length(group[group=="LUAD"])/2
        
        detected_genes = data.frame(
          NC <- c(sum(rowSums(nc>2)>= n.nc),
                  sum(rowSums(nc>5)>= n.nc),
                  sum(rowSums(nc>10)>= n.nc)),
          CRC <- c(sum(rowSums(crc>2)>= n.crc),
                   sum(rowSums(crc>5)>= n.crc),
                   sum(rowSums(crc>10)>= n.crc)),
          LUAD <- c(sum(rowSums(luad>2)>= n.luad),
                    sum(rowSums(luad>5)>= n.luad),
                    sum(rowSums(luad>10)>= n.luad))
          )
        colnames(detected_genes) = c("NC","CRC","LUAD")
        rownames(detected_genes) = c("TPM>2(50%)","TPM>5(50%)","TPM>10(50%)")
        detected <- detected_genes %>% mutate(group=rownames(detected_genes))
        detected <- melt(detected)
        detected$group = factor(detected$group, levels=c("TPM>2(50%)","TPM>5(50%)","TPM>10(50%)"))
        {
          p=ggplot(detected, aes(x=variable, y=value, fill=group))+
            #geom_bar(stat="identity",width=0.7, col="black")+
            geom_bar(stat="identity",position="dodge",width=0.7, col="black")+
            #scale_y_continuous(limits=c(0,20000))+
            scale_fill_brewer()+
            theme_bw()+
            guides(fill=guide_legend(title=NULL, ncol=1))+
            theme(
              plot.margin = unit(c(1, 1, 0, 1),"cm"),
              panel.border=element_rect(size=1, fill="transparent"),
              panel.grid=element_blank(),
              #legend.position=c(0.25,0.85),
              legend.position="top",
              legend.background = element_blank(),
              legend.title = element_blank(),
              legend.text= element_text(face="bold", color="black", size=12),
              axis.text.x = element_text(face="bold",color="black", size=12),
              axis.text.y = element_text(face="bold",color="black", size=12),
              axis.title.x = element_text(face="bold",color="black", size=12),
              axis.title.y = element_text(face="bold",color="black", size=12))+
            ylab("Number of detected genes")+xlab("")
          #geom_vline(aes(xintercept=6.5))+
          p
        }   
        ggsave("./output/result/unpaired/gene_number_cancer_TPMnormLibsize.pdf",p,width=4.5,height=4.5) 
      }
      filter = (rowSums(nc>5) >= n.nc) | (rowSums(crc>5) >= n.crc) | (rowSums(luad>5) >= n.luad)
      table(filter)
      y <- y[filter,]
      # 10794
  }
}
  ### edgeR: unPaired differential expression
  {
    # not normalized with lib.size
    y <- DGEList(counts=data, group=cancer)
    y <- y[filter,]   # filter: tpm>10: 5114 用5k个基因还是1w个基因差别不大
    dim(y)
    # 14112,73
    # 10794,68
    
    # TMM normalization
    #EDASeq::plotRLE(edgeR::cpm(y)) 
    y <- calcNormFactors(y, method="TMM")
    #EDASeq::plotRLE(edgeR::cpm(y))  
    
    limma::plotMDS(y)
    limma::plotMDS(edgeR::cpm(y, log=T), labels = class)
    limma::plotMDS(edgeR::cpm(y, log=T), labels = cancer)
    ########################################################
    ###### design: unpaired test in all samples
    cancer
    class
    
    ########################################################
    ###### design: unpaired test in cell-free and EV dataset
    cancer
    class
    
    test="EV"
    cutoff="pValue"
    difexp_cancerNC <- function(test){
      data.1 <- data[,class==test]
      cancer.1 <- cancer[class==test]
      cancer.1 <- ifelse(cancer.1=="NC","NC","cancer")
      cancer.1 <- factor(cancer.1, levels=c("NC","cancer"))
      class.1 <- class[class==test]
      
      design <- model.matrix(~cancer.1)
      
      y <- DGEList(counts=data.1, group=cancer.1)
      y <- y[filter,] 
      dim(y)
      y <- calcNormFactors(y, method="TMM")
      #limma::plotMDS(edgeR::cpm(y, log=T), labels = class.1)
      
      rownames(design) <- colnames(y)
      y <- estimateDisp(y, design)
      y$common.dispersion
      plotBCV(y)
      
      fit.ql <- glmQLFit(y, design)
      qlf <- glmQLFTest(fit.ql, coef="cancer.1cancer")
      de <- topTags(qlf, n=Inf)$table
      
      de.cancer.down <- filter(de, logFC < -1, PValue < 0.01)
      de.cancer.up <- filter(de, logFC > 1, PValue < 0.01)
      
      #de.cancer.down <- filter(de, logFC < -1, FDR < 0.1)
      #de.cancer.up <- filter(de, logFC > 1, FDR < 0.1)

      print(paste0("cancer upregulated RNA: ",nrow(de.cancer.up)))
      print(paste0("cancer downregulated RNA: ",nrow(de.cancer.down)))
      return(list(de=de, 
                  de.cancer.up=de.cancer.up, 
                  de.cancer.down=de.cancer.down))
    }
    difexp_eachcancer <- function(test){
      data.1 <- data[,class==test]
      cancer.1 <- factor(cancer[class==test], levels=c("NC","CRC","LUAD"))
      data.1 <- data.1[,!is.na(cancer.1)]
      cancer.1 <- cancer.1[!is.na(cancer.1)]
      class.1 <- class[class==test]
      
      design <- model.matrix(~cancer.1)
      
      y <- DGEList(counts=data.1, group=cancer.1)
      y <- y[filter,] 
      dim(y)
      y <- calcNormFactors(y, method="TMM")
      #limma::plotMDS(edgeR::cpm(y, log=T), labels = class.1)
      
      rownames(design) <- colnames(y)
      y <- estimateDisp(y, design)
      y$common.dispersion
      plotBCV(y)
      
      fit.ql <- glmQLFit(y, design)
      qlf.crc <- glmQLFTest(fit.ql, coef="cancer.1CRC")
      de.crc <- topTags(qlf.crc, n=Inf)$table
      qlf.luad <- glmQLFTest(fit.ql, coef="cancer.1LUAD")
      de.luad <- topTags(qlf.luad, n=Inf)$table
      #qlf.hcc <- glmQLFTest(fit.ql, coef="cancer.1HCC")
      #de.hcc <- topTags(qlf.hcc, n=Inf)$table

      de.crc.up <- filter(de.crc, logFC > 1, PValue < 0.01)
      de.crc.down <- filter(de.crc, logFC < -1, PValue < 0.01)
      de.luad.up <- filter(de.luad, logFC > 1,PValue < 0.01)
      de.luad.down <- filter(de.luad, logFC < -1, PValue < 0.01)
      
      #de.crc.up <- filter(de.crc, logFC > 1, FDR < 0.1)
      #de.crc.down <- filter(de.crc, logFC < -1, FDR < 0.1)
      #de.luad.up <- filter(de.luad, logFC > 1, FDR < 0.1)
      #de.luad.down <- filter(de.luad, logFC < -1, FDR < 0.1)

      print(paste0("CRC upregulated RNA: ",nrow(de.crc.up)))
      print(paste0("CRC downregulated RNA: ",nrow(de.crc.down)))
      print(paste0("LUAD upregulated RNA: ",nrow(de.luad.up)))
      print(paste0("LUAD downregulated RNA: ",nrow(de.luad.down)))
      return(list(de.crc=de.crc,
                  de.luad=de.luad,
                  de.crc.up=de.crc.up,
                  de.crc.down=de.crc.down,
                  de.luad.up=de.luad.up,
                  de.luad.down=de.luad.down))
    }

    de.ev <- difexp_cancerNC("EV")
    de.cf <- difexp_cancerNC("cf")
    de.each.ev <- difexp_eachcancer("EV")
    de.each.cf <- difexp_eachcancer("cf")
    
  }
  ### venn
  {
    color = brewer.pal(6, "Set3")
    venn <- venn.diagram(
      x = list(      
        # EV.cancer.up=rownames(de.ev$de.cancer.up),
        # cf.cancer.up=rownames(de.cf$de.cancer.up)
        EV.cancer.down=rownames(de.ev$de.cancer.down),
        cf.cancer.down=rownames(de.cf$de.cancer.down)
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
      main = "",
      main.col = "black",
      main.cex = 1.5
      #main.fontface = "bold"
    );
    grid.draw(venn);
    #ggsave("result/unpaired/cancer_difgenes_up_venn.pdf",venn,width = 5,height = 5)
  }
  ### volcano plot
  {
    de.ev.plt <- de.ev$de
    de.ev.plt$threshold <- factor(ifelse(de.ev.plt$PValue < 0.01 & abs(de.ev.plt$logFC) >=1, ifelse(de.ev.plt$logFC > 1 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
    de.cf.plt <- de.cf$de
    de.cf.plt$threshold <- factor(ifelse(de.cf.plt$PValue < 0.01 & abs(de.cf.plt$logFC) >=1, ifelse(de.cf.plt$logFC > 1 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
   
    col = brewer.pal(3,"Set1")
    {
      p=ggplot(de.ev.plt, aes(x=logFC, y =-log10(PValue), color=threshold)) + 
        scale_color_manual(values=c(col[2],col[1],"grey"))+
        geom_point(size=2, alpha=0.8, shape=16, stroke = 0)+
        #scale_y_continuous(limits=c(0,2.8))+
        theme_bw(base_size = 12) +
        geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+
        geom_hline(yintercept = -log10(0.01),lty=4,col="grey",lwd=0.6)+
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
        labs(x="log2FoldChange",y="-log10 (p-value)",title="", face="bold")
      p
    }
    ggsave("output/result/unpaired/difexp_ev_volcano.pdf",p,width = 5,height = 3)
  }
  ### heatmap
  {
    difgenes <- c(rownames(de.ev$de.cancer.up),rownames(de.ev$de.cancer.down))
    difgenes <- c(rownames(de.cf$de.cancer.up),rownames(de.cf$de.cancer.down))
    
    logRPM <- edgeR::cpm(y, log=T)
    logRPM <- logRPM[difgenes,]
    logRPM.scale <- scale(t(logRPM), center = T, scale = T)
    logRPM.scale <- t(logRPM.scale)[,order(class,cancer)]
    class. <- class[order(class,cancer)]
    cancer. <- cancer[order(class,cancer)]
    
    ## pheatmap::pheatmap 
    {
      library(pheatmap)
      ann_col = data.frame(cancer. = as.character(cancer.),
                           class. = as.character(class.))
      rownames(ann_col) = colnames(logRPM.scale)
      #ann_row = data.frame(GeneClass = factor(rep(c("CRC.cf","LUAD.cf","NC.EV", "CRC.EV", "NC.CRC.EV","LUAD.EV"), c(5,148,8,0,1,5))))
      #rownames(ann_row) = difgenes
      #set colors of each group
      ann_colors = list(cancer. = brewer.pal(3,"Set1")[1:3], 
                        class. = brewer.pal(5,"Set1")[4:5])
      names(ann_colors$class.) = c("cf","EV")
      names(ann_colors$cancer.) = c("NC","CRC","LUAD")
      
      col = brewer.pal(3,"Set1")
      pdf("output/result/unpaired/pheatmap_cfgenes.pdf", width = 6, height = 7)
      par(mar=c(2,2,2,2))
      pheatmap(logRPM.scale, 
               color = colorRampPalette(c(col[2],"white",col[1]))(1000),
               #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
               cutree_col = 2, 
               #cutree_row = 3, #break up the heatmap by clusters you define
               cluster_rows = F, 
               cluster_cols = F, #by default, pheatmap clusters by both row and col
               show_rownames = F,
               show_colnames = F,
               fontsize_col = 12,
               angle_col = 0,
               annotation_col = ann_col,
               #annotation_row = ann_row,
               annotation_colors = ann_colors,
               annotation_names_row = F)
      dev.off()
    }
    
  }
  ### Heterogenity
  {
    RPM <- round(edgeR::cpm(y),3)
    geneid <- rownames(RPM)
    split <- strsplit(geneid,"|",fixed=T)
    length <- unlist(lapply(split,function(split) split[2]))
    y$genes <- data.frame(length=as.numeric(length))
    RPKM <- round(rpkm(y),3)
    TPM <- round((RPKM/colSums(RPKM))*10^6,3)
    cpm <- as.data.frame(RPM)
    tpm <- as.data.frame(TPM)
    rpkm <- as.data.frame(RPKM)
    
    actb = data.frame(sample_id=colnames(cpm),
                      actb=as.numeric(cpm["ENSG00000075624.17|4405|ACTB|protein_coding",]),
                      tmsb4x=as.numeric(cpm["ENSG00000205542.11|1703|TMSB4X|protein_coding",]),
                      cancer=cancer,
                      class=class)
    actb <- actb[order(class,cancer),]
    actb$sample_id <- factor(actb$sample_id, levels=actb$sample_id)
    actb <- melt(actb, id.vars = c("sample_id","class","cancer"))
    actb <- actb %>% mutate(actb, group=paste0(cancer,"-",class,"-",variable))
    {
      p=ggplot(actb,aes(x=sample_id,y=value,color=variable))+
        #facet_grid(~class)+
        geom_point(size=1)+
        geom_line(aes(group=group))+
        scale_color_aaas()+
        scale_y_continuous(limits=c(0,40000))+
        theme_bw()+
        theme(
          plot.margin = unit(c(1, 1, 1, 0.5),"cm"),
          panel.border=element_rect(size=1, fill="transparent"),
          panel.grid=element_blank(),
          legend.position = "right",
          legend.title = element_blank(),
          legend.text= element_text(color="black", size=10),
          plot.title = element_text(hjust = 0.5, size=10, face="bold"),
          axis.text.x = element_text(color="black", size=7, angle=90),
          axis.text.y = element_text(color="black", size=10),
          axis.title.x = element_text(color="black", size=10),
          axis.title.y = element_text(color="black", size=10))+
        xlab("")+ylab("cpm")
      p
    }
    ggsave("result/actb_tmsb4x_variance.pdf",width=8,height=3.5)
    
  }
  ### PCA for all samples
  {
    library("preprocessCore")
    library("BBmisc")
    library("Rtsne")
    # normalize.quantiles: Using a normalization based upon quantiles, this function normalizes a matrix of probe level intensities.
    # normalize: Normalizes numeric data to a given scale.
    
    ###### X <- t(X)
    X <- t(edgeR::cpm(y))
    #X <- t(edgeR::cpm(y, log=T))
    
    # z-score scaling
    X_norm.4 <- scale(X, center = T, scale = T)
    
    ## Visualization
    {  
      Y <- cancer
      colors <- rainbow(length(unique(Y)))
      colors_plot <- as.character(Y)
      colors_plot[which(colors_plot == "NC")] <- colors[1]
      colors_plot[which(colors_plot == "CRC")] <- colors[2]
      colors_plot[which(colors_plot == "LUAD")] <- colors[3]
    }
    {
      Y <- class
      colors <- rainbow(length(unique(Y)))
      colors_plot <- as.character(Y)
      colors_plot[which(colors_plot == "cf")] <- colors[1]
      colors_plot[which(colors_plot == "EV")] <- colors[2]
    }
    
    # PCA
    pca <- prcomp(X_norm.4, center = F, scale = F, rank. = 2)  
    pca.var <- pca$sdev^2  ## sdev是标准偏差，十个样本，就有十个标准偏差，平方是避免负数的干扰
    pca.var.per <- round(pca.var/sum(pca.var)*100, 1)  ##求每个样本的variation
    barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")  ##用柱状图可视化
    summary(pca)
    plot(pca$x, col = colors_plot, asp = 3, pch = 20, xlab = "component_1", ylab = "component_2", main = "PCA plot")
    ## ggplot
    pca.data <- data.frame(pca$x, class=class, cancer=cancer)
    pca.data <- pca.data %>% mutate(group=paste0(class,"-",cancer))
    pca.data$group <- factor(pca.data$group, levels=c("cf-NC","EV-NC","cf-CRC","EV-CRC","cf-LUAD","EV-LUAD"))
    percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
    percentage <- paste(colnames(pca.data)[1:9],"(", paste(as.character(percentage)), "%", ")", sep="")
    {
      p=ggplot(pca.data,aes(x=PC1,y=PC2,color=cancer,shape=class)) + 
        geom_point(size=3,alpha=1) +
        #scale_color_aaas() +
        scale_color_brewer(palette = "Set1") + 
        theme_bw()+
        theme(
          plot.margin = unit(c(1, 0, 1, 1),"cm"),
          panel.border=element_rect(size=1.5, fill="transparent"),
          panel.grid=element_blank(),
          legend.position = "right",
          legend.title = element_blank(),
          legend.text= element_text(face="bold", color="black", size=15),
          plot.title = element_text(hjust = 0.5, size=12, face="bold"),
          axis.text.x = element_text(face="bold", color="black", size=15),
          axis.text.y = element_text(face="bold",  color="black", size=15),
          axis.title.x = element_text(face="bold", color="black", size=15),
          axis.title.y = element_text(face="bold",color="black", size=15))+
        xlab(percentage[1]) + ylab(percentage[2])
      p
    }
    ggsave("result/unpaired/all_PCA_cancer.pdf",p,width = 5, height = 4.5)
    
    # tSNE
    tsne <- Rtsne(X_norm.4, dims=2, check_duplicates = FALSE, perplexity = 7) 
    # dims参数设置降维之后的维度，默认值为2
    # perplexity参数的取值必须小于(nrow(data) - 1 )/ 3
    # theta参数取值越大，结果的准确度越低，默认值为0.5
    # max_iter参数设置最大迭代次数。
    plot(tsne$Y, col = colors_plot, asp = 3, pch = 20,
         xlab = "tSNE_1", ylab = "tSNE_2", main = "tSNE plot")
    
  }
  ### correlation for all samples
  {
    y <- DGEList(counts=data, group=cancer)
    y <- y[filter,]  
    dim(y)
    RPM <- round(edgeR::cpm(y),3)
    geneid <- rownames(RPM)
    split <- strsplit(geneid,"|",fixed=T)
    length <- unlist(lapply(split,function(split) split[2]))
    y$genes <- data.frame(length=as.numeric(length))
    RPKM <- round(rpkm(y),3)
    TPM <- round((RPKM/colSums(RPKM))*10^6,3)
    cpm <- as.data.frame(RPM)
    tpm <- as.data.frame(TPM)
    rpkm <- as.data.frame(RPKM)
    
    # two sample comparision
    { 
      a = "CRC-5"
      b = "CRC-5-EV"
      df <- cpm[cpm[,a]>0 & cpm[,b]>0, c(a,b)]  #CPM>0
      df <- tpm[tpm[,a]>0 & tpm[,b]>0, c(a,b)]  #TPM>0
      df <- rpkm[rpkm[,a]>0 & rpkm[,b]>0, c(a,b)] #RPKM>0
      
      df <- log10(df+1)
      # pearson correlation
      {
        cor.test(df[,a],
                 df[,b],
                 method = "pearson") # pearson比spearman算出的相关系数更高
        #label=paste("R = ",round(t$estimate,2),", p = ",signif(t$p.value,2),sep="")
        p=ggscatter(df,
                    x = a,
                    y = b,
                    #add = "reg.line", #拟合曲线
                    add.params = list(color=col[2],size=1,alpha=0.6),
                    conf.int = TRUE, #置信区间阴影带
                    cor.coef = TRUE, #系数
                    cor.method = "pearson", #方法
                    xlab = paste("log10(TPM+1), ", a, sep=""),
                    ylab = paste("log10(TPM+1), ", b, sep=""),
                    shape=16,size=2,color=col[2])
        #annotate("text",x=0.8,y=3.5,label="CPM>0",size=4,angle=0,color="black")
        #ggsave("",p,width=4,height=4)
      }
      # linear regression
      {
        fit <- lm(df[,b] ~ df[,a])
        summary <- summary(fit)
        rsq <- summary$adj.r.squared
      }
      # ggplot
      {
        ggplot(df, aes(x=`CRC-5`, y=`CRC-5-EV`))+
          geom_point(size = 2, shape = 16, color = col[2], alpha = 0.7)+
          geom_smooth(method = lm, size = 0.5, color = col[2])+
          #geom_abline(aes(intercept=0,slope=1),linetype=5,size=0.5, color=col[2])+
          #scale_x_continuous(breaks=c(0,3,6,9,12))+
          #scale_y_continuous(breaks=c(0,3,6,9,12))+
          theme_bw()+
          theme(
            #plot.margin = unit(c(1, 1, 1, 1),"cm"),
            legend.position=c(0.1,0.6),
            legend.title = element_text(color="black", size=12),
            legend.text= element_text(color="black", size=12),
            panel.grid = element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_text(color="black", size=12),
            axis.text.y = element_text(color="black", size=12),
            axis.title.x = element_text(color="black", size=12),
            axis.title.y = element_text(color="black", size=12),
            axis.line = element_line(size=0.5, colour = "black"))+
          labs(x="pico(logTPM)",y="ipico(logTPM)",face="bold")
        #annotate("text",x=3.5,y=12,label="R = 0.76, p < 2.2e-16",size=5,angle=0,color="black")
        annotate("text",x=1.3,y=6.5,label=bquote(R^2 ~ " = " ~.(rsq)),size=5,angle=0,color="black",fontface="bold")
        p
      }
    }
    
    # corrplot
    {
      # plot
      library(corrplot)
      library(psych)
      library(Hmisc)
      
      #df <- tpm
      df <- cpm
      #df <- rpkm
      #df <- log10(df+1)
      
      df <- df[,order(class, cancer)]
      class. <- class[order(class,cancer)]
      cancer. <- cancer[order(class,cancer)]
      
      #class.=class
      #cancer.=cancer
      
      # correlation
      rcorr <- rcorr(as.matrix(df), type="pearson")
      # R square
      {
        rsq=list()
        n=1
        for (i in seq(1,9,1)){
          for (j in seq(1,9,1)){
            fit = lm(df[,i] ~ df[,j])
            rsq[n] = summary(fit)$adj.r.squared
            n=n+1}}
        rsq <- matrix(unlist(rsq), nrow = 9, ncol = 9)
        colnames(rsq) <- c("H4R1","H4R2","H4R3","C1","C2","C3","N1","N2","N3")
        rownames(rsq) <- c("H4R1","H4R2","H4R3","C1","C2","C3","N1","N2","N3")
      }
      # corrplot
      {
        #col1 <- colorRampPalette(c("#6D9EC1", "white", "#E46726"))
        col1 <- colorRampPalette(c("blue", "white", "red"))
        pdf("result/paired/corrplot-cpm.pdf",width = 5,height = 5)
        par(oma=c(1,0,1,1),mar=c(2,0,2,2))
        # rcorr$r
        corrplot(rcorr$r, 
                 #method = "square",
                 type="lower",
                 #col=col[2],
                 is.corr=F,
                 order="original",
                 addrect = 4,     
                 insig="p-value",
                 tl.pos="lt",
                 tl.col="black", 
                 tl.cex=0.5,
                 tl.srt=90, 
                 cl.pos="r",
                 #cl.lim=c(0,1),
                 #addCoef.col="black",
                 addCoefasPercent=F) 
        #没有办法同时自定义颜色+自定义范围
        dev.off()
      }
      # corrplot-rsq
      {
        pdf("reproducibility/corrplot/corrplot-rsq-rpkm.no.lib.size.pdf",width = 5,height = 5)
        par(oma=c(1,0,1,1),mar=c(2,0,2,2))
        # rsq
        corrplot(rsq, 
                 method = "color",
                 type="lower",
                 is.corr=F, 
                 order="original",
                 addrect = 4,
                 tl.pos="lt",
                 tl.col="black", 
                 tl.srt=90, 
                 cl.pos="b",
                 #cl.lim=c(0,1),
                 addCoef.col="black",
                 addCoefasPercent=F)  
        #没有办法同时自定义颜色+自定义范围
        dev.off()
      }
      ## pheatmap::pheatmap 
      {
        library(pheatmap)
        ann_col = data.frame(cancer = as.character(cancer.),
                             class = as.character(class.))
        rownames(ann_col) = colnames(rcorr$r)
        ann_row = data.frame(cancer = as.character(cancer.),
                             class = as.character(class.))
        rownames(ann_row) = rownames(rcorr$r)
        ann_colors = list(class = brewer.pal(6,"Dark2")[1:2],
                          cancer = brewer.pal(6,"Dark2")[3:5])
        names(ann_colors$class) = c("cf","EV")
        names(ann_colors$cancer) = c("NC","CRC","LUAD")
        
        col = brewer.pal(8,"Blues")
        pdf("result/unpaired/corr-heatmap-cpm.pdf", width = 6, height = 5)
        par(mar=c(2,2,2,2))
        pheatmap(rcorr$r, 
                 color = colorRampPalette(c("grey","white",col[8]))(1000),
                 #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                 #cutree_col = 2, 
                 #cutree_row = 3, #break up the heatmap by clusters you define
                 cluster_rows = F, 
                 cluster_cols = F, #by default, pheatmap clusters by both row and col
                 show_rownames = F,
                 show_colnames = F,
                 fontsize_col = 12,
                 angle_col = 0,
                 annotation_col = ann_col,
                 annotation_row = ann_row,
                 annotation_colors = ann_colors,
                 annotation_names_row = F)
        dev.off()
      }
    }
    
    
    
    
    
  }
  ### enrichment analysis
  {
    library(clusterProfiler)
    library(enrichplot)
    #### input table
    #de.order <- de.ev$de
    de.order <- de.cf$de
    
    de.order$ENSG <- rownames(de.order)
    de.order$ENSG <- as.character(lapply(strsplit(de.order$ENSG,"|",fixed=TRUE), function(x) x[1]))
    de.order$ENSG <- as.character(lapply(strsplit(de.order$ENSG,".",fixed=TRUE), function(x) x[1]))
    de.order <- de.order[!duplicated(de.order$ENSG),]
    
    ##### get all geneList (for GSEA)
    #gene.list.ncbiid <- bitr(de.order$name, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    gene.list.ncbiid <- read.table("data/gencode.v38.txt",sep = "\t",header = T)
    gene.list.ncbiid <- gene.list.ncbiid[!is.na(gene.list.ncbiid$ENTREZID),]
    gene.list.ncbiid <- gene.list.ncbiid[,c("ensg","ENTREZID")]
    gene.list.ncbiid <- gene.list.ncbiid[!duplicated(gene.list.ncbiid$ensg) & !duplicated(gene.list.ncbiid$ENTREZID),]
    colnames(gene.list.ncbiid)[1] <- "ENSG"
    
    #de.order <- merge(x = de.order, y = gene.list.ncbiid, by.x="name", by.y="ENSG")
    de.order <- dplyr::left_join(x = de.order, y = gene.list.ncbiid)
    de.order <- de.order[order(de.order$logFC, -de.order$PValue, decreasing = T),] # order by logFC then PValue
    
    gene.list.ncbiid <- de.order$logFC[!duplicated(de.order$ENTREZID)]
    names(gene.list.ncbiid) <- de.order$ENTREZID[!duplicated(de.order$ENTREZID)]
    message(paste0("remaining ENTREZID geneList lenth: ",length(gene.list.ncbiid)))
    gene.list.ensgid <- de.order$logFC #[!duplicated(de.order$ENTREZID)]
    names(gene.list.ensgid) <- de.order$ENSG #[!duplicated(de.order$ENSG)]
    message(paste0("remaining ENSEMBL geneList lenth: ",length(gene.list.ensgid)))
    
    
    ##### get top geneList (for ORA)
    method="cutoff.FDR" # for all samples
    method="cutoff.pvalue" # for each cancer group
    {
      if (method=="cutoff.FDR"){
        cutoff.FDR = 0.1
        cutoff.logFC = 1
        message(paste0("FDR cutoff: ", cutoff.FDR, "  logFC cutoff: ", cutoff.logFC))
        gene.top.ncbiid <- as.character(na.omit(de.order$ENTREZID[abs(de.order$logFC) > cutoff.logFC & de.order$FDR <= cutoff.FDR]))   # FDR cutoff
        gene.top.ncbiid.up <- as.character(na.omit(de.order$ENTREZID[de.order$logFC > cutoff.logFC & de.order$FDR <= cutoff.FDR]))
        gene.top.ncbiid.down <- as.character(na.omit(de.order$ENTREZID[de.order$logFC < -cutoff.logFC & de.order$FDR <= cutoff.FDR]))
        message(paste0("remaining ENTREZID top gene lenth: ",length(gene.top.ncbiid)))
        message(paste0("remaining ENTREZID top up gene lenth: ",length(gene.top.ncbiid.up)))
        message(paste0("remaining ENTREZID top down gene lenth: ",length(gene.top.ncbiid.down)))
        gene.top.ensgid <- de.order$ENSG[abs(de.order$logFC) > cutoff.logFC & de.order$FDR <= cutoff.FDR]   # FDR cutoff
        gene.top.ensgid.up <- de.order$ENSG[de.order$logFC > cutoff.logFC & de.order$FDR <= cutoff.FDR] 
        gene.top.ensgid.down <- de.order$ENSG[de.order$logFC < -cutoff.logFC & de.order$FDR <= cutoff.FDR]
        message(paste0("remaining ENSEMBL top gene lenth: ",length(gene.top.ensgid)))
        message(paste0("remaining ENSEMBL top up gene lenth: ",length(gene.top.ensgid.up)))
        message(paste0("remaining ENSEMBL top down gene lenth: ",length(gene.top.ensgid.down)))
      }
      else if (method=="cutoff.pvalue"){
        cutoff.pvalue = 0.01   #??????
        cutoff.logFC = 1
        message(paste0("p-value cutoff: ", cutoff.pvalue, "  logFC cutoff: ", cutoff.logFC))
        gene.top.ncbiid <- as.character(na.omit(de.order$ENTREZID[abs(de.order$logFC) > cutoff.logFC & de.order$PValue <= cutoff.pvalue]))   # p-value cutoff
        gene.top.ncbiid.up <- as.character(na.omit(de.order$ENTREZID[de.order$logFC > cutoff.logFC & de.order$PValue <= cutoff.pvalue])) 
        gene.top.ncbiid.down <- as.character(na.omit(de.order$ENTREZID[de.order$logFC < -cutoff.logFC & de.order$PValue <= cutoff.pvalue]))
        message(paste0("remaining ENTREZID top gene lenth: ",length(gene.top.ncbiid)))
        message(paste0("remaining ENTREZID top up gene lenth: ",length(gene.top.ncbiid.up)))
        message(paste0("remaining ENTREZID top down gene lenth: ",length(gene.top.ncbiid.down)))
        gene.top.ensgid <- de.order$ENSG[abs(de.order$logFC) > cutoff.logFC & de.order$PValue <= cutoff.pvalue]   # FDR cutoff
        gene.top.ensgid.up <- de.order$ENSG[de.order$logFC > cutoff.logFC & de.order$PValue <= cutoff.pvalue] 
        gene.top.ensgid.down <- de.order$ENSG[de.order$logFC < -cutoff.logFC & de.order$PValue <= cutoff.pvalue]
        message(paste0("remaining ENSEMBL top gene lenth: ",length(gene.top.ensgid)))
        message(paste0("remaining ENSEMBL top up gene lenth: ",length(gene.top.ensgid.up)))
        message(paste0("remaining ENSEMBL top down gene lenth: ",length(gene.top.ensgid.down)))
      }
      else if (method=="FDR.top"){
        FDR.top = 100
        message(paste0("FDR top: ", FDR.top))
        de.order <- de.order[order(-de.order$FDR, decreasing = T),]  # order by FDR
        gene.top.ncbiid <- as.character(na.omit(de.order$ENTREZID[1:FDR.top]))
        gene.top.ncbiid.up <- as.character(na.omit(de.order$ENTREZID[de.order$logFC > 0][1:FDR.top]))
        gene.top.ncbiid.down <- as.character(na.omit(de.order$ENTREZID[de.order$logFC < 0][1:FDR.top]))
        message(paste0("remaining ENTREZID top gene lenth: ",length(gene.top.ncbiid)))
        message(paste0("remaining ENTREZID top up gene lenth: ",length(gene.top.ncbiid.up)))
        message(paste0("remaining ENTREZID top down gene lenth: ",length(gene.top.ncbiid.down)))
        gene.top.ensgid <- de.order$ENSG[1:FDR.top]   
        gene.top.ensgid.up <- de.order$ENSG[1:FDR.top]
        gene.top.ensgid.down <- de.order$ENSG[1:FDR.top]
        message(paste0("remaining ENSEMBL top gene lenth: ",length(gene.top.ensgid)))
        message(paste0("remaining ENSEMBL top up gene lenth: ",length(gene.top.ensgid.up)))
        message(paste0("remaining ENSEMBL top down gene lenth: ",length(gene.top.ensgid.down)))
      }
      else{
        message(paste0("No methods specified"))
      }
    }
    
    pValue=0.1
    qValue=0.1
    pAdjust="BH"
    # 1. Over-Representation Analysis 超几何检验或Fisher精确检验
    message("start ORA enrich...")
    # 1.2 ORA up gene
    {
      if(FALSE){
        enrich.GO.up <- enrichGO(
          gene = gene.top.ncbiid.up,
          #universe = names(gene.list.ncbiid),  #选择background gene?
          OrgDb = org.Hs.eg.db,
          ont = "ALL",
          keyType = "ENTREZID",
          pvalueCutoff = pValue,
          pAdjustMethod = pAdjust,
          qvalueCutoff = qValue)
      }
      enrich.GO.up <- enrichGO(
        gene = gene.top.ensgid.up,
        #universe = names(gene.list.ncbiid),  #选择background gene?
        OrgDb = org.Hs.eg.db,
        ont = "ALL",
        keyType = "ENSEMBL",
        pvalueCutoff = pValue,
        pAdjustMethod = pAdjust,
        qvalueCutoff = qValue
      )
      enrich.DO.up <- enrichDO(
        gene = gene.top.ncbiid.up,
        ont = "DO",
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust,
        qvalueCutoff = 1,
        #universe = names(gene.list.ncbiid),
        minGSSize = 5,
        maxGSSize = 500,
        readable = FALSE
      )
      enrich.KEGG.up <- enrichKEGG(
        gene = gene.top.ncbiid.up,
        keyType = "ncbi-geneid",
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust,
        qvalueCutoff = 1
      )
    }
    # 1.3 ORA down gene
    {
      if(FALSE){
        enrich.GO.down <- enrichGO(
          gene = gene.top.ncbiid.down,
          #universe = names(gene.list.ncbiid),  #选择background gene?
          OrgDb = org.Hs.eg.db,
          ont = "ALL",
          keyType = "ENTREZID",
          pvalueCutoff = pValue,
          pAdjustMethod = pAdjust,
          qvalueCutoff = qValue)
      }
      enrich.GO.down <- enrichGO(
        gene = gene.top.ensgid.down,
        #universe = names(gene.list.ncbiid),  #选择background gene?
        OrgDb = org.Hs.eg.db,
        ont = "ALL",
        keyType = "ENSEMBL",
        pvalueCutoff = pValue,
        pAdjustMethod = pAdjust,
        qvalueCutoff = qValue
      )
      enrich.DO.down <- enrichDO(
        gene = gene.top.ncbiid.down,
        ont = "DO",
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust,
        qvalueCutoff = 1,
        #universe = names(gene.list.ncbiid),
        minGSSize = 5,
        maxGSSize = 500,
        readable = FALSE
      )
      enrich.KEGG.down <- enrichKEGG(
        gene = gene.top.ncbiid.down,
        keyType = "ncbi-geneid",
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust,
        qvalueCutoff = 1
      )
    }
    # output table
    {
      types <- c(
                 enrich.GO.up,enrich.DO.up,enrich.KEGG.up,
                 enrich.GO.down,enrich.DO.down,enrich.KEGG.down)
      names(types) <-c(
                       "enrich.GO.up","enrich.DO.up","enrich.KEGG.up",
                       "enrich.GO.down","enrich.DO.down","enrich.KEGG.down")
      
      outdir <- paste0("output/result/unpaired/enrich/cf-table")
      dir.create(paste0(outdir), showWarnings = T)
      message("start writing tables...")
      for(i in 1:length(types)){
        write.table(types[[i]],paste0(outdir,"/", names(types)[i], ".txt"),row.names = F,col.names = T, quote = F,sep = "\t")}
    }
    # enrichment plot of all samples
    {
      options(scipen=999) 
      outdir="output/result/unpaired/enrich/cf-plot/"
      col = brewer.pal(3,"Set1")
      
      pdf(paste0(outdir,"dotplot-cf.KEGG.up",".pdf"), width = 6, height = 4)
      enrichplot::dotplot(enrich.KEGG.up, font.size=12, title="cf.KEGG.up", showCategory = 7) +
        scale_color_continuous(low=col[1],high=col[2]) 
        scale_x_continuous(limits = c(0.03,0.1), breaks = c(0.04,0.06,0.08))
      dev.off()
      pdf(paste0(outdir,"dotplot-cf.KEGG.down",".pdf"), width = 6, height = 4)
      enrichplot::dotplot(enrich.KEGG.down, font.size=14, title="cf.KEGG.down", showCategory = 3) +
        scale_color_continuous(low=col[1],high=col[2]) 
        scale_x_continuous(limits = c(0.04,0.125), breaks = c(0.05,0.1))
      dev.off()
      pdf(paste0(outdir,"dotplot-cf.DO.up",".pdf"), width = 6, height = 6)
      enrichplot::dotplot(enrich.DO.up, font.size=12, title="cf.DO.up", showCategory = 20) +
        scale_color_continuous(low=col[1],high=col[2]) 
      scale_x_continuous(limits = c(0.03,0.1), breaks = c(0.04,0.06,0.08))
      dev.off()
      
      {
        ev.go.up <- read.table("output/result/unpaired/enrich/ev-table/enrich.GO.down.txt",sep="\t",header=T,quote="")
        ev.go.up.ID <- ev.go.up[order(ev.go.up$p.adjust),][1:10,"ID"]
        
        cf.go.up <- read.table("output/result/unpaired/enrich/cf-table/enrich.GO.down.txt",sep="\t",header=T,quote="")
        cf.go.up.ID <- cf.go.up[order(cf.go.up$p.adjust),][1:10,"ID"]
        
        go.up <- rbind(ev.go.up[ev.go.up$ID %in% ev.go.up.ID,] %>% mutate(class="EV"),
                       cf.go.up[cf.go.up$ID %in% cf.go.up.ID,] %>% mutate(class="cf"))
        
        go.up$logp <- -log10(go.up$p.adjust)
        go.up$ONTOLOGY <- factor(go.up$ONTOLOGY, levels=c("CC","BP","MF"))
        go.up <- go.up[order(go.up$ONTOLOGY,go.up$p.adjust),]
        des <- go.up$Description[!duplicated(go.up$Description)]
        go.up$Description <- factor(go.up$Description, levels=rev(des))
      }
      {        
        ev.kegg.up <- read.table("output/result/unpaired/enrich/ev-table/enrich.KEGG.down.txt",sep="\t",header=T,quote="")
        ev.kegg.up.ID <- ev.kegg.up[order(ev.kegg.up$p.adjust),][1:4,"ID"]
        
        cf.kegg.up <- read.table("output/result/unpaired/enrich/cf-table/enrich.KEGG.down.txt",sep="\t",header=T,quote="")
        cf.kegg.up.ID <- cf.kegg.up[order(cf.kegg.up$p.adjust),][1:5,"ID"]
        
        kegg.up <- rbind(ev.kegg.up[ev.kegg.up$ID %in% ev.kegg.up.ID,] %>% mutate(class="EV"),
                       cf.kegg.up[cf.kegg.up$ID %in% cf.kegg.up.ID,] %>% mutate(class="cf"))
        
        kegg.up$logp <- -log10(kegg.up$p.adjust)
        kegg.up <- kegg.up[order(kegg.up$p.adjust),]
        des <- kegg.up$Description[!duplicated(kegg.up$Description)]
        kegg.up$Description <- factor(kegg.up$Description, levels=rev(des))
      }

      {
        p=ggplot(kegg.up, aes(class, Description)) +
          geom_point(aes(color = p.adjust, size = Count)) +
          #scale_color_aaas()+
          #scale_x_continuous(limits=c(0.9,1.6))+
          scale_color_gradient(low=col[1],high=col[2],guide="colorbar")+
          theme_bw()+
          guides(fill=guide_legend(title=NULL, ncol=1))+
          #coord_fixed(ratio=0.3)+
          theme(
            plot.margin = unit(c(1, 1, 1, 1),"cm"),
            panel.border=element_rect(size=1, fill="transparent"),
            panel.grid=element_line(size=rel(0.5)),
            #legend.position=c(0.25,0.85),
            legend.position="right",
            legend.background = element_blank(),
            legend.title = element_text(color="black", size=10),
            legend.text= element_text(color="black", size=10),
            axis.text.x = element_text(color="black", size=10),
            axis.text.y = element_text(color="black", size=10),
            axis.title.x = element_text(color="black", size=10),
            axis.title.y = element_blank())+
          labs(title="KEGG")+xlab("")
        p
      }
      outdir="output/result/unpaired/enrich/"
      ggsave(paste0(outdir,"dotplot-cf-ev-KEGG.down",".pdf"), width = 5, height = 5)
      
      
      # DAG 有向无环图
      plotGOgraph(enrich.GO)
      # igraph布局的DAG
      goplot(enrich.GO)
      
      barplot(enrich.GO, showCategory = 10)
      barplot(enrich.DO, showCategory = 10)
      barplot(enrich.KEGG, showCategory = 10)
    }

    # 2.Gene Set Enrichment Analysis online
    message("start gsea...")
    {
      gse.GO <- gseGO(
        geneList = gene.list.ensgid, 
        OrgDb = org.Hs.eg.db,
        ont = "ALL", 
        keyType = "ENSEMBL",  
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust
        #by = "fgsea", #'DOSE'
      ) 
      if(FALSE){
        gse.GO <- gseGO(
          geneList = gene.list.ncbiid, 
          OrgDb = org.Hs.eg.db,
          ont = "ALL", 
          keyType = "ENTREZID",  
          pvalueCutoff = pValue,
          pAdjustMethod = pAdjust
          #by = "fgsea", #'DOSE'
        ) 
      }
      gse.KEGG <- gseKEGG(
        geneList = gene.list.ncbiid,
        keyType = 'kegg', 
        organism = 'hsa',
        nPerm  = 1000,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust
        #by = "fgsea",
        #seed = T  
      )
      gse.DO <- gseDO(
        geneList = gene.list.ncbiid,
        nPerm  = 1000,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust
        #by = "fgsea",
        #seed = T
      )
    }
    # 2.1 KEGG Gene Set Enrichment Analysis 
    message("start gsea at specific geneset")
    {
      m_t <- read.table("data/PATH_ID_NAME_KEGGplusHallmark.txt",sep = "\t",header = T, stringsAsFactors = F,row.names = 1)
      m_t2g <- m_t[,c("short_DESCRPTION","ensembl_gene_id")]
      #m_t2g <- m_t2g[!duplicated(m_t2g$short_DESCRPTION),]
      m_t2e <- m_t[,c("short_DESCRPTION","ENTREZID")]
      #m_t2e <- m_t2e[!duplicated(m_t2e$short_DESCRPTION),]
      gsea <- GSEA(gene.list.ncbiid,
                   TERM2GENE = m_t2e, 
                   nPerm  = 1000,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = pValue,
                   pAdjustMethod = pAdjust)
      gsea2 <- setReadable(gsea.kegg,
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID")
    }
    # output table
    {
      types <- c(gse.GO, gse.KEGG, gse.DO)
      names(types) <-c("gse.GO","gse.KEGG","gse.DO")
      
      outdir <- paste0("output/result/unpaired/enrich/cf-table")
      dir.create(paste0(outdir), showWarnings = T)
      message("start writing tables...")
      for(i in 1:length(types)){
        write.table(types[[i]],paste0(outdir,"/", names(types)[i], ".txt"),row.names = F,col.names = T, quote = F,sep = "\t")}
    }
    # GSEA plot of all samples
    {
      outdir="output/result/unpaired/enrich/ev-plot/"
      col = brewer.pal(8,"Set1")
      # 气泡图，展示geneset被激活还是抑制
      {
        pdf(paste0(outdir,"dotplot-gse.GO",".pdf"), width = 12, height = 5.8)
        enrichplot::dotplot(gse.GO, split=".sign", font.size=10, title="gsea.GO", showCategory = 9) +
          facet_grid(~.sign)# + scale_color_continuous(low=col[1],high=col[2]) 
        scale_x_continuous(limits = c(0.18,1.85),breaks = c(0.2,0.4,0.6,0.8))
        dev.off()
        
        pdf(paste0(outdir,"dotplot-gse.KEGG",".pdf"), width = 7, height = 6)
        enrichplot::dotplot(gse.KEGG, split=".sign", font.size=10, title="gsea.KEGG", showCategory = 10) +
          facet_grid(~.sign) #+ scale_color_continuous(low=col[1],high=col[2]) 
        scale_x_continuous(limits = c(0.25,0.55),breaks = c(0.3,0.4,0.5))
        dev.off()
      }
      # multi-group 
      {
        ev.go <- read.table("output/result/unpaired/enrich/ev-table/gse.GO.txt", header=T, sep="\t")
        ev.go <- ev.go %>% mutate(Count = str_count(core_enrichment,"/")+1)
        ev.go.up <- ev.go[ev.go$NES > 1,]
        ev.go.up.order <- ev.go.up[order(-ev.go.up$p.adjust, ev.go.up$NES, decreasing = T),]
        ev.go.up.order <- ev.go.up.order[1:10,]

        ev.go.down <- ev.go[ev.go$NES < -1,]
        ev.go.down.order <- ev.go.down[order(-ev.go.down$p.adjust, -ev.go.down$NES, decreasing = T),]
        ev.go.down.order <- ev.go.down.order[1:10,]
 
        cf.go <- read.table("output/result/unpaired/enrich/cf-table/gse.GO.txt", header=T, sep="\t")
        cf.go <- cf.go %>% mutate(Count = str_count(core_enrichment,"/")+1)
        cf.go.up <- cf.go[cf.go$NES > 1,]
        cf.go.up.order <- cf.go.up[order(-cf.go.up$p.adjust, cf.go.up$NES, decreasing = T),]
        cf.go.up.order <- cf.go.up.order[1:10,]

        cf.go.down <- cf.go[cf.go$NES < -1,]
        cf.go.down.order <- cf.go.down[order(-cf.go.down$p.adjust, -cf.go.down$NES, decreasing = T),]
        cf.go.down.order <- cf.go.down.order[1:10,]

        go.up <- rbind(ev.go.up.order %>% mutate(class = "EV"),
                       cf.go.up.order %>% mutate(class = "cf"))
        go.up <- rbind(ev.go.down.order %>% mutate(class = "EV"),
                       cf.go.down.order %>% mutate(class = "cf"))
        go.up$ONTOLOGY <- factor(go.up$ONTOLOGY, levels=c("CC","BP","MF"))
        go.up <- go.up[order(go.up$ONTOLOGY,go.up$p.adjust),]
        des <- go.up$Description[!duplicated(go.up$Description)]
        go.up$Description <- factor(go.up$Description, levels=rev(des))
    
        
        
        {
          p=ggplot(go.up, aes(class, Description)) +
            geom_point(aes(color = p.adjust, size = Count)) +
            #scale_color_aaas()+
            #scale_x_continuous(limits=c(0.9,1.6))+
            scale_color_gradient(low=col[1],high=col[2],guide="colorbar")+
            theme_bw()+
            guides(fill=guide_legend(title=NULL, ncol=1))+
            #coord_fixed(ratio=0.3)+
            theme(
              plot.margin = unit(c(1, 1, 1, 1),"cm"),
              panel.border=element_rect(size=1, fill="transparent"),
              panel.grid=element_line(size=rel(0.5)),
              #legend.position=c(0.25,0.85),
              legend.position="right",
              legend.background = element_blank(),
              legend.title = element_text(color="black", size=10),
              legend.text= element_text(color="black", size=10),
              axis.text.x = element_text(color="black", size=10),
              axis.text.y = element_text(color="black", size=10),
              axis.title.x = element_text(color="black", size=10),
              axis.title.y = element_blank())+
            labs(title="GO")+xlab("")
          p
        }
        outdir="output/result/unpaired/enrich/"
        ggsave(paste0(outdir,"dotplot-cf-ev-GSEA-go.down",".pdf"), width = 12, height = 5)
        
      }
      {
        ev.kegg <- read.table("output/result/unpaired/enrich/ev-table/gse.KEGG.txt", header=T, sep="\t")
        ev.kegg <- ev.kegg %>% mutate(Count = str_count(core_enrichment,"/")+1)
        ev.kegg.up <- ev.kegg[ev.kegg$NES > 1,]
        ev.kegg.up.order <- ev.kegg.up[order(-ev.kegg.up$p.adjust, ev.kegg.up$NES, decreasing = T),]
        ev.kegg.up.order <- ev.kegg.up.order[1:10,]
        
        ev.kegg.down <- ev.kegg[ev.kegg$NES < -1,]
        ev.kegg.down.order <- ev.kegg.down[order(-ev.kegg.down$p.adjust, -ev.kegg.down$NES, decreasing = T),]
        ev.kegg.down.order <- ev.kegg.down.order[1:10,]
        
        cf.kegg <- read.table("output/result/unpaired/enrich/cf-table/gse.KEGG.txt", header=T, sep="\t")
        cf.kegg <- cf.kegg %>% mutate(Count = str_count(core_enrichment,"/")+1)
        cf.kegg.up <- cf.kegg[cf.kegg$NES > 1,]
        cf.kegg.up.order <- cf.kegg.up[order(-cf.kegg.up$p.adjust, cf.kegg.up$NES, decreasing = T),]
        cf.kegg.up.order <- cf.kegg.up.order[1:10,]
        
        cf.kegg.down <- cf.kegg[cf.kegg$NES < -1,]
        cf.kegg.down.order <- cf.kegg.down[order(-cf.kegg.down$p.adjust, -cf.kegg.down$NES, decreasing = T),]
        cf.kegg.down.order <- cf.kegg.down.order[1:10,]
        
        kegg.up <- rbind(ev.kegg.up.order %>% mutate(class = "EV"),
                       cf.kegg.up.order %>% mutate(class = "cf"))
        kegg.up <- rbind(ev.kegg.down.order %>% mutate(class = "EV"),
                       cf.kegg.down.order %>% mutate(class = "cf"))
        kegg.up <- kegg.up[order(kegg.up$p.adjust),]
        des <- kegg.up$Description[!duplicated(kegg.up$Description)]
        kegg.up$Description <- factor(kegg.up$Description, levels=rev(des))
        
        
        
        {
          p=ggplot(kegg.up, aes(class, Description)) +
            geom_point(aes(color = p.adjust, size = Count)) +
            #scale_color_aaas()+
            #scale_x_continuous(limits=c(0.9,1.6))+
            scale_color_gradient(low=col[1],high=col[2],guide="colorbar")+
            theme_bw()+
            guides(fill=guide_legend(title=NULL, ncol=1))+
            #coord_fixed(ratio=0.3)+
            theme(
              plot.margin = unit(c(1, 1, 1, 1),"cm"),
              panel.border=element_rect(size=1, fill="transparent"),
              panel.grid=element_line(size=rel(0.5)),
              #legend.position=c(0.25,0.85),
              legend.position="right",
              legend.background = element_blank(),
              legend.title = element_text(color="black", size=10),
              legend.text= element_text(color="black", size=10),
              axis.text.x = element_text(color="black", size=10),
              axis.text.y = element_text(color="black", size=10),
              axis.title.x = element_text(color="black", size=10),
              axis.title.y = element_blank())+
            labs(title="kegg")+xlab("")
          p
        }
        outdir="output/result/unpaired/enrich/"
        ggsave(paste0(outdir,"dotplot-cf-ev-GSEA-kegg.down",".pdf"), width = 10, height = 5)
        
      }
      
      # GSEA plot
      {
        col = brewer.pal(8,"Dark2")

        up.id = c("GO:0031091","GO:0030168","GO:0007599",
                  "GO:0005840","GO:0000184","GO:0005844")
        
        pdf(paste0(outdir,"gseaplot.cf.GO",".pdf"), width = 8, height = 6)
        gseaplot2(gse.GO, geneSetID = up.id, color = col[1:6], pvalue_table = F, 
                  base_size = 12, rel_heights = c(1.5, 0.4, 0.6), subplots = 1:3)
        dev.off()
        
        
        pdf(paste0(outdir,"gseaplot.KEGG",".pdf"), width = 8, height = 6)
        gseaplot2(gse.KEGG, geneSetID = 1:5, color = col[1:5], pvalue_table = F, 
                  base_size = 10, rel_heights = c(1.5, 0.5, 0.8), subplots = 1:3)
        dev.off()
        
        pdf(paste0(outdir,"gseaplot.DO",".pdf"), width = 8, height = 6)
        gseaplot2(gse.DO, geneSetID = 1:5, color = col[1:5], pvalue_table = T, 
                  base_size = 10, rel_heights = c(1.5, 0.5, 0.8), subplots = 1:3)
        dev.off()
      }

    }

    # 2.2 Cancer gene set Enrichment Analysis
    message("start gsea at cancer gene set")
    {
      dir="/BioII/lulab_b/zhanqing/cancer_gene/"
      CGC_1 <- read.table(paste0(dir,"v3/id/CGC_1.txt"), sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      CGC_2 <- read.table(paste0(dir,"v3/id/CGC_2.txt"), sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      CGC_all <- read.table(paste0(dir,"v3/id/CGC_all.txt"), sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      NCG_1 <- read.table(paste0(dir,"v3/id/NCG_1.txt"), sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      NCG_2 <- read.table(paste0(dir,"v3/id/NCG_2.txt"), sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      Uniprot <- read.table(paste0(dir,"v3/id/Uniprot.txt"), sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      Pansoft <- read.table(paste0(dir,"v3/id/PanSoftware.txt"), sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      IntOGene <- read.table(paste0(dir,"v3/id/IntOGen.txt"), sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      
      # immune gene set
      setTcellExaustion <- read.table(paste0(dir,"v3/ImmuneGene/T cell Exaustion genes.csv"), sep=",", header=T)
      setTcellExaustion <- setTcellExaustion$ensembl
      setRPmRNA <- read.table(paste0(dir,"v3/ImmuneGene/RP_mRNA.csv"), sep=",", header=T, check.names = F)
      setRPmRNA <- setRPmRNA$`Gene ID`
      setImmuneModu <- read.table(paste0(dir,"v3/ImmuneGene/78_ImmuneModulatory_Genes_full.csv"), sep=",", header=T, check.names = F)
      setImmuneModu <- setImmuneModu$ensembl
      setEVtop <- read.table(paste0(dir,"v3/EVGenes/EV_TOP_100_ensembl.txt"), sep=",", header=F, check.names = F)[,1]
      
      set1 <- read.table("v2/public/ensemblid/set/Set1.txt")[,1]
      set2 <- read.table("v2/public/ensemblid/set/Set2.txt")[,1]
      set3 <- read.table("v2/public/ensemblid/set/Set3.txt")[,1]
      set4 <- read.table("v2/public/ensemblid/set/Set4.txt")[,1]
      set5 <- read.table("v2/public/ensemblid/set/Set5.txt")[,1]
      set6 <- read.table("v2/public/ensemblid/set/Set6.txt")[,1]
      set7 <- read.table("v2/public/ensemblid/set/Set7.txt")[,1]
      
      geneset1 <- rbind(data.frame(Description = "CGC_1", Entrez = CGC_1$Entrez),
                       #data.frame(Description = "CGC_all", Entrez = CGC_all$Entrez),
                       data.frame(Description = "NCG_1", Entrez = NCG_1$ENTREZID),
                       #data.frame(Description = "NCG_all", Entrez = NCG_2$ENTREZID),
                       data.frame(Description = "Uniprot", Entrez = Uniprot$ENTREZID),
                       data.frame(Description = "Pansoft", Entrez = Pansoft$ENTREZID),
                       data.frame(Description = "IntOGene", Entrez = IntOGene$ENTREZID))
      geneset2 <- rbind(data.frame(Description = "set1", ensg = set1),
                        data.frame(Description = "set2", ensg = set2),
                        data.frame(Description = "set3", ensg = set3),
                        data.frame(Description = "set4", ensg = set4),
                        data.frame(Description = "set5", ensg = set5),
                        data.frame(Description = "set6", ensg = set6),
                        data.frame(Description = "set7", ensg = set7))
      geneset3 <- rbind(data.frame(Description = "TcellExaustion", ensg = setTcellExaustion),
                        data.frame(Description = "RPmRNA", ensg = setRPmRNA),
                        data.frame(Description = "ImmuneModu", ensg = setImmuneModu),
                        data.frame(Description = "EVtop", ensg = setEVtop)
      )
      geneset <- rbind(geneset2,geneset3)
      pValue = 1
      pAdjust = "BH"
      gsea <- clusterProfiler::GSEA(gene.list.ensgid,    #  ensg for geneset2, ncbi for geneset1
                                    TERM2GENE = geneset,
                                    nPerm  = 1000,
                                    minGSSize = 10,
                                    maxGSSize = 500,
                                    pvalueCutoff = pValue,
                                    pAdjustMethod = pAdjust)
      #write.table(gsea,  "result/unpaired/enrich/cf-table/gsea.geneset.txt",row.names = F,col.names = T, quote = F,sep = "\t")
    
      outdir="result/unpaired/enrich/"
      col = brewer.pal(8,"Set1")
      pdf(paste0(outdir,"gseaplot.ev.geneset",".pdf"), width = 8, height = 6)
      gseaplot2(gsea, geneSetID = c("set1","set3","set5","RPmRNA","TcellExaustion"), color = col[1:5], pvalue_table = F, 
                base_size = 12, rel_heights = c(1.5, 0.4, 0.7), subplots = 1:3)
      dev.off()
      
      gse.cf <- read.table("result/unpaired/enrich/cf-table/gsea.geneset.txt",sep="\t",header = T)
      gse.ev <- read.table("result/unpaired/enrich/ev-table/gsea.geneset.txt",sep="\t",header = T)
      gse <- rbind(gse.cf %>% mutate(class = "cf", Count = str_count(core_enrichment,"/")+1),
                   gse.ev %>% mutate(class = "EV", Count = str_count(core_enrichment,"/")+1))
      gse <- gse[gse$Description!="EVtop",]
      gse$Description <- factor(gse$Description, levels=c("RPmRNA","set5","set4","set3","set2","set1","ImmuneModu","TcellExaustion"))
      {
        p=ggplot(gse, aes(class, Description)) +
          geom_point(aes(color = p.adjust, size = Count)) +
          #scale_color_aaas()+
          #scale_x_continuous(limits=c(0.9,1.6))+
          scale_color_gradient(low=col[1],high=col[2],guide="colorbar")+
          theme_bw()+
          guides(fill=guide_legend(title=NULL, ncol=1))+
          #coord_fixed(ratio=0.3)+
          theme(
            plot.margin = unit(c(1, 1, 1, 1),"cm"),
            panel.border=element_rect(size=1, fill="transparent"),
            panel.grid=element_line(size=rel(0.5)),
            #legend.position=c(0.25,0.85),
            legend.position="right",
            legend.background = element_blank(),
            legend.title = element_text(color="black", size=10),
            legend.text= element_text(color="black", size=10),
            axis.text.x = element_text(color="black", size=10),
            axis.text.y = element_text(color="black", size=10),
            axis.title.x = element_text(color="black", size=10),
            axis.title.y = element_blank())+
          labs(title="")+xlab("")
        p
      }
      outdir="result/unpaired/enrich/"
      ggsave(paste0(outdir,"dotplot-cf-ev-GSEA-geneset",".pdf"), width = 3.5, height = 5)
      
    }
    # fgsea
    {
      rank = gene.list.ensgid
      geneset = list(CGC_1 = CGC_1$Ensg,
                     CGC_all = CGC_2$Ensg,
                     NCG_1 = NCG_1$ENSEMBL,
                     NCG_2 = NCG_2$ENSEMBL,
                     Uniprot = Uniprot$ENSEMBL,
                     Pansoft = Pansoft$ENSEMBL,
                     IntOGene = IntOGene$ENSEMBL)
      library(fgsea)
      fgsea <- fgsea(geneset, 
                     rank, 
                     nperm = 10,
                     minSize = 0, 
                     maxSize = 500)
      summary(fgsea)
      
      plotGseaTable(fgsea$pathway,
                    rank,
                    fgsea, 
                    gseaParam=0.5)
      plotEnrichment(fgsea$pathway,rank,)
      
    }
    {
      # emapplot
      enrichplot::emapplot(gse.GO)
      # cnetplot
      cnetplot(gsea)
      # 山峦图，展示每个geneset的基因logFC分布
      ridgeplot(gsea)
      # 选择单个gene set
      gsea3 <- data.frame(gsea)
      gseaplot2(gsea, geneSetID = 1, title=gsea3$Description[1])
      
    }
    
  }
}


########### Paired cf and EV in different cancer group 
{
  gencode <- read.delim2("count_matrix/gencode.cutoff_filter.24v24.txt", header = T, row.names = 1, check.names = F)
  dim(gencode)
  data = gencode
  ### paired group
  {
    # put cell-free as reference level
    cancer <- unlist(lapply(strsplit(colnames(data),"-",fixed=T),function(x) x[1]))
    cancer <- factor(cancer, levels = c("NC","CRC","LUAD"))
    # 24v24
    class <- factor(rep(c("cf","EV"),24), levels=c("cf","EV"))
    patient <- factor(rep(colnames(data)[seq(1,47,2)],each=2))
  }
  ### edgeR: tpm cutoff to filter genes
  {
    # lib.size to be specific as lib.size(hg38_dedup)
    # remove genes that are expressed less than 50% samples(tpm)
    {
      lib.size <- read.delim2("count_matrix/libsize.txt",sep="\t",header=T,row.names = 1)
      rownames(lib.size) <- annotate(rownames(lib.size))
      lib.size <- lib.size[colnames(data),]
      y <- DGEList(counts=data, lib.size=lib.size)
      #此时对lib.size进行normalized目的是通过tpm filter genes
      RPM <- round(edgeR::cpm(y, normalzied.lib.sizes = T),3)
      geneid <- rownames(RPM)
      split <- strsplit(geneid,"|",fixed=T)
      length <- unlist(lapply(split,function(split) split[2]))
      y$genes <- data.frame(length=as.numeric(length))
      RPKM <- round(rpkm(y),3)
      TPM <- round((RPKM/colSums(RPKM))*10^6,3)
      cpm <- as.data.frame(RPM)
      tpm <- as.data.frame(TPM)
      ###################
      ### group = class
      {
        # class
        cf = tpm[,class=="cf"]
        ev = tpm[,class=="EV"]
        n.cf = length(class[class=="cf"])/2
        n.ev = length(class[class=="EV"])/2
        
        detected_genes = data.frame(
          cellfree <- c(
            sum(rowSums(cf>10)>= n.cf),
            sum(rowSums(cf>2)>= n.cf),
            sum(rowSums(cf>5)>= n.cf)),
          EV <- c(
            sum(rowSums(ev>10)>= n.ev),
            sum(rowSums(ev>2)>= n.ev),
            sum(rowSums(ev>5)>= n.ev)))
        colnames(detected_genes) = c("cf","ev")
        rownames(detected_genes) = c("TPM>10(50%)","TPM>2(50%)","TPM>5(50%)")
        # stacked bar plot
        detected <- detected_genes %>% mutate(group=rownames(detected_genes))
        detected <- melt(detected)
        detected$group = factor(detected$group, levels=c("TPM>2(50%)","TPM>5(50%)","TPM>10(50%)"))
        {
          p=ggplot(detected, aes(x=variable, y=value, fill=group))+
            #geom_bar(stat="identity",width=0.7, col="black")+
            geom_bar(stat="identity",position="dodge",width=0.7, col="black")+
            #scale_y_continuous(limits=c(0,20000))+
            scale_fill_brewer()+
            theme_bw()+
            guides(fill=guide_legend(title=NULL, ncol=1))+
            theme(
              plot.margin = unit(c(1, 1, 0, 1),"cm"),
              panel.border=element_rect(size=1, fill="transparent"),
              panel.grid=element_blank(),
              #legend.position=c(0.25,0.85),
              legend.position="top",
              legend.background = element_blank(),
              legend.title = element_blank(),
              legend.text= element_text(face="bold", color="black", size=12),
              axis.text.x = element_text(face="bold",color="black", size=12),
              axis.text.y = element_text(face="bold",color="black", size=12),
              axis.title.x = element_text(face="bold",color="black", size=12),
              axis.title.y = element_text(face="bold",color="black", size=12))+
            ylab("Number of detected genes")+xlab("")
          #geom_vline(aes(xintercept=6.5))+
          p
        }             
        #ggsave("result/genenumber_cf-EV_TPMnormalized.pdf",p,width=4.5,height=4.5)               
        
        cf.genes <- rownames(cf)[rowSums(cf>5) >= n.cf]
        ev.genes <- rownames(ev)[rowSums(ev>5) >= n.ev]
        # venn plot
        {
          library(VennDiagram)
          color = brewer.pal(5, "Set3")
          venn <- venn.diagram(
            x = list(cf=cf.genes,
                     EV=ev.genes),
            filename = NULL,
            fill = color[1:2],
            alpha = 0.7,
            lwd = 2,
            lty = "blank",
            col = "transparent", #线条色
            cex = 1.5,
            fontfamily = "sans",
            fontface = "bold", #加粗
            cat.col = "black",
            cat.default.pos = "text", #位置, outer内 text外
            cat.cex = 1.5,
            cat.fontfamily = "sans",     
            cat.fontface = "bold",     
            cat.dist = c(0.1, 0.1), #位置，用圆的距离
            cat.pos = c(3,3), #位置，用圆的度数
            main = "",
            main.col = "black",
            main.cex = 1.5
            #main.fontface = "bold"
          );
          #grid.draw(venn);
          #ggsave("L1-L21/diff_exp_24v24/normalized/cf-ev_tpm>10_gene_number_venn.pdf",venn,width = 5,height = 5)
        }  
      }
      filter = (rowSums(cf>5) >= n.cf) | (rowSums(ev>5) >= n.ev)
      filter = (rowSums(cf>10) >= n.cf) | (rowSums(ev>10) >= n.ev)
      table(filter)
      y <- y[filter,]
      # tpm>10: 5114
      # tpm>5: 11478
    }
  }
  ### edgeR: Paired differential expression
  {
    # not normalized with lib.size
    y <- DGEList(counts=data, group=cancer, sample=patient)
    y <- y[filter,]   # filter: tpm>10: 5114 用5k个基因还是1w个基因差别不大
    dim(y)
    # 5114, 48
    # 11478, 48
    
    # TMM normalization
    EDASeq::plotRLE(edgeR::cpm(y)) 
    y <- calcNormFactors(y, method="TMM")
    EDASeq::plotRLE(edgeR::cpm(y))  
    
    limma::plotMDS(y)
    limma::plotMDS(edgeR::cpm(y, log=T), labels = class)
    limma::plotMDS(edgeR::cpm(y, log=T), labels = cancer)
    ########################################################
    ###### design: paired test in all samples
    cancer
    class
    patient
    
    difexp_cfEV_all <- function(){
      design <- model.matrix(~patient+class)
      # paired sample 和 EV cell-free class影响，不考虑癌症cancer之间的差异。
      
      y <- DGEList(counts=data, group=class, sample=patient)
      y <- y[filter,] 
      dim(y)
      y <- calcNormFactors(y, method="TMM")
      
      rownames(design) <- colnames(y)
      y <- estimateDisp(y, design)
      y$common.dispersion
      
      fit.ql <- glmQLFit(y, design)
      qlf <- glmQLFTest(fit.ql, coef = "classEV")
      de <- topTags(qlf, n=Inf)$table
      
      #edgeR::cpm(y)[rownames(de)[1:10],]
      #summary(decideTests(qlf))
      #plotMD(qlf)
      #abline(h=c(-1, 1), col="blue")
      
      de.cf <- filter(de, logFC <= -1, FDR <= 0.1)
      de.ev <- filter(de, logFC >= 1, FDR <= 0.1)
      print(paste0("cell-free enrichend RNA: ",nrow(de.cf)))
      print(paste0("EV enriched RNA: ",nrow(de.ev)))
      return(list(de=de, 
                  de.cf=de.cf, 
                  de.ev=de.ev))
    }
    
    de.all <- difexp_cfEV_all()
    nrow(de.all$de.ev)
    
    a=sort(de.all$de$logFC,decreasing = T)
    pdf("result/paired/all_barplot.pdf",width=5,height=3)
    barplot(a, axisnames=F, ylim=c(-4,4), ylab="log2FC",main="")
    dev.off()
    
    ########################################################
    ###### design: paired test in each cancer group
    cancer
    class
    patient
    
    test="NC"
    difexp_cfEV_each <- function(test){
      data.1 <- data[,cancer==test]
      cancer.1 <- unlist(lapply(strsplit(colnames(data.1),"-",fixed=T),function(x) x[1]))
      cancer.1 <- factor(cancer.1, levels = c(test))
      class.1 <- factor(rep(c("cf","EV"),8), levels=c("cf","EV"))
      patient.1 <- factor(rep(colnames(data.1)[seq(1,16,2)],each=2))
      
      design <- model.matrix(~patient.1+class.1)
      # patient放在第一项，关心class对一群患者的影响，class/treatment放在最后一项
      
      y <- DGEList(counts=data.1, group=class.1, sample=patient.1)
      y <- y[filter,] 
      dim(y)
      y <- calcNormFactors(y, method="TMM")
      #limma::plotMDS(edgeR::cpm(y, log=T), labels = class.1)
      
      rownames(design) <- colnames(y)
      y <- estimateDisp(y, design)
      y$common.dispersion
      plotBCV(y)
      
      fit.ql <- glmQLFit(y, design)
      qlf <- glmQLFTest(fit.ql)
      de <- topTags(qlf, n=Inf)$table
      
      #edgeR::cpm(y)[rownames(de)[1:10],]
      #summary(decideTests(qlf))
      #plotMD(qlf)
      #abline(h=c(-1, 1), col="blue")
      
      de.cf <- filter(de, logFC < -1, FDR < 0.1)
      de.ev <- filter(de, logFC > 1, FDR < 0.1)
      de.cf <- filter(de, logFC < -1, PValue < 0.05)
      de.ev <- filter(de, logFC > 1, PValue < 0.05)
      print(paste0("cell-free enrichend RNA: ",nrow(de.cf)))
      print(paste0("EV enriched RNA: ",nrow(de.ev)))
      return(list(de=de, 
                  de.cf=de.cf, 
                  de.ev=de.ev))
    }
    
    de.nc <- difexp_cfEV_each("NC")
    de.crc <- difexp_cfEV_each("CRC")
    de.luad <- difexp_cfEV_each("LUAD")
    
    a=sort(de.crc$de$logFC,decreasing = T)
    pdf("result/paired/barplot_crc.pdf",width=6,height=3)
    barplot(a, axisnames=F, ylim=c(-4,4), ylab="log2FC",main="CRC")
    dev.off()
    
    a=Reduce(intersect, list(rownames(de.nc$de.ev),
                             rownames(de.crc$de.ev),
                             rownames(de.luad$de.ev)))
    
    
  }
  ### correlation for all samples
  {
    # not normalized with lib.size
    y <- DGEList(counts=data, group=cancer, sample=patient)
    y <- y[filter,]   # filter: tpm>5: 11478 用5k个基因还是1w个基因差别不大
    dim(y)
    RPM <- round(edgeR::cpm(y),3)
    geneid <- rownames(RPM)
    split <- strsplit(geneid,"|",fixed=T)
    length <- unlist(lapply(split,function(split) split[2]))
    y$genes <- data.frame(length=as.numeric(length))
    RPKM <- round(rpkm(y),3)
    TPM <- round((RPKM/colSums(RPKM))*10^6,3)
    cpm <- as.data.frame(RPM)
    tpm <- as.data.frame(TPM)
    rpkm <- as.data.frame(RPKM)
    
    # 2 sample comparision
    { 
      a = "CRC-5"
      b = "CRC-5-EV"
      df <- cpm[cpm[,a]>0 & cpm[,b]>0, c(a,b)]  #CPM>0
      df <- tpm[tpm[,a]>0 & tpm[,b]>0, c(a,b)]  #TPM>0
      df <- rpkm[rpkm[,a]>0 & rpkm[,b]>0, c(a,b)] #RPKM>0
      
      df <- log10(df+1)
      # pearson correlation
      {
        cor.test(df[,a],
                 df[,b],
                 method = "pearson") # pearson比spearman算出的相关系数更高
        #label=paste("R = ",round(t$estimate,2),", p = ",signif(t$p.value,2),sep="")
        p=ggscatter(df,
                    x = a,
                    y = b,
                    #add = "reg.line", #拟合曲线
                    add.params = list(color=col[2],size=1,alpha=0.6),
                    conf.int = TRUE, #置信区间阴影带
                    cor.coef = TRUE, #系数
                    cor.method = "pearson", #方法
                    xlab = paste("log10(TPM+1), ", a, sep=""),
                    ylab = paste("log10(TPM+1), ", b, sep=""),
                    shape=16,size=2,color=col[2])
        #annotate("text",x=0.8,y=3.5,label="CPM>0",size=4,angle=0,color="black")
        #ggsave("",p,width=4,height=4)
      }
      # linear regression
      {
        fit <- lm(df[,b] ~ df[,a])
        summary <- summary(fit)
        rsq <- summary$adj.r.squared
      }
      # ggplot
      {
        ggplot(df, aes(x=`CRC-5`, y=`CRC-5-EV`))+
          geom_point(size = 2, shape = 16, color = col[2], alpha = 0.7)+
          geom_smooth(method = lm, size = 0.5, color = col[2])+
          #geom_abline(aes(intercept=0,slope=1),linetype=5,size=0.5, color=col[2])+
          #scale_x_continuous(breaks=c(0,3,6,9,12))+
          #scale_y_continuous(breaks=c(0,3,6,9,12))+
          theme_bw()+
          theme(
            #plot.margin = unit(c(1, 1, 1, 1),"cm"),
            legend.position=c(0.1,0.6),
            legend.title = element_text(color="black", size=12),
            legend.text= element_text(color="black", size=12),
            panel.grid = element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_text(color="black", size=12),
            axis.text.y = element_text(color="black", size=12),
            axis.title.x = element_text(color="black", size=12),
            axis.title.y = element_text(color="black", size=12),
            axis.line = element_line(size=0.5, colour = "black"))+
          labs(x="pico(logTPM)",y="ipico(logTPM)",face="bold")
        #annotate("text",x=3.5,y=12,label="R = 0.76, p < 2.2e-16",size=5,angle=0,color="black")
        annotate("text",x=1.3,y=6.5,label=bquote(R^2 ~ " = " ~.(rsq)),size=5,angle=0,color="black",fontface="bold")
        p
      }
    }
    
    # corrplot
    {
      # plot
      library(corrplot)
      library(psych)
      library(Hmisc)
      
      df <- tpm
      df <- cpm
      df <- rpkm
      #df <- log10(df+1)
      
      df <- df[,order(class, cancer)]
      class. <- class[order(class,cancer)]
      cancer. <- cancer[order(class,cancer)]
      
      # correlation
      rcorr <- rcorr(as.matrix(df), type="pearson")
      # R square
      {
        rsq=list()
        n=1
        for (i in seq(1,9,1)){
          for (j in seq(1,9,1)){
            fit = lm(df[,i] ~ df[,j])
            rsq[n] = summary(fit)$adj.r.squared
            n=n+1}}
        rsq <- matrix(unlist(rsq), nrow = 9, ncol = 9)
        colnames(rsq) <- c("H4R1","H4R2","H4R3","C1","C2","C3","N1","N2","N3")
        rownames(rsq) <- c("H4R1","H4R2","H4R3","C1","C2","C3","N1","N2","N3")
      }
      
      #col1 <- colorRampPalette(c("#6D9EC1", "white", "#E46726"))
      col1 <- colorRampPalette(c("blue", "white", "red"))
      
      a <- rcorr$r[c("CRC-15-EV","CRC-17-EV","CRC-2-EV","CRC-5-EV","CRC-8-EV",
                     "CRC-15","CRC-17","CRC-2","CRC-5","CRC-8"),
                   c("CRC-15-EV","CRC-17-EV","CRC-2-EV","CRC-5-EV","CRC-8-EV",
                     "CRC-15","CRC-17","CRC-2","CRC-5","CRC-8")]
      {
        pdf("reproducibility/corrplot/corrplot-rpkm.lib.size.pdf",width = 5,height = 5)
        par(oma=c(1,0,1,1),mar=c(2,0,2,2))
        # rcorr$r
        corrplot(a, 
                 method = "color",
                 type="lower",
                 #col=col[2],
                 is.corr=F,
                 order="original",
                 addrect = 4,
                 insig="p-value",
                 tl.pos="lt",
                 tl.col="black", 
                 tl.srt=90, 
                 cl.pos="b",
                 #addCoef.col="black",
                 addCoefasPercent=F) 
        
        #没有办法同时自定义颜色+自定义范围
        dev.off()
      }
      {
        pdf("reproducibility/corrplot/corrplot-rsq-rpkm.no.lib.size.pdf",width = 5,height = 5)
        par(oma=c(1,0,1,1),mar=c(2,0,2,2))
        # rsq
        corrplot(rsq, 
                 method = "color",
                 type="lower",
                 is.corr=F, 
                 order="original",
                 addrect = 4,
                 tl.pos="lt",
                 tl.col="black", 
                 tl.srt=90, 
                 cl.pos="b",
                 #cl.lim=c(0,1),
                 addCoef.col="black",
                 addCoefasPercent=F)  
        #没有办法同时自定义颜色+自定义范围
        dev.off()
      }
      {
        library(ggcorrplot)
        ggcorrplot(rcorr$r,
                   method = "square",
                   outline.col = "white",
                   type = "lower",
                   hc.order = F,
                   lab = T,
                   colors = c("#6D9EC1", "white", "#E46726"),
                   # pvalue
                   #p.mat = rcorr$P,
                   insig = "blank",
                   #ggtheme = ggplot2::theme_void,
                   #ggtheme = ggplot2::theme_bw(),
                   tl.col = "black",
                   tl.cex = 12)
        # ggcorrplot无法设置颜色范围？？
        ggsave("output/result/ggcorrplot.pdf",p,width=9,height=9)
      }
      ## pheatmap::pheatmap 
      {
        library(pheatmap)
        ann_col = data.frame(cancer = as.character(cancer.),
                             class = as.character(class.))
        rownames(ann_col) = colnames(rcorr$r)
        ann_row = data.frame(cancer = as.character(cancer.),
                             class = as.character(class.))
        rownames(ann_row) = rownames(rcorr$r)
        ann_colors = list(class = brewer.pal(5,"Set1")[4:5],
                          cancer = brewer.pal(5,"Set1")[1:3])
        names(ann_colors$class) = c("cf","EV")
        names(ann_colors$cancer) = c("NC","CRC","LUAD")
        
        col = brewer.pal(3,"Set1")
        pdf("result/paired/cf-EV-corr.pdf", width = 6, height = 5)
        par(mar=c(2,2,2,2))
        pheatmap(rcorr$r, 
                 color = colorRampPalette(c(col[2],"white",col[1]))(1000),
                 #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                 #cutree_col = 2, 
                 #cutree_row = 3, #break up the heatmap by clusters you define
                 cluster_rows = F, 
                 cluster_cols = F, #by default, pheatmap clusters by both row and col
                 show_rownames = F,
                 show_colnames = F,
                 fontsize_col = 12,
                 angle_col = 0,
                 annotation_col = ann_col,
                 annotation_row = ann_row,
                 annotation_colors = ann_colors,
                 annotation_names_row = F)
        dev.off()
      }
    }
    
    
    
    
    
  }
  ### PCA for all samples
  {
    library("preprocessCore")
    library("BBmisc")
    library("Rtsne")
    # normalize.quantiles: Using a normalization based upon quantiles, this function normalizes a matrix of probe level intensities.
    # normalize: Normalizes numeric data to a given scale.
    
    ###### X <- t(X)
    X <- t(edgeR::cpm(y))
    #X <- t(edgeR::cpm(y, log=T))
    
    # min_max_norm
    min_max_norm <- function(x) {
      (x - min(x)) / (max(x) - min(x))
    }
    X_norm.1 <- as.data.frame(apply(X, 2, min_max_norm))
    # min_max_norm
    X_norm.2 <- log2(normalize.quantiles(as.matrix(X), copy = FALSE)+0.001)   # preprocessCore
    X_norm.3 <- normalize(X_norm.2, method = "range", range=c(0,1), margin = 2L)   # BBmisc
    # min_max_norm (有问题？)
    X_norm.3 <- apply(X, 1, rescale)
    X_norm.3 <- t(X_norm.3)
    # z-score scaling
    X_norm.4 <- scale(X, center = T, scale = T)
    
    ## Visualization
    {  
      Y <- cancer
      colors <- rainbow(length(unique(Y)))
      colors_plot <- as.character(Y)
      colors_plot[which(colors_plot == "NC")] <- colors[1]
      colors_plot[which(colors_plot == "CRC")] <- colors[2]
      colors_plot[which(colors_plot == "LUAD")] <- colors[3]
    }
    {
      Y <- class
      colors <- rainbow(length(unique(Y)))
      colors_plot <- as.character(Y)
      colors_plot[which(colors_plot == "cf")] <- colors[1]
      colors_plot[which(colors_plot == "EV")] <- colors[2]
    }
    
    # PCA
    pca <- prcomp(X_norm.4, center = F, scale = F, rank. = 2)  
    pca.var <- pca$sdev^2  ## sdev是标准偏差，十个样本，就有十个标准偏差，平方是避免负数的干扰
    pca.var.per <- round(pca.var/sum(pca.var)*100, 1)  ##求每个样本的variation
    barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")  ##用柱状图可视化
    summary(pca)
    plot(pca$x, col = colors_plot, asp = 3, pch = 20, xlab = "component_1", ylab = "component_2", main = "PCA plot")
    ## ggplot
    pca.data <- data.frame(pca$x, class=class, cancer=cancer)
    percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
    percentage <- paste(colnames(pca.data)[1:9],"(", paste(as.character(percentage)), "%", ")", sep="")
    {
      p=ggplot(pca.data,aes(x=PC1,y=PC2,color=class,shape=cancer)) + 
        geom_point(size=3,alpha=1) +
        #scale_color_aaas() +
        scale_color_brewer(palette = "Set2") + 
        theme_bw()+
        theme(
          plot.margin = unit(c(1, 0, 1, 1),"cm"),
          panel.border=element_rect(size=1, fill="transparent"),
          panel.grid=element_blank(),
          legend.position = "right",
          legend.title = element_blank(),
          legend.text= element_text(face="bold", color="black", size=15),
          plot.title = element_text(hjust = 0.5, size=12, face="bold"),
          axis.text.x = element_text(face="bold", color="black", size=15),
          axis.text.y = element_text(face="bold",  color="black", size=15),
          axis.title.x = element_text(face="bold", color="black", size=15),
          axis.title.y = element_text(face="bold",color="black", size=15))+
        xlab(percentage[1]) + ylab(percentage[2])
      p
    }
    ggsave("result/paired/all_PCA.pdf",p,width = 5, height = 4.5)
    
    # tSNE
    tsne <- Rtsne(X_norm.4, dims=2, check_duplicates = FALSE, perplexity = 7) 
    # dims参数设置降维之后的维度，默认值为2
    # perplexity参数的取值必须小于(nrow(data) - 1 )/ 3
    # theta参数取值越大，结果的准确度越低，默认值为0.5
    # max_iter参数设置最大迭代次数。
    plot(tsne$Y, col = colors_plot, asp = 3, pch = 20,
         xlab = "tSNE_1", ylab = "tSNE_2", main = "tSNE plot")
    
  }
  ### volcano plot
  {
    de.plt <- de.all$de
    de.plt$threshold <- factor(ifelse(de.plt$FDR < 0.1 & abs(de.plt$logFC) >=1, ifelse(de.plt$logFC > 1 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
    
    de.nc.plt <- de.nc$de
    de.nc.plt$threshold <- factor(ifelse(de.nc.plt$PValue < 0.05 & abs(de.nc.plt$logFC) >=1, ifelse(de.nc.plt$logFC > 1 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
    de.crc.plt <- de.crc$de
    de.crc.plt$threshold <- factor(ifelse(de.crc.plt$PValue < 0.05 & abs(de.crc.plt$logFC) >=1, ifelse(de.crc.plt$logFC > 1 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
    de.luad.plt <- de.luad$de
    de.luad.plt$threshold <- factor(ifelse(de.luad.plt$PValue < 0.05 & abs(de.luad.plt$logFC) >=1, ifelse(de.luad.plt$logFC > 1 ,'Up','Down'),'Not'), levels = c("Down","Up","Not"))
    
    col = brewer.pal(3,"Set1")
    {
      p=ggplot(de.plt, aes(x=logFC, y =-log10(FDR), color=threshold)) + 
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
        labs(x="log2FoldChange",y="-log10 (adjusted p-value)",title="", face="bold")
      p
    }
    {
      p=ggplot(de.luad.plt, aes(x=logFC, y =-log10(PValue), color=threshold)) + 
        scale_color_manual(values=c(col[2],col[1],"grey"))+
        geom_point(size=2, alpha=0.8, shape=16, stroke = 0)+
        #scale_y_continuous(limits=c(0,2.8))+
        theme_bw(base_size = 12) +
        geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+
        geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
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
        labs(x="log2FoldChange",y="-log10 (p-value)",title="", face="bold")
      p
    }
    ggsave("result/paired/difexp_luad_volcano.pdf",p,width = 4.5,height = 4.5)
  }
  ### venn
  {
    color = brewer.pal(6, "Set3")
    venn <- venn.diagram(
      x = list(      
        #NC.cf=rownames(de.nc$de.cf),
        #CRC.cf=rownames(de.crc$de.cf),
        #LUAD.cf=rownames(de.luad$de.cf)
        NC.EV=rownames(de.nc$de.ev),
        CRC.EV=rownames(de.crc$de.ev),
        LUAD.EV=rownames(de.luad$de.ev)
      ),
      filename = NULL,
      fill = color[1:3],
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
      cat.dist = c(0.15, 0.15, 0.15), 
      cat.pos = c(330,30,180),
      main = "",
      main.col = "black",
      main.cex = 1.5
      #main.fontface = "bold"
    );
    grid.draw(venn);
    #ggsave("result/paired/EV_genes_venn.pdf",venn,width = 5,height = 5)
  } 
  ### heatmap
  {
    difgenes <- c(rownames(de.all$de.cf),rownames(de.all$de.ev))
    difgenes <- c(rownames(de.nc$de.cf),rownames(de.nc$de.ev))
    difgenes <- c(rownames(de.crc$de.cf),rownames(de.crc$de.ev))
    difgenes <- c(rownames(de.luad$de.cf),rownames(de.luad$de.ev))
    
    logRPM <- edgeR::cpm(y, log=T)
    logRPM <- logRPM[difgenes,]
    logRPM.scale <- scale(t(logRPM), center = T, scale = T)
    logRPM.scale <- t(logRPM.scale)[,c(seq(1,47,2),seq(2,48,2))]
    
    cancer. <- rep(rep(c("CRC","LUAD","NC"),each=8),2)
    class. <- rep(c("cf","EV"),each=24)
    
    ## pheatmap::pheatmap 
    {
      library(pheatmap)
      ann_col = data.frame(cancer. = as.character(cancer.),
                           class. = as.character(class.))
      rownames(ann_col) = colnames(logRPM.scale)
      #ann_row = data.frame(GeneClass = factor(rep(c("CRC.cf","LUAD.cf","NC.EV", "CRC.EV", "NC.CRC.EV","LUAD.EV"), c(5,148,8,0,1,5))))
      #rownames(ann_row) = difgenes
      #set colors of each group
      ann_colors = list(cancer. = brewer.pal(3,"Set1")[1:3], 
                        class. = brewer.pal(3,"Set2")[1:2])
      #GeneClass = c(Path1 = "#807DBA", Path2 = "#9E9AC8", Path3 = "#BCBDDC"))
      
      col = brewer.pal(3,"Set1")
      pdf("result/paired/pheatmap_ncgenes.pdf", width = 6, height = 7)
      par(mar=c(2,2,2,2))
      pheatmap(logRPM.scale, 
               color = colorRampPalette(c(col[2],"white",col[1]))(1000),
               #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
               cutree_col = 2, 
               #cutree_row = 3, #break up the heatmap by clusters you define
               cluster_rows = F, cluster_cols = F, #by default, pheatmap clusters by both row and col
               show_rownames = F,
               show_colnames = F,
               fontsize_col = 12,
               angle_col = 0,
               annotation_col = ann_col,
               #annotation_row = ann_row,
               #annotation_colors = ann_colors,
               annotation_names_row = F)
      dev.off()
    }
    
  }
  ### enrichment analysis
  {
    library(clusterProfiler)
    library(enrichplot)
    {
      ref <- read.table("/BioII/lulab_b/zhanqing/cfRNAseq/reference/gtf/genome/gene.gtf",sep = "\t",header = F)
      colnames(ref) <- c("chr","source","region","start","end","score","strand","phase","meta")
      a <- unlist(lapply(strsplit(ref[,9],".",fixed=T), function(x) x[1]))
      ref$ensg <- unlist(lapply(strsplit(a," ",fixed=T), function(x) x[2]))
      library(gprofiler2)
      tmp <- gprofiler2::gconvert(query = ref$ensg, organism = "hsapiens", target = "ENTREZGENE_ACC" )
      tmp <- tmp[,c("input","target")]
      colnames(tmp) <- c("ensg","ENTREZID")
      ref$ENTREZID <- tmp$ENTREZID[match(ref$ensg,tmp$ensg)]
      ref <- ref[!duplicated(ref$ensg),]
      table(duplicated(ref$ENTREZID))
      #write.table(ref,"/BioII/lulab_b/zhanqing/cfRNAseq/reference/geneset/gencode.v38.txt",quote = F,sep = "\t",col.names = T,row.names = F)
    }
    #### input table
    #de.order <- de.all$de
    de.order <- de.nc$de
    #de.order <- de.crc$de
    #de.order <- de.luad$de
    
    de.order$ENSG <- rownames(de.order)
    de.order$ENSG <- as.character(lapply(strsplit(de.order$ENSG,"|",fixed=TRUE), function(x) x[1]))
    de.order$ENSG <- as.character(lapply(strsplit(de.order$ENSG,".",fixed=TRUE), function(x) x[1]))
    de.order <- de.order[!duplicated(de.order$ENSG),]
    
    ##### get all geneList (for GSEA)
    #gene.list.ncbiid <- bitr(de.order$name, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    gene.list.ncbiid <- read.table("/BioII/lulab_b/zhanqing/cfRNAseq/reference/geneset/gencode.v38.txt",sep = "\t",header = T)
    gene.list.ncbiid <- gene.list.ncbiid[!is.na(gene.list.ncbiid$ENTREZID),]
    gene.list.ncbiid <- gene.list.ncbiid[,c("ensg","ENTREZID")]
    gene.list.ncbiid <- gene.list.ncbiid[!duplicated(gene.list.ncbiid$ensg) & !duplicated(gene.list.ncbiid$ENTREZID),]
    colnames(gene.list.ncbiid)[1] <- "ENSG"
    
    #de.order <- merge(x = de.order, y = gene.list.ncbiid, by.x="name", by.y="ENSG")
    de.order <- dplyr::left_join(x = de.order, y = gene.list.ncbiid)
    de.order <- de.order[order(de.order$logFC, -de.order$PValue, decreasing = T),] # order by logFC then PValue
    
    gene.list.ncbiid <- de.order$logFC[!duplicated(de.order$ENTREZID)]
    names(gene.list.ncbiid) <- de.order$ENTREZID[!duplicated(de.order$ENTREZID)]
    message(paste0("remaining ENTREZID geneList lenth: ",length(gene.list.ncbiid)))
    gene.list.ensgid <- de.order$logFC
    names(gene.list.ensgid) <- de.order$ENSG
    message(paste0("remaining ENSEMBL geneList lenth: ",length(gene.list.ensgid)))
    
    
    ##### get top geneList (for ORA)
    method="cutoff.FDR" # for all samples
    method="cutoff.pvalue" # for each cancer group
    {
      if (method=="cutoff.FDR"){
        cutoff.FDR = 0.1
        cutoff.logFC = 1
        message(paste0("FDR cutoff: ", cutoff.FDR, "  logFC cutoff: ", cutoff.logFC))
        gene.top.ncbiid <- as.character(na.omit(de.order$ENTREZID[abs(de.order$logFC) > cutoff.logFC & de.order$FDR <= cutoff.FDR]))   # FDR cutoff
        gene.top.ncbiid.up <- as.character(na.omit(de.order$ENTREZID[de.order$logFC > cutoff.logFC & de.order$FDR <= cutoff.FDR]))
        gene.top.ncbiid.down <- as.character(na.omit(de.order$ENTREZID[de.order$logFC < -cutoff.logFC & de.order$FDR <= cutoff.FDR]))
        message(paste0("remaining ENTREZID top gene lenth: ",length(gene.top.ncbiid)))
        message(paste0("remaining ENTREZID top up gene lenth: ",length(gene.top.ncbiid.up)))
        message(paste0("remaining ENTREZID top down gene lenth: ",length(gene.top.ncbiid.down)))
        gene.top.ensgid <- de.order$ENSG[abs(de.order$logFC) > cutoff.logFC & de.order$FDR <= cutoff.FDR]   # FDR cutoff
        gene.top.ensgid.up <- de.order$ENSG[de.order$logFC > cutoff.logFC & de.order$FDR <= cutoff.FDR] 
        gene.top.ensgid.down <- de.order$ENSG[de.order$logFC < -cutoff.logFC & de.order$FDR <= cutoff.FDR]
        message(paste0("remaining ENSEMBL top gene lenth: ",length(gene.top.ensgid)))
        message(paste0("remaining ENSEMBL top up gene lenth: ",length(gene.top.ensgid.up)))
        message(paste0("remaining ENSEMBL top down gene lenth: ",length(gene.top.ensgid.down)))
      }else if (method=="cutoff.pvalue"){
        cutoff.pvalue = 0.05
        cutoff.logFC = 1
        message(paste0("p-value cutoff: ", cutoff.pvalue, "  logFC cutoff: ", cutoff.logFC))
        gene.top.ncbiid <- as.character(na.omit(de.order$ENTREZID[abs(de.order$logFC) > cutoff.logFC & de.order$PValue <= cutoff.pvalue]))   # p-value cutoff
        gene.top.ncbiid.up <- as.character(na.omit(de.order$ENTREZID[de.order$logFC > cutoff.logFC & de.order$PValue <= cutoff.pvalue])) 
        gene.top.ncbiid.down <- as.character(na.omit(de.order$ENTREZID[de.order$logFC < -cutoff.logFC & de.order$PValue <= cutoff.pvalue]))
        message(paste0("remaining ENTREZID top gene lenth: ",length(gene.top.ncbiid)))
        message(paste0("remaining ENTREZID top up gene lenth: ",length(gene.top.ncbiid.up)))
        message(paste0("remaining ENTREZID top down gene lenth: ",length(gene.top.ncbiid.down)))
        gene.top.ensgid <- de.order$ENSG[abs(de.order$logFC) > cutoff.logFC & de.order$PValue <= cutoff.pvalue]   # FDR cutoff
        gene.top.ensgid.up <- de.order$ENSG[de.order$logFC > cutoff.logFC & de.order$PValue <= cutoff.pvalue] 
        gene.top.ensgid.down <- de.order$ENSG[de.order$logFC < -cutoff.logFC & de.order$PValue <= cutoff.pvalue]
        message(paste0("remaining ENSEMBL top gene lenth: ",length(gene.top.ensgid)))
        message(paste0("remaining ENSEMBL top up gene lenth: ",length(gene.top.ensgid.up)))
        message(paste0("remaining ENSEMBL top down gene lenth: ",length(gene.top.ensgid.down)))
      }else if (method=="FDR.top"){
        FDR.top = 100
        message(paste0("FDR top: ", FDR.top))
        de.order <- de.order[order(-de.order$FDR, decreasing = T),]  # order by FDR
        gene.top.ncbiid <- as.character(na.omit(de.order$ENTREZID[1:FDR.top]))
        gene.top.ncbiid.up <- as.character(na.omit(de.order$ENTREZID[de.order$logFC > 0][1:FDR.top]))
        gene.top.ncbiid.down <- as.character(na.omit(de.order$ENTREZID[de.order$logFC < 0][1:FDR.top]))
        message(paste0("remaining ENTREZID top gene lenth: ",length(gene.top.ncbiid)))
        message(paste0("remaining ENTREZID top up gene lenth: ",length(gene.top.ncbiid.up)))
        message(paste0("remaining ENTREZID top down gene lenth: ",length(gene.top.ncbiid.down)))
        gene.top.ensgid <- de.order$ENSG[1:FDR.top]   
        gene.top.ensgid.up <- de.order$ENSG[1:FDR.top]
        gene.top.ensgid.down <- de.order$ENSG[1:FDR.top]
        message(paste0("remaining ENSEMBL top gene lenth: ",length(gene.top.ensgid)))
        message(paste0("remaining ENSEMBL top up gene lenth: ",length(gene.top.ensgid.up)))
        message(paste0("remaining ENSEMBL top down gene lenth: ",length(gene.top.ensgid.down)))
      }else{
        message(paste0("No methods specified"))
      }
    }
    
    pValue=0.1
    qValue=0.1
    pAdjust="BH"
    # 1. Over-Representation Analysis 超几何检验或Fisher精确检验
    message("start ORA enrich...")
    # 1.1 ORA all gene
    {
      if(FALSE){
        enrich.GO <- enrichGO(
          gene = gene.top.ensgid,  # gene.top.ncbiid
          #universe = names(gene.list.ncbiid),  #选择background gene?
          OrgDb = org.Hs.eg.db,
          ont = "ALL",
          keyType = "ENSEMBL",  # ENTREZID
          pvalueCutoff = pValue,
          pAdjustMethod = pAdjust,
          qvalueCutoff = qValue)
      }
      enrich.GO <- enrichGO(
        gene = gene.top.ncbiid,
        #universe = names(gene.list.ncbiid),  #选择background gene?
        OrgDb = org.Hs.eg.db,
        ont = "ALL",
        keyType = "ENTREZID",
        pvalueCutoff = pValue,
        pAdjustMethod = pAdjust,
        qvalueCutoff = qValue
      )
      enrich.DO <- enrichDO(
        gene = gene.top.ncbiid, 
        ont = "DO",
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust,
        qvalueCutoff = 1,
        #universe = names(gene.list.ncbiid),
        minGSSize = 5,
        maxGSSize = 500,
        readable = FALSE
      )
      enrich.KEGG <- enrichKEGG(
        gene = gene.top.ncbiid,
        keyType = "ncbi-geneid",
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust,
        qvalueCutoff = 1
      )
    }
    # 1.2 ORA up gene
    {
      if(FALSE){
        enrich.GO.up <- enrichGO(
          gene = gene.top.ncbiid.up,
          #universe = names(gene.list.ncbiid),  #选择background gene?
          OrgDb = org.Hs.eg.db,
          ont = "ALL",
          keyType = "ENTREZID",
          pvalueCutoff = 1,
          pAdjustMethod = pAdjust,
          qvalueCutoff = 1)
      }
      enrich.GO.up <- enrichGO(
        gene = gene.top.ensgid.up,
        #universe = names(gene.list.ncbiid),  #选择background gene?
        OrgDb = org.Hs.eg.db,
        ont = "ALL",
        keyType = "ENSEMBL",
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust,
        qvalueCutoff = 1
      )
      enrich.DO.up <- enrichDO(
        gene = gene.top.ncbiid.up,
        ont = "DO",
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust,
        qvalueCutoff = 1,
        #universe = names(gene.list.ncbiid),
        minGSSize = 5,
        maxGSSize = 500,
        readable = FALSE
      )
      enrich.KEGG.up <- enrichKEGG(
        gene = gene.top.ncbiid.up,
        keyType = "ncbi-geneid",
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust,
        qvalueCutoff = 1
      )
    }
    # 1.3 ORA down gene
    {
      if(FALSE){
        enrich.GO.down <- enrichGO(
          gene = gene.top.ncbiid.down,
          #universe = names(gene.list.ncbiid),  #选择background gene?
          OrgDb = org.Hs.eg.db,
          ont = "ALL",
          keyType = "ENTREZID",
          pvalueCutoff = pValue,
          pAdjustMethod = pAdjust,
          qvalueCutoff = qValue)
      }
      enrich.GO.down <- enrichGO(
        gene = gene.top.ensgid.down,
        #universe = names(gene.list.ncbiid),  #选择background gene?
        OrgDb = org.Hs.eg.db,
        ont = "ALL",
        keyType = "ENSEMBL",
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust,
        qvalueCutoff = 1
      )
      enrich.DO.down <- enrichDO(
        gene = gene.top.ncbiid.down,
        ont = "DO",
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust,
        qvalueCutoff = 1,
        #universe = names(gene.list.ncbiid),
        minGSSize = 5,
        maxGSSize = 500,
        readable = FALSE
      )
      enrich.KEGG.down <- enrichKEGG(
        gene = gene.top.ncbiid.down,
        keyType = "ncbi-geneid",
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust,
        qvalueCutoff = 1
      )
    }
    # output table
    {
      types <- c(
        enrich.GO.up,enrich.DO.up,enrich.KEGG.up,
        enrich.GO.down,enrich.DO.down,enrich.KEGG.down)
      names(types) <-c(
        "enrich.GO.up","enrich.DO.up","enrich.KEGG.up",
        "enrich.GO.down","enrich.DO.down","enrich.KEGG.down")
      
      outdir <- paste0("result/paired/enrich/nc-table-5114genes/")
      dir.create(paste0(outdir), showWarnings = T)
      message("start writing tables...")
      for(i in 1:length(types)){
        write.table(types[[i]],paste0(outdir,"/", names(types)[i], ".txt"),row.names = F,col.names = T, quote = F,sep = "\t")}
    }
    # enrichment plot of all samples
    {
      options(scipen=999) 
      outdir="result/paired/enrich/nc-table-5114genes/"
      col = brewer.pal(3,"Set1")
      pdf(paste0(outdir,"dotplot-enrich.GO.down",".pdf"), width = 6, height = 5)
      enrichplot::dotplot(enrich.GO.down, font.size=10, title="enrich.GO.down", showCategory = 20) +
        scale_color_continuous(low=col[1],high=col[2]) +
        scale_x_continuous(limits = c(0.01,0.09))
      dev.off()
      pdf(paste0(outdir,"dotplot-enrich.GO.up",".pdf"), width = 6.5, height = 5)
      enrichplot::dotplot(enrich.GO.up, font.size=8, title="enrich.GO.up", showCategory = 20) +
        scale_color_continuous(low=col[1],high=col[2]) +
        scale_x_continuous(limits = c(0.034,0.109))
      dev.off()
      
      pdf(paste0(outdir,"dotplot-enrich.KEGG.up",".pdf"), width = 5, height = 4)
      enrichplot::dotplot(enrich.KEGG.up, font.size=12, title="enrich.KEGG.up", showCategory = 10) +
        scale_color_continuous(low=col[1],high=col[2]) +
        scale_x_continuous(limits = c(0.06,0.15), breaks=c(0.06,0.1,0.14))
      dev.off()
      pdf(paste0(outdir,"dotplot-enrich.KEGG.down",".pdf"), width = 5.3, height = 4)
      enrichplot::dotplot(enrich.KEGG.down, font.size=12, title="enrich.KEGG.down", showCategory = 10) +
        scale_color_continuous(low=col[1],high=col[2]) +
        scale_x_continuous(limits = c(0.015,0.13), breaks = c(0.05,0.1))
      dev.off()
      
      go.up <- read.table("result/paired/enrich/nc-table-5114genes/enrich.GO.down.txt",sep="\t",header=T,quote="")
      go.up <- go.up[order(go.up$p.adjust, decreasing = F),][1:20,]
      go.up$Description <- factor(go.up$Description, levels=go.up$Description[order(go.up$ONTOLOGY,-go.up$p.adjust)])
      go.up$logp <- -log10(go.up$p.adjust)
      
      
      
      {
        p=ggplot(go.up, aes(logp, Description)) +
          geom_point(aes(color = ONTOLOGY, size = Count)) +
          scale_color_aaas()+
          scale_x_continuous(limits=c(0.3,0.72))+
          #scale_color_gradient(low=col[1],high=col[2],guide="colorbar")+
          theme_bw()+
          guides(fill=guide_legend(title=NULL, ncol=1))+
          #coord_fixed(ratio=0.5)+
          theme(
            plot.margin = unit(c(1, 1, 1, 1),"cm"),
            panel.border=element_rect(size=1, fill="transparent"),
            panel.grid=element_line(size=rel(0.5)),
            #legend.position=c(0.25,0.85),
            legend.position="right",
            legend.background = element_blank(),
            legend.title = element_text(color="black", size=10),
            legend.text= element_text(color="black", size=10),
            axis.text.x = element_text(color="black", size=10),
            axis.text.y = element_text(color="black", size=8),
            axis.title.x = element_text(color="black", size=10),
            axis.title.y = element_blank())+
          labs(title="GO.down")+
          xlab("p.adjust")
        p
        ggsave(paste0(outdir,"dotplot-enrich.GO.down-2",".pdf"), width = 6.5, height = 5)
      }
      
      
      
      # DAG 有向无环图
      plotGOgraph(enrich.GO)
      # igraph布局的DAG
      goplot(enrich.GO)
      
      barplot(enrich.GO, showCategory = 10)
      barplot(enrich.DO, showCategory = 10)
      barplot(enrich.KEGG, showCategory = 10)
      
    }
    # enrichment plot of each cancer samples
    {
      options(scipen=-1) 
      nc.go.up <- read.table("result/enrich/nc-table/enrich.GO.up.txt",sep="\t",header=T,quote="")
      nc.kegg.up <- read.table("result/enrich/nc-table/enrich.KEGG.up.txt",sep="\t",header=T,quote="")
      nc.do.up <- read.table("result/enrich/nc-table/enrich.DO.up.txt",sep="\t",header=T,quote="")
      nc.go.down <- read.table("result/enrich/nc-table/enrich.GO.down.txt",sep="\t",header=T,quote="")
      nc.kegg.down <- read.table("result/enrich/nc-table/enrich.KEGG.down.txt",sep="\t",header=T,quote="")
      nc.do.down <- read.table("result/enrich/nc-table/enrich.DO.down.txt",sep="\t",header=T,quote="")
      
      crc.go.up <- read.table("result/enrich/crc-table/enrich.GO.up.txt",sep="\t",header=T,quote="")
      crc.kegg.up <- read.table("result/enrich/crc-table/enrich.KEGG.up.txt",sep="\t",header=T,quote="")
      crc.do.up <- read.table("result/enrich/crc-table/enrich.DO.up.txt",sep="\t",header=T,quote="")
      crc.go.down <- read.table("result/enrich/crc-table/enrich.GO.down.txt",sep="\t",header=T,quote="")
      crc.kegg.down <- read.table("result/enrich/crc-table/enrich.KEGG.down.txt",sep="\t",header=T,quote="")
      crc.do.down <- read.table("result/enrich/crc-table/enrich.DO.down.txt",sep="\t",header=T,quote="")
      
      luad.go.up <- read.table("result/enrich/luad-table/enrich.GO.up.txt",sep="\t",header=T,quote="")
      luad.kegg.up <- read.table("result/enrich/luad-table/enrich.KEGG.up.txt",sep="\t",header=T,quote="")
      luad.do.up <- read.table("result/enrich/luad-table/enrich.DO.up.txt",sep="\t",header=T,quote="")
      luad.go.down <- read.table("result/enrich/luad-table/enrich.GO.down.txt",sep="\t",header=T,quote="")
      luad.kegg.down <- read.table("result/enrich/luad-table/enrich.KEGG.down.txt",sep="\t",header=T,quote="")
      luad.do.down <- read.table("result/enrich/luad-table/enrich.DO.down.txt",sep="\t",header=T,quote="")
      
      id.go.up <- Reduce(union,list(nc.go.up[order(nc.go.up$p.adjust, decreasing = F),"ID"][1:10],
                                    crc.go.up[order(crc.go.up$p.adjust, decreasing = F),"ID"][1:10],
                                    luad.go.up[order(luad.go.up$p.adjust, decreasing = F),"ID"][1:5]))
      df <- rbind(nc.go.up[nc.go.up$ID %in% id.go.up,] %>% mutate(group = "NC"),
                  crc.go.up[crc.go.up$ID %in% id.go.up,] %>% mutate(group = "CRC"),
                  luad.go.up[luad.go.up$ID %in% id.go.up,] %>% mutate(group = "LUAD"))
      df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
      df$Description <- factor(df$Description, levels= rev(df$Description[match(id.go.up, df$ID)]))
      
      id.go.down <- Reduce(union,list(nc.go.down[order(nc.go.down$p.adjust, decreasing = F),"ID"][1:3],
                                      crc.go.down[order(crc.go.down$p.adjust, decreasing = F),"ID"][1:10],
                                      luad.go.down[order(luad.go.down$p.adjust, decreasing = F),"ID"][1:10]))
      df <- rbind(nc.go.down[nc.go.down$ID %in% id.go.down,] %>% mutate(group = "NC"),
                  crc.go.down[crc.go.down$ID %in% id.go.down,] %>% mutate(group = "CRC"),
                  luad.go.down[luad.go.down$ID %in% id.go.down,] %>% mutate(group = "LUAD"))
      df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
      df$Description <- factor(df$Description, levels= rev(df$Description[match(id.go.down, df$ID)]))
      
      id.kegg.up <- Reduce(union,list(nc.kegg.up[order(nc.kegg.up$p.adjust, decreasing = F),"ID"][1:10]))
      #crc.kegg.up[order(crc.kegg.up$p.adjust, decreasing = F),"ID"][1:10],
      #luad.kegg.up[order(luad.kegg.up$p.adjust, decreasing = F),"ID"][1:10]))
      df <- rbind(nc.kegg.up[nc.kegg.up$ID %in% id.kegg.up,] %>% mutate(group = "NC"),
                  crc.kegg.up[crc.kegg.up$ID %in% id.kegg.up,] %>% mutate(group = "CRC"),
                  luad.kegg.up[luad.kegg.up$ID %in% id.kegg.up,] %>% mutate(group = "LUAD"))
      df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
      df$Description <- factor(df$Description, levels= rev(df$Description[match(id.kegg.up, df$ID)]))
      
      id.kegg.down <- Reduce(union,list(#nc.kegg.down[order(nc.kegg.down$p.adjust, decreasing = F),"ID"][1:10],
        #crc.kegg.down[order(crc.kegg.down$p.adjust, decreasing = F),"ID"][1:10],
        luad.kegg.down[order(luad.kegg.down$p.adjust, decreasing = F),"ID"][1:10]))
      df <- rbind(nc.kegg.down[nc.kegg.down$ID %in% id.kegg.down,] %>% mutate(group = "NC"),
                  crc.kegg.down[crc.kegg.down$ID %in% id.kegg.down,] %>% mutate(group = "CRC"),
                  luad.kegg.down[luad.kegg.down$ID %in% id.kegg.down,] %>% mutate(group = "LUAD"))
      df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
      df$Description <- factor(df$Description, levels= rev(df$Description[match(id.kegg.down, df$ID)]))
      
      id.do.up <- Reduce(union,list(nc.do.up[order(nc.do.up$p.adjust, decreasing = F),"ID"][1:10],
                                    crc.do.up[order(crc.do.up$p.adjust, decreasing = F),"ID"][1:10],
                                    luad.do.up[order(luad.do.up$p.adjust, decreasing = F),"ID"][1:10]))
      df <- rbind(nc.do.up[nc.do.up$ID %in% id.do.up,] %>% mutate(group = "NC"),
                  crc.do.up[crc.do.up$ID %in% id.do.up,] %>% mutate(group = "CRC"),
                  luad.do.up[luad.do.up$ID %in% id.do.up,] %>% mutate(group = "LUAD"))
      df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
      df$Description <- factor(df$Description, levels= rev(df$Description[match(id.do.up, df$ID)]))
      
      id.do.down <- Reduce(union,list(nc.do.down[order(nc.do.down$p.adjust, decreasing = F),"ID"][1:10],
                                      crc.do.down[order(crc.do.down$p.adjust, decreasing = F),"ID"][1:10],
                                      luad.do.down[order(luad.do.down$p.adjust, decreasing = F),"ID"][1:10]))
      df <- rbind(nc.do.down[nc.do.down$ID %in% id.do.down,] %>% mutate(group = "NC"),
                  crc.do.down[crc.do.down$ID %in% id.do.down,] %>% mutate(group = "CRC"),
                  luad.do.down[luad.do.down$ID %in% id.do.down,] %>% mutate(group = "LUAD"))
      df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
      df$Description <- factor(df$Description, levels= rev(df$Description[match(id.do.down, df$ID)]))
      
      {
        kegg.up.1 <- nc.kegg.up[order(nc.kegg.up$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "NC")
        kegg.up.2 <- crc.kegg.up[order(crc.kegg.up$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "CRC")
        kegg.up.3 <- luad.kegg.up[order(luad.kegg.up$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "LUAD")
        kegg.down.1 <- nc.kegg.down[order(nc.kegg.down$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "NC")
        kegg.down.2 <- crc.kegg.down[order(crc.kegg.down$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "CRC")
        kegg.down.3 <- luad.kegg.down[order(luad.kegg.down$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "LUAD")
        df <- rbind(kegg.up.1, kegg.up.2, kegg.up.3)
        df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
        des <- df$Description[!duplicated(df$Description)]
        df$Description <- factor(df$Description, levels=rev(des))
        df <- rbind(kegg.down.1, kegg.down.2, kegg.down.3)
        df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
        des <- df$Description[!duplicated(df$Description)]
        df$Description <- factor(df$Description, levels=rev(des))
        
        do.up.1 <- nc.do.up[order(nc.do.up$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "NC")
        do.up.2 <- crc.do.up[order(crc.do.up$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "CRC")
        do.up.3 <- luad.do.up[order(luad.do.up$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "LUAD")
        do.down.1 <- nc.do.down[order(nc.do.down$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "NC")
        do.down.2 <- crc.do.down[order(crc.do.down$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "CRC")
        do.down.3 <- luad.do.down[order(luad.do.down$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "LUAD")
        df <- rbind(do.up.1, do.up.2, do.up.3)
        df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
        des <- df$Description[!duplicated(df$Description)]
        df$Description <- factor(df$Description, levels=rev(des))
        df <- rbind(do.down.1, do.down.2, do.down.3)
        df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
        des <- df$Description[!duplicated(df$Description)]
        df$Description <- factor(df$Description, levels=rev(des))
        
      }
      
      col = brewer.pal(8,"Set1")
      {
        p=ggplot(df, aes(group, Description), showCategory = 20) +
          geom_point(aes(color = p.adjust, size = Count)) +
          scale_color_gradient(low=col[1],high=col[2],guide="colorbar")+
          theme_bw()+
          guides(fill=guide_legend(title=NULL, ncol=1))+
          coord_fixed(ratio=0.5)+
          theme(
            plot.margin = unit(c(1, 1, 1, 1),"cm"),
            panel.border=element_rect(size=1, fill="transparent"),
            panel.grid=element_line(size=rel(0.5)),
            #legend.position=c(0.25,0.85),
            legend.position="right",
            legend.background = element_blank(),
            legend.title = element_text(color="black", size=12),
            legend.text= element_text(color="black", size=12),
            axis.text.x = element_text(color="black", size=8),
            axis.text.y = element_text(color="black", size=10),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())+
          labs(title="DO.up")
        p
      }
      #ggsave("result/enrich/cancer-plot/DO.up.pdf",p,width=5.5,height=5.5)
      
      
      # DAG 有向无环图
      plotGOgraph(enrich.GO)
      # igraph布局的DAG
      goplot(enrich.GO)
      
      barplot(enrich.GO, showCategory = 10)
      barplot(enrich.DO, showCategory = 10)
      barplot(enrich.KEGG, showCategory = 10)
      
    }
    # 2.Gene Set Enrichment Analysis online
    message("start gsea...")
    {
      gse.GO <- gseGO(
        geneList = gene.list.ensgid, 
        OrgDb = org.Hs.eg.db,
        ont = "ALL", 
        keyType = "ENSEMBL",  
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust
        #by = "fgsea", #'DOSE'
      ) 
      if(FALSE){
        gse.GO <- gseGO(
          geneList = gene.list.ncbiid, 
          OrgDb = org.Hs.eg.db,
          ont = "ALL", 
          keyType = "ENTREZID",  
          pvalueCutoff = pValue,
          pAdjustMethod = pAdjust
          #by = "fgsea", #'DOSE'
        ) 
      }
      gse.KEGG <- gseKEGG(
        geneList = gene.list.ncbiid,
        keyType = 'kegg', 
        organism = 'hsa',
        nPerm  = 1000,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust
        #by = "fgsea",
        #seed = T  
      )
      gse.DO <- gseDO(
        geneList = gene.list.ncbiid,
        nPerm  = 1000,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1,
        pAdjustMethod = pAdjust
        #by = "fgsea",
        #seed = T
      )
    }
    # 2.1 KEGG Gene Set Enrichment Analysis 
    message("start gsea at specific geneset")
    {
      m_t <- read.table("/BioII/lulab_b/zhanqing/cfRNAseq/reference/geneset/PATH_ID_NAME_KEGGplusHallmark.txt",sep = "\t",header = T, stringsAsFactors = F,row.names = 1)
      m_t2g <- m_t[,c("short_DESCRPTION","ensembl_gene_id")]
      #m_t2g <- m_t2g[!duplicated(m_t2g$short_DESCRPTION),]
      m_t2e <- m_t[,c("short_DESCRPTION","ENTREZID")]
      #m_t2e <- m_t2e[!duplicated(m_t2e$short_DESCRPTION),]
      gsea <- GSEA(gene.list.ncbiid,
                   TERM2GENE = m_t2e, 
                   nPerm  = 1000,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = pValue,
                   pAdjustMethod = pAdjust)
      gsea2 <- setReadable(gsea.kegg,
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID")
    }
    # output table
    {
      types <- c(gse.GO, gse.KEGG, gse.DO)
      names(types) <-c("gse.GO","gse.KEGG","gse.DO")
      
      outdir <- paste0("result/paired/enrich/nc-table-5114genes/")
      dir.create(paste0(outdir), showWarnings = T)
      message("start writing tables...")
      for(i in 1:length(types)){
        write.table(types[[i]],paste0(outdir,"/", names(types)[i], ".txt"),row.names = F,col.names = T, quote = F,sep = "\t")}
    }
    # GSEA plot of all samples
    {
      outdir="result/paired/enrich/nc-table-5114genes/"
      col = brewer.pal(8,"Set1")
      # 气泡图，展示geneset被激活还是抑制
      
      pdf(paste0(outdir,"dotplot-gse.GO",".pdf"), width = 8, height = 5.8)
      enrichplot::dotplot(gse.GO, split=".sign", font.size=10, title="gsea.GO", showCategory = 9) +
        facet_grid(~.sign) + 
        scale_color_continuous(low=col[1],high=col[2]) +
        scale_x_continuous(limits = c(0.18,0.85),breaks = c(0.2,0.4,0.6,0.8))
      dev.off()
      
      pdf(paste0(outdir,"dotplot-gse.KEGG",".pdf"), width = 6, height = 5.8)
      enrichplot::dotplot(gse.KEGG, split=".sign", font.size=10, title="gsea.KEGG", showCategory = 10) +
        facet_grid(~.sign) + 
        scale_color_continuous(low=col[1],high=col[2]) 
      scale_x_continuous(limits = c(0.25,0.55),breaks = c(0.3,0.4,0.5))
      dev.off()
      
      
      
      # GSEA plot
      col = brewer.pal(8,"Dark2")
      pdf(paste0(outdir,"gseaplot.GO.up",".pdf"), width = 8, height = 6)
      a = gse.GO@result[gse.GO@result$p.adjust<0.1,]
      up.id = a$ID[order(a$NES, decreasing = T)][1:6]
      #up.id = c("GO:0071352","GO:0006482","GO:0006968","GO:0031430")
      gseaplot2(gse.GO, geneSetID = up.id, color = col[1:6], pvalue_table = F, 
                base_size = 12, rel_heights = c(1.5, 0.4, 0.6), subplots = 1:3)
      dev.off()
      
      pdf(paste0(outdir,"gseaplot.GO.down",".pdf"), width = 8, height = 6)
      a = gse.GO@result[gse.GO@result$p.adjust<0.1,]
      down.id = a$ID[order(a$NES, decreasing = F)][1:6]
      gseaplot2(gse.GO, geneSetID = down.id, color = col[1:6], pvalue_table = F, 
                base_size = 12, rel_heights = c(1.5, 0.4, 0.6), subplots = 1:3)
      dev.off()
      
      pdf(paste0(outdir,"gseaplot.KEGG",".pdf"), width = 8, height = 6)
      gseaplot2(gse.KEGG, geneSetID = 1:5, color = col[1:5], pvalue_table = F, 
                base_size = 10, rel_heights = c(1.5, 0.5, 0.8), subplots = 1:3)
      dev.off()
      
      pdf(paste0(outdir,"gseaplot.DO",".pdf"), width = 8, height = 6)
      gseaplot2(gse.DO, geneSetID = 1:5, color = col[1:5], pvalue_table = T, 
                base_size = 10, rel_heights = c(1.5, 0.5, 0.8), subplots = 1:3)
      dev.off()
    }
    {
      nc.go <- read.table("result/paired/enrich/nc-table-5114genes/gse.GO.txt", header=T, sep="\t")
      nc.go <- nc.go %>% mutate(Count = str_count(core_enrichment,"/")+1)
      nc.go.up <- nc.go[nc.go$NES > 1,]
      nc.go.up.order <- nc.go.up[order(-nc.go.up$p.adjust, nc.go.up$NES, decreasing = T),]
      nc.go.up.order <- nc.go.up.order[1:20,]
      
      nc.go.down <- nc.go[nc.go$NES < -1,]
      nc.go.down.order <- nc.go.down[order(-nc.go.down$p.adjust, -nc.go.down$NES, decreasing = T),]
      nc.go.down.order <- nc.go.down.order[1:20,]
      
      go.up <- nc.go.down.order %>% mutate(logp = -log10(p.adjust))
      go.up$ONTOLOGY <- factor(go.up$ONTOLOGY, levels=c("CC","BP","MF"))
      go.up <- go.up[order(go.up$ONTOLOGY,go.up$p.adjust),]
      des <- go.up$Description[!duplicated(go.up$Description)]
      go.up$Description <- factor(go.up$Description, levels=rev(des))
      
      nc.kegg <- read.table("result/paired/enrich/nc-table-5114genes/gse.KEGG.txt", header=T, sep="\t")
      nc.kegg <- nc.kegg %>% mutate(Count = str_count(core_enrichment,"/")+1)
      nc.kegg.up <- nc.kegg[nc.kegg$NES > 1,]
      nc.kegg.up.order <- nc.kegg.up[order(-nc.kegg.up$p.adjust, nc.kegg.up$NES, decreasing = T),]
      nc.kegg.up.order <- nc.kegg.up.order[1:20,]
      
      nc.kegg.down <- nc.kegg[nc.kegg$NES < -1,]
      nc.kegg.down.order <- nc.kegg.down[order(-nc.kegg.down$p.adjust, -nc.kegg.down$NES, decreasing = T),]
      nc.kegg.down.order <- nc.kegg.down.order[1:10,]
      
      kegg.up <- nc.kegg.down.order %>% mutate(logp = -log10(p.adjust))
      kegg.up <- kegg.up[order(kegg.up$p.adjust),]
      des <- kegg.up$Description[!duplicated(kegg.up$Description)]
      kegg.up$Description <- factor(kegg.up$Description, levels=rev(des))
      kegg.up$ONTOLOGY = "kegg"
      
      
      {
        p=ggplot(kegg.up, aes(logp, Description)) +
          geom_point(aes(color = ONTOLOGY, size = Count)) +
          scale_color_aaas()+
          scale_x_continuous(limits=c(0.2,0.83))+
          #scale_color_gradient(low=col[1],high=col[2],guide="colorbar")+
          theme_bw()+
          guides(fill=guide_legend(title=NULL, ncol=1))+
          #coord_fixed(ratio=0.5)+
          theme(
            plot.margin = unit(c(1, 1, 1, 1),"cm"),
            panel.border=element_rect(size=1, fill="transparent"),
            panel.grid=element_line(size=rel(0.5)),
            #legend.position=c(0.25,0.85),
            legend.position="right",
            legend.background = element_blank(),
            legend.title = element_text(color="black", size=10),
            legend.text= element_text(color="black", size=10),
            axis.text.x = element_text(color="black", size=10),
            axis.text.y = element_text(color="black", size=8),
            axis.title.x = element_text(color="black", size=10),
            axis.title.y = element_blank())+
          labs(title="KEGG.down")+
          xlab("p.adjust")
        p
        ggsave(paste0(outdir,"dotplot-gse.kegg.down",".pdf"), width = 5.5, height = 5)
      }
    }
    # GSEA plot of each cancer samples
    {
      outdir="result/enrich/cancer-plot/"
      col = brewer.pal(8,"Set1")
      
      nc.go.gse <- read.table("result/enrich/nc-table/gse.GO.txt",sep="\t",header=T,quote="")
      nc.kegg.gse <- read.table("result/enrich/nc-table/gse.KEGG.txt",sep="\t",header=T,quote="")
      nc.do.gse <- read.table("result/enrich/nc-table/gse.DO.txt",sep="\t",header=T,quote="")
      crc.go.gse <- read.table("result/enrich/crc-table/gse.GO.txt",sep="\t",header=T,quote="")
      crc.kegg.gse <- read.table("result/enrich/crc-table/gse.KEGG.txt",sep="\t",header=T,quote="")
      crc.do.gse <- read.table("result/enrich/crc-table/gse.DO.txt",sep="\t",header=T,quote="")
      luad.go.gse <- read.table("result/enrich/luad-table/gse.GO.txt",sep="\t",header=T,quote="")
      luad.kegg.gse <- read.table("result/enrich/luad-table/gse.KEGG.txt",sep="\t",header=T,quote="")
      luad.do.gse <- read.table("result/enrich/luad-table/gse.DO.txt",sep="\t",header=T,quote="")
      
      {
        go.up.1 <- nc.go.gse[nc.go.gse$NES > 0 & order(nc.go.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "NC")
        go.up.2 <- crc.go.gse[crc.go.gse$NES > 0 & order(crc.go.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "CRC")
        go.up.3 <- luad.go.gse[luad.go.gse$NES > 0 & order(luad.go.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "LUAD")
        go.down.1 <- nc.go.gse[nc.go.gse$NES < 0 & order(nc.go.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "NC")
        go.down.2 <- crc.go.gse[crc.go.gse$NES < 0 & order(crc.go.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "CRC")
        go.down.3 <- luad.go.gse[luad.go.gse$NES < 0 & order(luad.go.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "LUAD")
        df <- rbind(go.up.1, go.up.2, go.up.3)
        df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
        df$Description <- factor(df$Description, levels=rev(df$Description))
        df <- df %>% mutate(Count = str_count(df$core_enrichment,"/")+1)
        df <- rbind(go.down.1, go.down.2, go.down.3)
        df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
        des <- df$Description[!duplicated(df$Description)]
        df$Description <- factor(df$Description, levels=rev(des))
        df <- df %>% mutate(Count = str_count(df$core_enrichment,"/")+1)
        
        kegg.up.1 <- nc.kegg.gse[nc.kegg.gse$NES > 0 & order(nc.kegg.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "NC")
        kegg.up.2 <- crc.kegg.gse[crc.kegg.gse$NES > 0 & order(crc.kegg.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "CRC")
        kegg.up.3 <- luad.kegg.gse[luad.kegg.gse$NES > 0 & order(luad.kegg.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "LUAD")
        kegg.down.1 <- nc.kegg.gse[nc.kegg.gse$NES < 0 & order(nc.kegg.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "NC")
        kegg.down.2 <- crc.kegg.gse[crc.kegg.gse$NES < 0 & order(crc.kegg.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "CRC")
        kegg.down.3 <- luad.kegg.gse[luad.kegg.gse$NES < 0 & order(luad.kegg.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "LUAD")
        df <- rbind(kegg.up.1, kegg.up.2, kegg.up.3)
        df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
        des <- df$Description[!duplicated(df$Description)]
        df$Description <- factor(df$Description, levels=rev(des))
        df <- df %>% mutate(Count = str_count(df$core_enrichment,"/")+1)
        df <- rbind(kegg.down.1, kegg.down.2, kegg.down.3)
        df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
        des <- df$Description[!duplicated(df$Description)]
        df$Description <- factor(df$Description, levels=rev(des))
        df <- df %>% mutate(Count = str_count(df$core_enrichment,"/")+1)
        
        do.up.1 <- nc.do.gse[nc.do.gse$NES > 0 & order(nc.do.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "NC")
        do.up.2 <- crc.do.gse[crc.do.gse$NES > 0 & order(crc.do.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "CRC")
        do.up.3 <- luad.do.gse[luad.do.gse$NES > 0 & order(luad.do.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "LUAD")
        do.down.1 <- nc.do.gse[nc.do.gse$NES < 0 & order(nc.do.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "NC")
        do.down.2 <- crc.do.gse[crc.do.gse$NES < 0 & order(crc.do.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "CRC")
        do.down.3 <- luad.do.gse[luad.do.gse$NES < 0 & order(luad.do.gse$p.adjust, decreasing = F),][1:10,] %>% mutate(group = "LUAD")
        df <- rbind(do.up.1, do.up.2, do.up.3)
        df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
        des <- df$Description[!duplicated(df$Description)]
        df$Description <- factor(df$Description, levels=rev(des))
        df <- df %>% mutate(Count = str_count(df$core_enrichment,"/")+1)
        df <- rbind(do.down.1, do.down.2, do.down.3)
        df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
        des <- df$Description[!duplicated(df$Description)]
        df$Description <- factor(df$Description, levels=rev(des))
        df <- df %>% mutate(Count = str_count(df$core_enrichment,"/")+1)
      }
      {
        id.go.gse.up <- Reduce(union,list(nc.go.gse[nc.go.gse$NES > 0,"ID"][1:10],
                                          crc.go.gse[crc.go.gse$NES > 0,"ID"][1:10],
                                          luad.go.gse[luad.go.gse$NES > 0,"ID"][1:10]))
        df <- rbind(nc.go.gse[nc.go.gse$ID %in% id.go.gse.up,] %>% mutate(group = "NC"),
                    crc.go.gse[crc.go.gse$ID %in% id.go.gse.up,] %>% mutate(group = "CRC"),
                    luad.go.gse[luad.go.gse$ID %in% id.go.gse.up,] %>% mutate(group = "LUAD"))
        df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
        df$Description <- factor(df$Description, levels= rev(df$Description[match(id.go.gse.up, df$ID)]))
        df <- df %>% mutate(Count = str_count(df$core_enrichment,"/")+1)
        df <- df[df$p.adjust < 0.4,]
        
        id.go.gse.down <- Reduce(union,list(nc.go.gse[nc.go.gse$NES < 0,"ID"][1:10],
                                            crc.go.gse[crc.go.gse$NES < 0,"ID"][1:10],
                                            luad.go.gse[luad.go.gse$NES < 0,"ID"][1:10]))
        df <- rbind(nc.go.gse[nc.go.gse$ID %in% id.go.gse.down,] %>% mutate(group = "NC"),
                    crc.go.gse[crc.go.gse$ID %in% id.go.gse.down,] %>% mutate(group = "CRC"),
                    luad.go.gse[luad.go.gse$ID %in% id.go.gse.down,] %>% mutate(group = "LUAD"))
        df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
        df$Description <- factor(df$Description, levels= rev(df$Description[match(id.go.gse.down, df$ID)]))
        df <- df %>% mutate(Count = str_count(df$core_enrichment,"/")+1)
        df <- df[df$p.adjust < 0.4,]
        
        id.kegg.gse.up <- Reduce(union,list(nc.kegg.gse[nc.kegg.gse$NES > 0,"ID"][1:10],
                                            crc.kegg.gse[crc.kegg.gse$NES > 0,"ID"][1:10],
                                            luad.kegg.gse[luad.kegg.gse$NES > 0,"ID"][1:10]))
        df <- rbind(nc.kegg.gse[nc.kegg.gse$ID %in% id.kegg.gse.up,] %>% mutate(group = "NC"),
                    crc.kegg.gse[crc.kegg.gse$ID %in% id.kegg.gse.up,] %>% mutate(group = "CRC"),
                    luad.kegg.gse[luad.kegg.gse$ID %in% id.kegg.gse.up,] %>% mutate(group = "LUAD"))
        df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
        df$Description <- factor(df$Description, levels= rev(df$Description[match(id.kegg.gse.up, df$ID)]))
        df <- df %>% mutate(Count = str_count(df$core_enrichment,"/")+1)
        df <- df[df$p.adjust < 0.4,]
        
        id.kegg.gse.down <- Reduce(union,list(nc.kegg.gse[nc.kegg.gse$NES < 0,"ID"][1:10],
                                              crc.kegg.gse[crc.kegg.gse$NES < 0,"ID"][1:10],
                                              luad.kegg.gse[luad.kegg.gse$NES < 0,"ID"][1:10]))
        df <- rbind(nc.kegg.gse[nc.kegg.gse$ID %in% id.kegg.gse.down,] %>% mutate(group = "NC"),
                    crc.kegg.gse[crc.kegg.gse$ID %in% id.kegg.gse.down,] %>% mutate(group = "CRC"),
                    luad.kegg.gse[luad.kegg.gse$ID %in% id.kegg.gse.down,] %>% mutate(group = "LUAD"))
        df$group <- factor(df$group, levels=c("NC","CRC","LUAD"))
        df$Description <- factor(df$Description, levels= rev(df$Description[match(id.kegg.gse.down, df$ID)]))
        df <- df %>% mutate(Count = str_count(df$core_enrichment,"/")+1)
        df <- df[df$p.adjust < 0.4,]
      }
      
      {
        p=ggplot(df, aes(group, Description), showCatekeggry = 20) +
          geom_point(aes(color = p.adjust,size=Count)) +
          scale_color_gradient(low=col[1],high=col[2])+
          theme_bw()+
          guides(fill=guide_legend(title=NULL, ncol=1))+
          coord_fixed(ratio=0.7)+
          theme(
            plot.margin = unit(c(1, 1, 1, 1),"cm"),
            panel.border=element_rect(size=1, fill="transparent"),
            panel.grid=element_line(size=rel(0.5)),
            #legend.position=c(0.25,0.85),
            legend.position="right",
            legend.background = element_blank(),
            legend.title = element_text(color="black", size=12),
            legend.text= element_text(color="black", size=12),
            axis.text.x = element_text(color="black", size=10),
            axis.text.y = element_text(color="black", size=10),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())+
          labs(title="GSEA.KEGG.down")
        p
      }
      ggsave("result/enrich/cancer-plot/GSEA.KEGG.down-pvalue<0.4.pdf",p,width=7,height=5)
      
      
    }
    # GSEA plot of each cancer samples.2
    {
      # GSEA plot
      outdir="result/enrich/cancer-plot/"
      
      
      col=brewer.pal(8,"Set1")
      col=col[-6]
      go.id = c("GO:0004386","GO:0038096","GO:0000932",
                "GO:0000978","GO:0001503","GO:0005667")
      pdf(paste0(outdir,"gseaplot.GO.luad",".pdf"), width = 8, height = 7)
      gseaplot2(gse.GO.nc, geneSetID = go.id, color = col[1:6], pvalue_table = F, 
                base_size = 14, rel_heights = c(1.5, 0.4, 0.6), subplots = 1:3)
      gseaplot2(gse.GO.crc, geneSetID = go.id, color = col[1:6], pvalue_table = F, 
                base_size = 14, rel_heights = c(1.5, 0.4, 0.6), subplots = 1:3)
      gseaplot2(gse.GO.luad, geneSetID = go.id, color = col[1:6], pvalue_table = F, 
                base_size = 14, rel_heights = c(1.5, 0.4, 0.6), subplots = 1:3)
      dev.off()
      
      kegg.id = c("hsa05168","hsa04510","hsa04666",
                  "hsa05226","hsa04310","hsa03010")
      pdf(paste0(outdir,"gseaplot.kegg.luad",".pdf"), width = 8, height = 7)
      gseaplot2(gse.KEGG.nc, geneSetID = kegg.id, color = col[1:6], pvalue_table = F, 
                base_size = 14, rel_heights = c(1.5, 0.4, 0.6), subplots = 1:3)
      gseaplot2(gse.KEGG.crc, geneSetID = kegg.id, color = col[1:6], pvalue_table = F, 
                base_size = 14, rel_heights = c(1.5, 0.4, 0.6), subplots = 1:3)
      gseaplot2(gse.KEGG.luad, geneSetID = kegg.id, color = col[1:6], pvalue_table = F, 
                base_size = 14, rel_heights = c(1.5, 0.4, 0.6), subplots = 1:3)
      dev.off()
    }
    # 2.2 Cancer gene set Enrichment Analysis
    message("start gsea at cancer gene set")
    {
      CGC_1 <- read.table("id/CGC_1.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      CGC_2 <- read.table("id/CGC_2.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      CGC_all <- read.table("id/CGC_all.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      NCG_1 <- read.table("id/NCG_1.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      NCG_2 <- read.table("id/NCG_2.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      Uniprot <- read.table("id/Uniprot.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      Pansoft <- read.table("id/PanSoftware.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      IntOGene <- read.table("id/IntOGen.txt", sep = "\t", header = T, stringsAsFactors = F, check.names = F, quote = "")
      
      
      set1 <- read.table("v2/public/ensemblid/set/Set1.txt")[,1]
      set2 <- read.table("v2/public/ensemblid/set/Set2.txt")[,1]
      set3 <- read.table("v2/public/ensemblid/set/Set3.txt")[,1]
      set4 <- read.table("v2/public/ensemblid/set/Set4.txt")[,1]
      set5 <- read.table("v2/public/ensemblid/set/Set5.txt")[,1]
      set6 <- read.table("v2/public/ensemblid/set/Set6.txt")[,1]
      set7 <- read.table("v2/public/ensemblid/set/Set7.txt")[,1]
      
      geneset <- rbind(data.frame(Description = "CGC_1", Entrez = CGC_1$Entrez),
                       #data.frame(Description = "CGC_all", Entrez = CGC_all$Entrez),
                       data.frame(Description = "NCG_1", Entrez = NCG_1$ENTREZID),
                       #data.frame(Description = "NCG_all", Entrez = NCG_2$ENTREZID),
                       data.frame(Description = "Uniprot", Entrez = Uniprot$ENTREZID),
                       data.frame(Description = "Pansoft", Entrez = Pansoft$ENTREZID),
                       data.frame(Description = "IntOGene", Entrez = IntOGene$ENTREZID)
      )
      geneset2 <- rbind(data.frame(Description = "set1", ensg = set1),
                        data.frame(Description = "set2", ensg = set2),
                        data.frame(Description = "set3", ensg = set3),
                        data.frame(Description = "set4", ensg = set4),
                        data.frame(Description = "set5", ensg = set5),
                        data.frame(Description = "set6", ensg = set6),
                        data.frame(Description = "set7", ensg = set7))
      
      pValue = 1
      pAdjust = "BH"
      gsea <- clusterProfiler::GSEA(gene.list.ncbiid,
                                    TERM2GENE = geneset,
                                    nPerm  = 1000,
                                    minGSSize = 10,
                                    maxGSSize = 500,
                                    pvalueCutoff = pValue,
                                    pAdjustMethod = pAdjust)
      gsea2 <- setReadable(gsea,
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENRTEZID")  # map geneID to geneSymbol
      
      outdir="result/enrich/plot/"
      col = brewer.pal(8,"Set1")
      pdf(paste0(outdir,"gseaplot.CancerGeneSet",".pdf"), width = 8, height = 6)
      gseaplot2(gsea, geneSetID = 1:5, color = col[1:5], pvalue_table = T, 
                base_size = 10, rel_heights = c(1.5, 0.4, 0.7), subplots = 1:3)
      dev.off()
      
    }
    # fgsea
    {
      rank = gene.list.ensgid
      geneset = list(CGC_1 = CGC_1$Ensg,
                     CGC_all = CGC_2$Ensg,
                     NCG_1 = NCG_1$ENSEMBL,
                     NCG_2 = NCG_2$ENSEMBL,
                     Uniprot = Uniprot$ENSEMBL,
                     Pansoft = Pansoft$ENSEMBL,
                     IntOGene = IntOGene$ENSEMBL)
      library(fgsea)
      fgsea <- fgsea(geneset, 
                     rank, 
                     nperm = 10,
                     minSize = 0, 
                     maxSize = 500)
      summary(fgsea)
      
      plotGseaTable(fgsea$pathway,
                    rank,
                    fgsea, 
                    gseaParam=0.5)
      plotEnrichment(fgsea$pathway,rank,)
      
    }
    {
      # emapplot
      enrichplot::emapplot(gse.GO)
      # cnetplot
      cnetplot(gsea)
      # 山峦图，展示每个geneset的基因logFC分布
      ridgeplot(gsea)
      # 选择单个gene set
      gsea3 <- data.frame(gsea)
      gseaplot2(gsea, geneSetID = 1, title=gsea3$Description[1])
      
    }
    
  }
}
