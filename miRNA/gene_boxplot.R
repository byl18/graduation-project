######################################################
## differential expression gene boxplot
## 20220218
library(ggplot2)
library(ggpubr)
# cancer <- c(rep('FTA',20),rep('FTC',21),rep('FTC',20),rep('FTA',17))
# class <- c(rep('EV',41),rep('cf',37))
cancer <- c(rep('FTA',20),rep('FTC',21),rep('FTA',18),rep('FTC',21))
class <- c(rep('EV',41),rep('cf',39))

geneset = read.table("C:/Users/86185/bioinf/lulab/ML_study/model/result/miRNA/EV.txt", sep = "\t", header = F, stringsAsFactors = F, check.names = F)
geneset = read.table("C:/Users/86185/bioinf/lulab/ML_study/model/result/miRNA/EV_name.txt", sep = "\t", header = F, stringsAsFactors = F, check.names = F)

gene. = c()
for(gene in geneset$V1[1:10]){
  if(is.na(unlist(lapply(strsplit(gene,"|",fixed=T),function(x) x[4])))){
    gene. = append(gene.,gene)
    
  }
  else{
    gene. = append(gene.,unlist(lapply(strsplit(gene,"|",fixed=T),function(x) x[4])))

    }

}
print(gene.)


for(gene in geneset$V1[1:10]){
gene= geneset$V1[1:10][10]
gene. = unlist(lapply(strsplit(gene,"::",fixed=T),function(x) x[1]))

df = data.frame(tpm = (logRPM[gene,]), 
                #cpm = RPM[gene,],
                cancer = class,
                class = cancer) #patient = patient
colnames(df)[1] <- 'tpm'
# unpaired wilcox test
{
  p=ggplot(df, aes(x=cancer,y=tpm,fill=class))+  # cpm or tpm
    geom_boxplot(notch = FALSE, alpha = 0.9, size=0.4, outlier.shape = NA, position=position_dodge(0.9)) +
    #geom_jitter(position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.9), alpha = 0.8, size=0.8)+
    #scale_color_brewer(palette="Set1")+
    scale_fill_brewer(palette="Set1")+
    #scale_color_aaas()+
    ##########修改图例内容
    theme_bw()+
    theme(#plot.margin = unit(c(1, 1, 1, 1),"cm"),
      panel.grid = element_blank(),
      panel.border=element_rect(size=1,fill="transparent"),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(face="bold", color="black", size=12),
      plot.title = element_text(hjust = 0.5, size=12, face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=12),
      axis.text.y = element_text(face="bold",  color="black", size=12),
      axis.title.x = element_text(face="bold", color="black", size=12),
      axis.title.y = element_text(face="bold",color="black", size=12))+
    labs(x="",y="log(TPM)",title=gene., face="bold")+
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)),method = "t.test")+
    #stat_compare_means(aes(label =  ..p.signif..),label.x =seq(from=0.8,to=1.8,by=1), label.y = 17,method = "wilcox.test")+
    geom_signif(annotations = c("ns","**"),
                y_position =rep(11.5,2),  # change
                xmin = seq(from=0.8,to=1.8,by=1),
                xmax = seq(from=1.2,to=2.2,by=1),
                tip_length = rep(c(0.02,0.02),2),
                textsize=4)
  p
  ggsave(p,filename = paste("./output/figure/top10_miRNA/",gene.,".png",sep=''),width = 4,height = 4)
  # p>0.05 : ns
  # p<0.05 : * 
  # p<0.01 : **
  # p<0.001 : ***
  # p<0.0001 : ****
}
}
# paired test
{
  p=ggplot(df, aes(x=class,y=tpm,fill=class))+  # cpm or tpm
    facet_wrap(~ cancer)+
    geom_boxplot(notch = FALSE, alpha = 0.9, size=0.4, outlier.shape = NA, position=position_dodge(0.9)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), alpha = 0.8, size=0.8)+
    geom_line(aes(group=patient), size=0.4, color="darkgrey", alpha=0.8,linetype="dashed")+
    #scale_color_brewer(palette="Set1")+
    scale_fill_brewer(palette="Set1")+
    #scale_color_aaas()+
    #scale_y_continuous(limits = c(0,850))+
    ##########修改图例内容
    theme_bw()+
    theme(
      plot.margin = unit(c(0.5, 1, 0.5, 0.5),"cm"),
      panel.grid = element_blank(),
      panel.border=element_rect(size=1,fill="transparent"),
      legend.position = "none",
      legend.title = element_blank(),
      legend.text = element_text(face="bold", color="black", size=12),
      plot.title = element_text(hjust = 0.5, size=12, face="bold"),
      axis.text.x = element_text(face="bold", color="black", size=12),
      axis.text.y = element_text(face="bold",  color="black", size=12),
      axis.title.x = element_text(face="bold", color="black", size=12),
      axis.title.y = element_text(face="bold",color="black", size=12))+
    labs(x="",y="TPM",title=gene, face="bold")  
  p
  p <- p+stat_compare_means(aes(label = paste0("p = ", ..p.format.., "\n    ", ..p.signif..)),
                            paired = TRUE,  ########## paired wilcox test
                            method = "wilcox.test",
                            label.x = 1.1,
                            label.y = 500)+   # change coordinate everytime
    geom_signif(annotations = c(""),
                y_position = 560,   # change coordinate everytime
                xmin = 1, 
                xmax = 2,
                tip_length = c(0.02,0.02),
                textsize=4)
  p
}
#ggsave(paste("L1-L21/diff_exp/unnormalized/paired/diff_gene_",gene,".pdf",sep=""),p,width = 4.3,height = 4)

