###TBS vs cancer####
library(psych)
setwd("C:/study/lulab_research/毕设/zq/miRNA")
TBS=read.csv("data/TBS.txt",sep='\t',header = TRUE)
colnames(TBS) <- gsub("\\.","-",colnames(TBS))

EV_TBS <- logRPM[,1:21]
cf_TBS <- logRPM[,22:42]
colnames(cf_TBS) <- gsub("cs","",colnames(cf_TBS))

EV_TBS <- EV_TBS[,colnames(TBS)]
cf_TBS <- cf_TBS[,colnames(TBS)]
cor_TBS <- EV_TBS[,1:2]
colnames(cor_TBS) <- c('EV','cf')

for(i in 1:dim(EV_TBS)[1]){
  cor_TBS[i,1] <- corr.test(data.frame(t( rbind(TBS[1,],EV_TBS[i,]) )),method="spearman",adjust="none")$r[1,2]
  cor_TBS[i,2] <- corr.test(data.frame(t( rbind(TBS[1,],cf_TBS[i,]) )),method="spearman",adjust="none")$r[1,2]
}

cor_TBS <- na.omit(cor_TBS)
cor_TBS <- data.frame(cor_TBS)
# cor_TBS[order(cor_TBS$EV,decreasing = TRUE),]
# cor_TBS[order(cor_TBS$cf,decreasing = TRUE),]
# common_EV_cf <- intersect(rownames(res.ev.sig),rownames(res.cf.sig))

cor_EV <- cor_TBS[order(cor_TBS$EV,decreasing = TRUE),][rownames(res.ev.sig),]
cor_cf <- cor_TBS[order(cor_TBS$cf,decreasing = TRUE),][rownames(res.cf.sig),]
cor_EV <- cor_EV[order(cor_EV$EV,decreasing = TRUE),]
cor_cf <- cor_cf[order(cor_cf$cf,decreasing = TRUE),]


# cor_TBS[order(cor_TBS$EV,decreasing = TRUE),][common_EV_cf,1]
# cor_TBS[order(cor_TBS$cf,decreasing = TRUE),][common_EV_cf,2]

for(i in 1:12){
  tmp <- data.frame(EV_TBS[rownames(cor_EV)[i],1:21],t(TBS[1,]))
  colnames(tmp) <- c('RNA','TBS')
  attach(tmp)
  p <- ggplot(tmp,aes(x = RNA,y=TBS))+geom_point(size = 4,shape=19)+
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
    labs(x="RNA reads",y="TBS",title= rownames(cor_EV)[i] , face="bold")+
    geom_text(x = 12,y=300,label=paste('spearman r = ',cor_EV[i,1]),size=3)
  p
  ggsave(p,filename = paste("./output/figure/TBS/EVpos/", rownames(cor_EV)[i],".png",sep=''),width = 4,height = 4)
}

for(i in (dim(cor_EV)[1]-5):dim(cor_EV)[1]){
  tmp <- data.frame(t(EV_TBS[rownames(cor_EV)[i],1:21]),t(TBS[1,]))
  colnames(tmp) <- c('RNA','TBS')
  attach(tmp)
  p <- ggplot(tmp,aes(x = RNA,y=TBS))+geom_point(size = 4,shape=19)+
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
    labs(x="RNA reads",y="TBS",title= rownames(cor_EV)[i] , face="bold")+
    geom_text(x = max(EV_TBS[rownames(cor_EV)[i],1:21]) * 0.8,y=300,label=paste('spearman r = ',cor_EV[i,1]),size=3)
  p
  ggsave(p,filename = paste("./output/figure/TBS/EVneg/", rownames(cor_EV)[i],".png",sep=''),width = 4,height = 4)
}

for(i in 1:3){
  tmp <- data.frame(cf_TBS[rownames(cor_cf)[i],1:21],t(TBS[1,]))
  colnames(tmp) <- c('RNA','TBS')
  attach(tmp)
  p <- ggplot(tmp,aes(x = RNA,y=TBS))+geom_point(size = 4,shape=19)+
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
    labs(x="RNA reads",y="TBS",title= rownames(cor_cf)[i] , face="bold")+
    geom_text(x = max(cf_TBS[rownames(cor_cf)[i],1:21]) * 0.8,y=300,label=paste('spearman r = ',cor_cf[i,2]),size=3)
  p
  ggsave(p,filename = paste("./output/figure/TBS/cfpos/", rownames(cor_cf)[i],".png",sep=''),width = 4,height = 4)
}

for(i in (dim(cor_cf)[1]-10) : dim(cor_cf)[1]){
  gene. <- unlist(lapply(strsplit(rownames(cor_cf)[i],"|",fixed=T),function(x) x[3]))
  tmp <- data.frame(t(cf_TBS[rownames(cor_cf)[i],1:21]),t(TBS[1,]))
  colnames(tmp) <- c('RNA','TBS')
  attach(tmp)
  p <- ggplot(tmp,aes(x = RNA,y=TBS))+geom_point(size = 4,shape=19)+
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
    labs(x="RNA reads",y="TBS",title= gene. , face="bold")+
    geom_text(x = max(cf_TBS[rownames(cor_cf)[i],1:21]) * 0.8,y=300,label=paste('spearman r = ',cor_cf[i,2]),size=3)
  p
  ggsave(p,filename = paste("./output/figure/TBS/cfneg/", gene.,".png",sep=''),width = 4,height = 4)
}

