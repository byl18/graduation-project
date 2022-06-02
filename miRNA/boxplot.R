library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsignif)

FTA_box <- c('FTA-1','FTA-2','FTA-3','FTA-4','FTA-6','FTA-9','FTA-10','FTA-11','FTA-12','FTA-13','FTA-15','FTA-16','FTA-17','FTA-18','FTA-20','FTA-23','FTA-24')
FTC_box <- c('FTC-1','FTC-2','FTC-3','FTC-4','FTC-5','FTC-6','FTC-7','FTC-8','FTC-10','FTC-11','FTC-13','FTC-14','FTC-15','FTC-16','FTC-17','FTC-18','FTC-19','FTC-20','FTC-21','FTC-22')
csFTA_box <- c('csFTA-1','csFTA-2','csFTA-3','csFTA-4','csFTA-6','csFTA-9','csFTA-10','csFTA-11','csFTA-12','csFTA-13','csFTA-15','csFTA-16','csFTA-17','csFTA-18','csFTA-20','csFTA-23','csFTA-24')
csFTC_box <- c('csFTC-1','csFTC-2','csFTC-3','csFTC-4','csFTC-5','csFTC-6','csFTC-7','csFTC-8','csFTC-10','csFTC-11','csFTC-13','csFTC-14','csFTC-15','csFTC-16','csFTC-17','csFTC-18','csFTC-19','csFTC-20','csFTC-21','csFTC-22')

exp = matrix(rnorm(60),nrow = 10)
colnames(exp) <- paste0("sample",1:6)
rownames(exp) <- paste0("gene",1:10)
exp[1:4,1:4]
dat = data.frame(t(exp))
dat = mutate(dat,group = rep(c("A","B"),each = 3))%>% mutate(pair = rep(c("AA","BB","cc"),each = 2))

dat2 = gather(dat,key = "gene",value = "expression",-group,-pair)
ggplot(data = dat2)+
  geom_boxplot(aes(x = group,y = expression,color = group))+
  theme_bw()+
  facet_wrap(~gene,nrow = 2)

remove(exp,)

compaired <- list(c("FTC", "FTA"),c("csFTC","csFTA"))
for(i in 1:length(pf_EV_up)){
  mydat <- (log(t(rawdata[pf_EV_up,])+1,2))
  mydat1 <- cbind(mydat[c(FTC_box,FTA_box),i],mydat[c(csFTC_box,csFTA_box),i])
  colnames(mydat1) <- c('EV',paste('Cell','-','free',sep=''))
  mydat1 <- data.frame(mydat1)
  mydat1 <- mutate(mydat1,group = c(rep("FTC",20),rep("FTA",17)))
  mydatg <- gather(mydat1,key = "EV_cf",value = "log2RPM",-group)
  mydatg$label <- c(rep(0.017,37),rep(0.38,37))
  
  mynew <- mydatg[,c(3,2)]
  mynew[,2] <-  c(rep("FTC",20),rep("FTA",17),rep("csFTC",20),rep("csFTA",17))
  colnames(mynew) <- c('log2RPM','group')
  mydatg$EV_cf[38:74] <- rep('Cell-free',37)
  # ggplot(mynew,aes(group,expression,fill=group))+
  #   geom_boxplot(width=0.5)+
  #   theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+
  #   labs(y= 'Expression')+
  #   geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)+
  #   ggtitle(colnames(mydat)[i]) +
  #   theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
  
  ggplot(data = mydatg)+
    geom_boxplot(aes(x = group,y = log2RPM,color = group))+
    geom_text(x='FTA', y=5, label=paste('FDR = ',mydatg$label,sep=''))+
   # geom_signif(comparisons = c('FTA','FTC'),step_increase = 0.1,map_signif_level = F,test = t.test)+
    theme_bw()+
    facet_wrap(~EV_cf,ncol = 2)+
    ggtitle(colnames(mydat)[i]) +
    theme(plot.title = element_text(hjust = 0.5))
    #+geom_signif(comparisons = c('FTA','FTC'),step_increase = 0.1,map_signif_level = F,test = t.test)+
}

uptxt <- read.table("data/up.txt",sep="\t",header=F,check.names=F)
rownames(uptxt) <- uptxt[,1]
uptxt <- uptxt[,-1]
uptxt[colnames(mydat)[i],]




mydat <- (log(t(rawdata[pf_EV_up,])+1,2))
mydat1 <- cbind(mydat[c(FTC_box,FTA_box),i],mydat[c(csFTC_box,csFTA_box),i])
colnames(mydat1) <- c('EV',paste('Cell','-','free',sep=''))
mydat1 <- data.frame(mydat1)
mydat1 <- mutate(mydat1,group = c(rep("FTC",20),rep("FTA",17)))
mydatg <- gather(mydat1,key = "EV_cf",value = "log2RPM",-group)
mydatg$label <- c(rep(0.072,37),rep(0.93,37))

mynew <- mydatg[,c(3,2)]
mynew[,2] <-  c(rep("FTC",20),rep("FTA",17),rep("csFTC",20),rep("csFTA",17))
colnames(mynew) <- c('log2RPM','group')
mydatg$EV_cf[38:74] <- rep('Cell-free',37)
# ggplot(mynew,aes(group,expression,fill=group))+
#   geom_boxplot(width=0.5)+
#   theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+
#   labs(y= 'Expression')+
#   geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)+
#   ggtitle(colnames(mydat)[i]) +
#   theme(plot.title = element_text(hjust = 0.5)) #设置标题居中

ggplot(data = mydatg)+
  geom_boxplot(aes(x = group,y = log2RPM,color = group))+
  geom_text(x='FTA', y=6, label=paste('FDR = ',mydatg$label,sep=''))+
  # geom_signif(comparisons = c('FTA','FTC'),step_increase = 0.1,map_signif_level = F,test = t.test)+
  theme_bw()+
  facet_wrap(~EV_cf,ncol = 2)+
  ggtitle(colnames(mydat)[i]) +
  theme(plot.title = element_text(hjust = 0.5))



