library(ggplot2)
library(reshape)
library(pROC)
library(ggpubr)
library(dplyr)

#source("./model/loaddata.r")
od="./figures/fig4"; dir.create(od)
od2=paste0(od,"/roc"); dir.create(od2,showWarnings=F)

dirt="./cox"
tars=list( aml1=list("ddAML",14),
           mds1=list("ddMDS",7),
           mpn1=list("ddMPN",6))
ctimes=list( X1=c(0,365),X2=c(365,1825),X3=c(1825,Inf))


###########################
########################### Fig4.a-c
###########################
pd=data.frame()
for(ty in names(tars)){
  tarcn=tars[[ty]][[1]] 
  dis=ifelse(ty=="aml1","AML",ifelse(ty=="mds1","MDS","MPN") )
  t=read.table(paste0(dirt,"/",ty,"/genopheno_pred.txt"),row.names=1,head=1); ids=row.names(t)
  cids=row.names(datfgeno)%in%ids
  dftemp=data.frame(tdd[ids,], ref=ifelse(tdd[ids,tarcn]<tdd[ids,"dcen"],1,0), 
                    n=nmutrow[cids], aml1=caml1[cids], mds1=cmds1[cids], mpn1=cmpn1[cids],cmm=cmm[cids])
  for(tn in names(ctimes)){
    low=ctimes[[tn]][1]; high=ctimes[[tn]][2]; ref2=ifelse(dftemp[,tarcn]>low&dftemp[,tarcn]<high,1,0)
    if(tn=="X1"){ pred=t$X365; }else if(tn=="X2"){pred=(t$X365+t$X1825)/2}else{pred=(t$X1825+t$max)/2}
    tntag=ifelse(tn=="X1","0-1y",ifelse(tn=="X2","1-5y",">5y"))
    disurv=dftemp[,tarcn]
    dfsel = data.frame(dftemp,pred,ref2)
    rr=suppressMessages(roc(dfsel$ref2,-dfsel$pred))
    png(paste0(od2,"/",ty,"_",tn,".png"),height=350,width=350); 
    plot.roc(rr, main=paste0(ty,"_",tntag,"_",leng(dftemp[,ty]&disurv>low&disurv<high),"\nAUC=",round(auc(rr),digits=3) ));
    dev.off()
    pd=rbind(pd,data.frame(Dis=dis, Interval=tntag, Sensitivity=rr$sensitivities, Specificity=rr$specificities))
  }
}
pd$Dis=factor(pd$Dis,levels=c("AML","MDS","MPN"))
pd$Interval=factor(pd$Interval,levels=c(">5y","1-5y","0-1y"))

for(dis in unique(pd$Dis)){
  pdt=subset(pd,Dis==dis)
  ggplot(pdt, aes(y=Sensitivity,x=Specificity, color=Interval))+
    geom_line(aes(linetype=Interval,size=Interval)) + 
    scale_color_manual(values=c(col3,col1,col2)) + 
    scale_linetype_manual(values=c("solid","solid","solid")) +
    scale_size_manual(values=c(1,0.8,0.6))+
    scale_x_reverse() + 
    geom_abline(slope=1,intercept=1, col="grey",size=0.6, linetype="dashed")+
    ggtitle(dis) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       axis.text=element_text(size=12), axis.title=element_text(size=14),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  ggsave(paste0(od,"/fig4_a-c_roc_",dis,".pdf"),width=4.7,height=3.6)
}


#####################
##################### Fig4.d-f
#####################
tars=list( 
  AML=list(c("HGB","PLT","MCV","RDW","CYS"),caml1,tdd$ddAML),
  MDS=list(c("HGB","PLT","MCV","RDW","CYS"),cmds1,tdd$ddMDS),
  MPN=list(c("HGB","PLT","MCV","RDW","CYS"),cmpn1,tdd$ddMPN))
ctimes=list( t1=c(Inf,1825),t2=c(1825,365),t3=c(365,0))

datf=datfpheno
for(ty in names(tars)){
  tarcns=tars[[ty]][[1]]; cbg=tars[[ty]][[2]]; dd=tars[[ty]][[3]]
  pls=c()
  for(feat in tarcns){
    vs=tpheno[,feat]; vs[is.na(vs)]=phemed[feat]; mat=data.frame()
    mat=rbind( data.frame( Time="Control", Z=vs[ (1:nrow(datf))%in%c(1:500) & !cmm ] )); ns=c()
    for(tn in names(ctimes)){
      tntag=ifelse(tn=="t1",">5y",ifelse(tn=="t2","1-5y",ifelse(tn=="t3","0-1y","xx")))
      high=ctimes[[tn]][1]; low=ctimes[[tn]][2]
      vst=vs[cbg&dd>low&dd<high]; ns=c(ns,length(vst))
      mat=rbind(mat,data.frame(Time=tntag,Z=vst))
    }
    mat$Time <- factor(mat$Time , levels=c("Control",">5y","1-5y","0-1y"))
    xx = mat %>% group_by(Time) %>% summarise(y = quantile(mat$Z,c(0, 1))[1], Time = first(Time))
    print(paste0(ty,"  ",feat,"  ",wilcox.test(mat$Z[mat$Time=="Control"],mat$Z[mat$Time==">5y"])[[3]]))
    print(paste0(ty,"  ",feat,"  ",wilcox.test(mat$Z[mat$Time=="Control"],mat$Z[mat$Time=="1-5y"])[[3]]))
    print(paste0(ty,"  ",feat,"  ",wilcox.test(mat$Z[mat$Time=="Control"],mat$Z[mat$Time=="0-1y"])[[3]]))
    pl=
      ggplot(mat, aes(x=Time, y=Z, fill=Time)) + 
      geom_boxplot(outlier.shape=NA, fatten = NULL) + ylab(feat) + 
      stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                   width = 0.75, size = 1, linetype = "solid")+
      coord_cartesian(ylim = quantile(mat$Z, c(0.005, 0.99)))+
      scale_fill_manual(values=c(col4,col3,col1,col2)) + xlab("")+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         axis.text=element_text(size=12),
                         axis.text.x=element_text(size=12,angle=45,hjust=1),legend.position = "none", 
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    pls=c(pls,list(pl))
  }
  ggarrange(pls[[1]],pls[[2]],pls[[3]],pls[[4]],pls[[5]] , ncol = 5, nrow = 1)
  ggsave(paste0(od,"/fig4_d-f_boxplot_",ty,".pdf"),width=7.5, height=4)
}

