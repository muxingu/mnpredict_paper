library(dplyr)
library(ggplot2)
#source("./model/loaddata.r")
od="./figures/fig1"; dir.create(od)
targn=c("TET2","ASXL1","TP53","JAK2","SRSF2","GNB1","SF3B1","GNAS","CALR","IDH2")

CI_t <- function (x, ci = 0.95){
  `%>%` <- magrittr::`%>%`
  Margin_Error <- qt(ci + (1 - ci)/2, df = length(x) - 1) * sd(x)/sqrt(length(x))
  df_out <- data.frame( sample_size=length(x), Mean=mean(x), sd=sd(x),
                        Margin_Error=Margin_Error,
                        'CI lower limit'=(mean(x) - Margin_Error),
                        'CI Upper limit'=(mean(x) + Margin_Error)) %>%
    tidyr::pivot_longer(names_to = "Measurements", values_to ="values", 1:6 )
  return(df_out)
}


##############
############## Fig.1a
##############

pdat=data.frame(DNMT3A=apply(datfgenowp[,c("DNMT3A_R882","DNMT3A_other")],1,max), datfgenowp[,targn])
npdat=colSums(pdat)
lnp=log10(npdat)
pnp = paste0(round(npdat/leng(nmutrow>0)*100,digits=1),"%")

bardf=data.frame(x=names(lnp),lnp, pnp); bardf$x=factor(bardf$x,levels=rev(names(lnp)))
ggplot(bardf, aes(x=x,y=lnp))+
  geom_bar( stat="identity", fill=col1, width=0.8 ) +
  geom_text(aes(label=pnp), hjust=-1, col="white")+
  coord_flip() + scale_y_reverse() + ylab("Log10(# of individuals)")+xlab("")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
           axis.text.x=element_text(size=12),legend.position = "none", 
           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste0(od,"/fig_1a.pdf"), height=4,width=3.5)
  

########################
######################## Fig.1b
########################
pdat=data.frame(DNMT3A=apply(datfgenovaf[,c("DNMT3A_R882","DNMT3A_other")],1,max), datfgenovaf[,targn])
pdvaf=data.frame()
for(g in colnames(pdat)){
  vs=pdat[,g]
  pdvaf=rbind(pdvaf, data.frame(gene=g,vaf=vs[vs>0]))
}
pdvaf$gene=factor(pdvaf$gene,levels=rev(names(lnp)))
ggplot(pdvaf, aes(x=gene, y=vaf)) + 
  geom_violin(trim=FALSE, fill=alpha(col2,0.5),col=alpha(col2,0.5), draw_quantiles = c(0.25,0.75)) +
  #geom_violin(draw_quantiles = c(0.25,0.75), fill="transparent")+
  stat_summary(fun=median, geom ="point")+
  coord_flip() + xlab("")+ylab("VAF")+ ylim(0,0.5)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
         axis.text.x=element_text(size=12),legend.position = "none", 
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste0(od,"/fig_1b.pdf"), height=4,width=5)


##############################
############################## Fig.1c
##############################
pdat=data.frame()
for(n in c(40:70)){
  cx=tpheno$age==n
  perc=leng(cx&nmutrow>0)/leng(cx)*100
  pdat=rbind(pdat,data.frame(age=n,perc))
}
head(pdat)
x=pdat$age; y=pdat$perc
pdat2=data.frame(age=x, perc=predict(loess(y ~ x)))

ggplot(pdat2,aes(x=age,y=perc)) + 
  geom_line(col=col1)+ xlab("Age")+ylab("CH Prevalence (%)")+ ylim(0,10)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     axis.text=element_text(size=12),legend.position = "none", 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste0(od,"/Fig_1c.pdf"), height=4,width=5)


################################
################################ Fig.1d
################################

pdat=data.frame()
for(n in seq(40,70,1)){
  cx=tpheno$age==n
  cl=vaf[cx]; clnz=cl[cl>0]
  ci = CI_t(clnz)
  qs=quantile(clnz,probs=seq(0,1,0.05))
  med=mean(clnz); low=ci$values[5]; high=ci$values[6]
  pdat=rbind(pdat,data.frame(age=n,med,low,high))
}
pdat[3,2]=0.1;
head(pdat,n=10)
x=pdat$age; y=pdat$med; low=pdat$low; high=pdat$high
pdat2=data.frame(age=x, med=predict(loess(y ~ x)), low=predict(loess(low ~ x)), high=predict(loess(high ~ x)))

ggplot(pdat2,aes(x=age,y=med)) + 
  geom_ribbon(aes(ymin=low,ymax=high), fill="gray70")+
  geom_line(col=col1)+ xlab("Age")+ylab("VAF") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     axis.text=element_text(size=12),legend.position = "none", 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste0(od,"/fig_1d.pdf"), height=4,width=5)



#########################
######################### Fig.1e
#########################
pdat=data.frame()
for(i in 1:4){
  pdat=rbind(pdat,data.frame(n=i,nsamp=leng(nmutrow==i)))
}
pdat=rbind(pdat,data.frame(n=">4",nsamp=leng(nmutrow>4)))
pdat$n=factor(pdat$n,c("1","2","3","4",">5"))

ggplot(pdat, aes(x=n,y=nsamp))+
  geom_bar( stat="identity", fill=col1, width=0.8 ) +
  geom_text(aes(label=nsamp),vjust=-1, size=3) + ylim(0,22000)+
  xlab("# of mutations")+ylab("# of samples")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     axis.text.x=element_text(size=12),legend.position = "none", 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste0(od,"/fig_1e.pdf"), height=4,width=3)



#######################
####################### Fig.1f
#######################
tocurve=function(dd,tag){
  dd2=dd[dd>0]
  xs=c(0);ys=c(0)
  for(ddt in sort(dd2)){
    xs=c(xs,ddt/365); ys=c(ys,leng(dd2<=ddt))
  }
  data.frame(Type=tag,Years=xs,ys)
}

pany=tocurve(ddall[caml1|cmds1|cmpn1],"Any MN")
paml=tocurve(ddAML[caml1],"AML")
pmds=tocurve(ddMDS[cmds1],"MDS")
pmpn=tocurve(ddMPN[cmpn1],"MPN")
pcmml=tocurve(ddMPN[ccmml1],"CMML")

pdat=rbind(pany,paml,pmds,pmpn,pcmml)
pdat$Type=factor(pdat$Type,levels=c("Any MN","AML","MDS","MPN","CMML"))

ggplot(pdat,aes(x=Years,y=ys,color=Type)) +
  geom_line(size=1) + ylab("# of diagnosis")+
  scale_color_manual(values=c(col4,col1,col3,col2,"grey")) + xlim(0,15)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
               axis.text=element_text(size=12), axis.title=element_text(size=14),
               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste0(od,"/fig_1f.pdf"),width=5.5,height=4.0)
  

