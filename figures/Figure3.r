library(ggplot2)
library(survival)
library(survminer)

#source("./model/loaddata.r")
od="./figures/fig3"; dir.create(od)

###################
################### Fig.3a
###################
tags=list(aml1=0,mds1=0,mpn1=0)
pd=data.frame(stringsAsFactors=F)
for(tag in names(tags)){
  tsel= read.table(paste0(d,"/cox/",tag,"/genopheno_coeff.txt"),row.names=1,head=1)
  tsel=tsel[!is.na(tsel$coef) & tsel$qlow>-10 & tsel$qhigh<10,]
  tabs=list( geno=tsel[row.names(tsel)%in%c(colnames(datfgeno),"vaf"),], 
             pheno=tsel[row.names(tsel)%in%colnames(datfpheno),])
  for(tyn in names(tabs)){
    tab=tabs[[tyn]]
    pd=rbind(pd, data.frame( Dis=tag, Type=tyn, Gene=row.names(tab), pr=tab$pr, 
                             qlow=tab$qlow, coef=tab$coef, qhigh=tab$qhigh,stringsAsFactors=F))
  }
}

plotforest=function(tab,ll=Inf,hh=Inf){
  agg= aggregate(tab$coef, by=list(Category=tab$Gene), FUN=median)
  sgenes = agg$Category[order(agg$x)];
  if("sex"%in%sgenes){sgenes=c(sgenes[sgenes!="sex"],"sex")}
  if("age"%in%sgenes){sgenes=c(sgenes[sgenes!="age"],"age")}
  if("clone"%in%sgenes){sgenes=c(sgenes[sgenes!="clone"],"clone")}
  if(ll==Inf){ll=floor(min(tab$qlow)/5)*5} 
  if(hh==Inf){hh=ceiling(max(tab$qhigh)/5)*5}
  plot(0, type="n", xlim=c(ll,hh), ylim=c(0,length(sgenes)),xlab="",ylab="",yaxt="n",bty="n"); 
  par(lwd=2,las=1,mar = c(2, 5, 2, 2))
  for(i in seq(1,length(sgenes),2)){ rect(ll,i-0.5,hh,i+0.5, col="grey", border=F)}
  axis(2, at=1:length(sgenes), labels=sgenes)
  abline(v=0, lty=2, lwd=1.5, col="black")
  for(i in 1:length(sgenes)){
    sgene=sgenes[i]; stab=tab[tab$Gene==sgene,]; stab=stab[order(stab$Dis),]
    medrow=median(1:nrow(stab))
    for(j in 1:nrow(stab)){
      diff = (medrow-j)*nrow(stab)/12; dis=stab$Dis[j]
      v=stab$coef[j]; l=stab$qlow[j]; h=stab$qhigh[j]
      xd=0.005*(hh-ll); yd=0.07
      col=ifelse(dis=="aml1",col1,ifelse(dis=="mds1",col3,col2))
      plotbar(i+diff,v,xd,yd,l,h,col)
    }
  }
}
plotbar=function(i,v,xd,yd,l,h,col){
  rect(v-xd,i-yd, v+xd,i+yd, col=col, border=col)
  segments(l,i, h,i,col=col)
}
plgenes=c("DNMT3A_R882","DNMT3A_other","TET2","JAK2","ASXL1",
          "IDH2","SRSF2","SF3B1","CALR","MPL","TP53","U2AF1",
          "vaf")
pdf( paste0(od,"/fig_3a_geno.pdf"), width=6.5, height=6.5); 
plotforest(pd[pd$Type=="geno" & pd$Gene%in%plgenes,], ll=-4,hh=6)
dev.off()

pdf( paste0(od,"/fig_3a_pheno.pdf"), width=6.5, height=7.5)
plotforest(pd[pd$Type=="pheno",],ll=-1,hh=1)
dev.off()


##########################
########################## Fig.3b-g
##########################
tars=list(
  list(aml1=caml1, tag2="AML", plgenes=c("IDH2","SRSF2")), 
  list(mds1=cmds1, tag2="MDS", plgenes=c("SF3B1","SRSF2")), 
  list(mpn1=cmpn1,tag="MPN", plgenes=c("JAK2","CALR")))
datf=data.frame(datfgenowp,datfpheno)
cbg=(!cmm | caml1|cmds1|cmpn1) & (nablood+nachem)<=2 & !cexcl

for(i in 1:length(tars)){
  tag=names(tars[[i]])[1]; tag2=tars[[i]][[2]]; ctar=tars[[i]][[1]]
  survtime = ifelse(tdd[,paste0("dd",tag2)]==-1, tdd$dcen, ifelse(tdd[,paste0("dd",tag2)]>0, tdd[,paste0("dd",tag2)], 0) )
  evs=ifelse( ctar, 1, 0); evtag="_comp"
  dft= data.frame(datf, time=survtime, status=evs)
  gns=tars[[i]][[3]]
  for(gn in gns){
    vaf2=datfgenovaf[,gn]
    vafcat=ifelse(vaf2>0.2,2,ifelse(vaf2>0,1,0))
    dft2=data.frame(dft);dft2[,"vaf"]=vafcat; dft2=dft2[cbg,]
    fit= survfit( as.formula(paste0("Surv(time,status) ~ vaf")), data=dft2 )
    ggsurvplot(fit, data=dft2, palette = c(rgb(0.6,0.6,0.6,1),col4,col1) )
    ggsave(filename=paste0(od,"/fig3_kaplan_",tag,"_",gn,".pdf"), width=3.5, height=3)
  }
}
