######
###### required library
######
library("survival")
library("survminer")

#######
####### load data
#######
source("./model/loaddata.r")

#######
####### read arguments: 1 - type of data; 2-type of MN 
#######
args = commandArgs(trailingOnly=TRUE)
tagty=args[1]; dis=args[2]
vars=read.table(paste0(d,"/cox_stepwise/final_",dis,".txt"),stringsAsFactors=F)[,1]

#######
####### read training set IDs
#######
trainids= as.character(read.table(paste0(d,"/resource/training_ids.txt"),stringsAsFactors=F)[,1])
cbg=(!cmm | caml1|cmds1|cmpn1|ccmml1) & (nablood+nachem)<=2
tartime=c(365,730,1825,3650)
datfs=list(geno=data.frame(datfgenowp,vaf), pheno=datfpheno, genopheno=data.frame(datfgenowp,vaf,datfpheno), 
    genophenocnv=data.frame(datfgenowp,vaf,datfpheno,datfcnv))
tars=list( all=list(all=caml1|cmds1|cmpn1, tag2="all"), 
   aml1=list(aml1=caml1, tag2="AML"), 
   mds1=list(mds1=cmds1, tag2="MDS"), 
   mpn1=list(mpn1=cmpn1,tag="MPN"))


closest=function(vs,tars){
  ret=c()
  for(tar in tars){
    mindis=Inf; minv=-1
    for(i in 1:length(vs)){
      v=vs[i]
      dis=abs(v-tar)
      if(dis<mindis){mindis=dis;minv=i}
    }
    ret=c(ret,minv)
  }
  ret=c(ret,length(vs))
  ret
}

#######
####### construct training set
#######
datf=datfs[[tagty]]; tag=names(tars[[dis]])[1]; tag2=tars[[dis]][[2]]; ctar=tars[[dis]][[1]]
### death as competing risk
evs=ifelse( ctar, 1, ifelse(tdd$ddead<tdd$dcen,2,0) )
od=paste0(d,"/cox/",tag); dir.create(od,showWarnings=F,recursive=T)
survtime = ifelse(tdd[,paste0("dd",tag2)]==Inf, tdd$dcen, ifelse(tdd[,paste0("dd",tag2)]>0, 
  tdd[,paste0("dd",tag2)], 0) )
ctrain=(row.names(datf)%in%trainids) & cbg
datftrain= data.frame(datf[ctrain,vars], time=survtime[ctrain], status=evs[ctrain])

#######
####### Cox regression
#######
cox=coxph( Surv(time,status==1) ~ ., data=datftrain, iter.max=10000)
xx=concordance(cox)
ocon=c(xx$concordance, xx$var, xx$count)
write.table(ocon, file=paste0(od,"/",tagty,"_concordance.txt"),sep="\t",col.names=NA,quote=F)

scox=summary(cox)$coefficients; conf=confint(cox,level=0.9); otabco= data.frame(scox,conf);
colnames(otabco)=c("coef","expcoef","secoef","z","pr","qlow","qhigh")
#### write coefficients
write.table(otabco, file=paste0(od,"/",tagty,"_coeff.txt"),sep="\t",col.names=NA,quote=F)

#######
####### write validation to output
#######
cval=!row.names(datf)%in%trainids & cbg;
datfval = datf[cval,vars]; survtimeval=survtime[cval]
dfpred=matrix(,nrow=0,ncol=length(tartime)+1); colnames(dfpred)=c(tartime,"max")
for(j in seq(1,nrow(datfval),1000)){
  low=j; high=min(j+999,nrow(datfval))
  fit = survfit(cox, datfval[low:high,vars])
  tis=closest(fit$time,tartime)
  pred=t(fit$surv[tis,])
  dfpred=rbind(dfpred,pred)
}
otab=data.frame( round(dfpred,digits=4), survtime=survtimeval)
write.table(otab, file=paste0(od,"/",tagty,"_pred.txt"),sep="\t",col.names=NA,quote=F)

