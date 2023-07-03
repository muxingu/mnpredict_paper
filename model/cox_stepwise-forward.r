######
###### required library
######
library("survival")
library("survminer")
library(pROC)

######
###### load data
######
source("./model/loaddata.r")
######
###### read arguments
######
args = commandArgs(trailingOnly=TRUE)
dis=args[1]; tag2=args[2]

######
###### read training set IDs
######
trainids= as.character(read.table(paste0(d,"/resource/training_ids.txt"),stringsAsFactors=F)[,1])
cbg=(!cmm | caml1|cmds1|cmpn1|ccmml1 ) & (nablood+nachem)<=2 & !cexcl

######
###### construct training dataframe and outcome
######
datf=data.frame(datfgenowp,datfpheno,vaf); vars=colnames(datf)
ctar=get(paste0("c",dis));tv=tdd[,paste0("dd",tag2)]
evs=ifelse( ctar, 1, ifelse(tdd$ddead<tdd$dcen,2,0) )
survtime = ifelse(tv==Inf, tdd$dcen, ifelse(tv>0, tv, 0) )
ctrain=(row.names(datf)%in%trainids) & cbg

od=paste0(d,"/cox_stepwise-forward/",tag2); dir.create(od,showWarnings=F,recursive=T)
od2=paste0(od,"/iter"); dir.create(od2,showWarnings=F,recursive=T)

######
###### stepwise regression
######
excl=read.table(paste0(d,"/resource/init.txt"),stringsAsFactors=F)[,1]
omat=matrix(,nrow=0,ncol=8); i=0
while(length(excl)<(length(vars))){
  i=i+1
  itermat= matrix(,nrow=0,ncol=8)
  bestconc=-1; bestxx=-1; bestgn=""
  for(gn in vars[!vars%in%excl]){
    exclt=c(excl,gn)
    datftrain= data.frame(datf[ctrain,vars[vars%in%exclt]], time=survtime[ctrain], status=evs[ctrain])
    cox=coxph( Surv(time,status==1) ~ ., data=datftrain, iter.max=10000)
    xx=concordance(cox); conc=xx$concordance
    itermat=rbind(itermat, c(gn,conc,xx$var,xx$count))
    if(conc>bestconc){bestconc=conc; bestxx=xx; bestgn=gn}
  }
  omat=rbind(omat, c(bestgn,bestxx$concordance,bestxx$var,bestxx$count) )
  excl=c(excl,bestgn)
  colnames(itermat)=c("Gene","Conc","Var","concordant","discordant","tied.x","tied.y","tied.xy")
  write.table(itermat,file=paste0(od2,"/",i,"_",bestgn,".txt"),row.names=F,sep="\t",quote=F)
}
colnames(omat)=c("Gene","Conc","Var","concordant","discordant","tied.x","tied.y","tied.xy")
write.table(omat,file=paste0(od,"/summary.txt"),row.names=F,sep="\t",quote=F)

