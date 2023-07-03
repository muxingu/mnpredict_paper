library(ggplot2)

#source("./model/loaddata.r")
od="./figures/fig5"; dir.create(od)
tartime=c(seq(100,4500,100))


###### functions
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
getcox = function(datf=data.frame(datfgenowp,datfpheno,vaf), ctar, tag2, vars){
  evs=ifelse( ctar, 1, ifelse(tdd$ddead<tdd$dcen,2,0) )
  tv=tdd[,paste0("dd",tag2)]
  survtime = ifelse(tv==Inf, tdd$dcen, ifelse(tv>0, tv, 0) )
  ctrain=(row.names(datf)%in%trainids) & cbg 
  datftrain= data.frame(datf[ctrain,vars], time=survtime[ctrain], status=evs[ctrain])
  cox=coxph( Surv(time,status==1) ~ ., data=datftrain, iter.max=10000)
  cox
}
customplot = function(samp){
  sfaml=survfit(coxaml, samp); ciaml=closest(sfaml$time,tartime); svaml=sfaml$surv[ciaml]
  sfmds = survfit(coxmds, samp);cimds=closest(sfmds$time,tartime); svmds=sfmds$surv[cimds]
  sfmpn = survfit(coxmpn, samp);cimpn=closest(sfmpn$time,tartime); svmpn=sfmpn$surv[cimpn]
  svdat=data.frame(svaml,svmds,svmpn)
  svnorm= matrix(,nrow=0,ncol=5)
  svnorm=c(0,1,1,1,1)
  for(i in 1:nrow(svdat)){
    x=as.numeric(as.vector(svdat[i,])) 
    free=prod(x); xsum=sum(1-x); ds=(1-free)*(1-x)/xsum
    vs=c(0,free,free+ds[1], free+ds[1]+ds[2],1)
    svnorm= rbind(svnorm,vs)
  }
  pdat=data.frame(svnorm,x=c(0,sfaml$time[ciaml]/365))
  pdat2=matrix(,nrow=0,ncol=4)
  for(i in 1:nrow(pdat)){
    pdat2=rbind(pdat2, c(pdat$x[i],pdat[i,1],pdat[i,2],"MN-free"))
    pdat2=rbind(pdat2, c(pdat$x[i],pdat[i,2],pdat[i,3],"AML"))
    pdat2=rbind(pdat2, c(pdat$x[i],pdat[i,3],pdat[i,4],"MDS"))
    pdat2=rbind(pdat2, c(pdat$x[i],pdat[i,4],pdat[i,5],"MPN"))
  }
  pdat2=data.frame(x=as.numeric(pdat2[,1]),min=as.numeric(pdat2[,2]),max=as.numeric(pdat2[,3]),Type=pdat2[,4])
  pdat2$Type=factor(pdat2$Type,levels=c("MN-free","AML","MDS","MPN"))
  c0=rgb(0,0,0,0.3);c1=rgb(253/255,213/255,130/255,0.7);c2=rgb(225/255,43/255,54/255,0.4);c3=rgb(87/255,187/255,171/255,0.4)
  ggplot(pdat2, aes(x=x)) +
    geom_ribbon(aes(ymin=min,ymax=max,fill=Type))+
    scale_fill_manual(values=c(c0,c2,c1,c3))+
    #geom_line(aes(y=X2),size=1)+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       axis.text=element_text(size=12),axis.title=element_text(size=14),
                       axis.line = element_line(colour = "black"))+
    xlab("Time / years")+ylab("Probabilby")+
    scale_x_continuous(limits = c(0.3,15), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0,0.05)))
}

######################
###################### load cox models
######################
trainids= as.character(read.table(paste0(d,"/resource/training_ids.txt"),stringsAsFactors=F)[,1])
cbg=(!cmm | caml1|cmds1|cmpn1) & (nablood+nachem)<=2 & !cexcl
varsaml=read.table(paste0(d,"/cox_stepwise/final_aml1.txt"),stringsAsFactors=F)[,1]
varsmds=read.table(paste0(d,"/cox_stepwise/final_mds1.txt"),stringsAsFactors=F)[,1]
varsmpn=read.table(paste0(d,"/cox_stepwise/final_mpn1.txt"),stringsAsFactors=F)[,1]
coxaml = getcox(ctar=caml1,tag2="AML",vars=varsaml) 
coxmds = getcox(ctar=cmds1,tag2="MDS",vars=varsmds)
coxmpn = getcox(ctar=cmpn1,tag2="MPN",vars=varsmpn)


######################
###################### Fig.5
######################
tarsamps=list(aml1="3849974", mds1="2529724", mpn1="5437196")
datf=data.frame(datfgenowp,datfpheno,vaf)
for(dis in names(tarsamps)){
  sampn=tarsamps[[dis]]
  customplot(datf[sampn,])
  ggsave(paste0(od,"/fig5_MN-predict_",dis,".pdf"), width=5.5, height=3.5)
  vars=read.table(paste0(d,"/cox_stepwise/final_",dis,".txt"))[,1]
  pdats=list(
    blood=data.frame(t(datf[sampn,intersect(vars,colnames(datfblood)),drop=F])),
    chem=data.frame(t(datf[sampn,intersect(vars,colnames(datfchem)),drop=F]))
    )
  for(pdatn in names(pdats)){
    pdat=pdats[[pdatn]]; colnames(pdat)="val";
    pdat$par=row.names(pdat); pdat=pdat[order(pdat$val,decreasing=T),]
    pdat$par=factor(pdat$par,levels=pdat$par)
    ggplot(pdat,aes(x=par,y=val)) + 
      geom_bar(stat="identity", fill=col1,width=0.8) + xlab("")+ylab("Normalized value")+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
           axis.text=element_text(size=12), axis.title=element_text(size=14),
           axis.text.x=element_text(angle=45,hjust=1),
           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    w=1+nrow(pdat)/5
    ggsave(paste0(od,"/fig5_",pdatn,"_",dis,".pdf"), width=w, height=3.5)
  }
}

