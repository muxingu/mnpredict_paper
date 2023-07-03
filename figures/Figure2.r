library(gplots)
library(ggplot2)
library(dplyr)
library(reshape)
library(GenVisR)

#source("./model/loaddata.r")
od="./figures/fig2"; dir.create(od)

tars=list(AML=caml1,MDS=cmds1,MPN=cmpn1_)
datf=datfgenowp
cbg=T; bgtag="bg-t"
targenes=c("DNMT3A_other","DNMT3A_R882","ASXL1","TP53","TET2",
           "JAK2","CALR","MPL","SRSF2","U2AF1","SF3B1",
           "IDH2")


##################### 
##################### Fig.2a
#####################

nmm=leng(cmm); nh=leng(!cmm)
mat=matrix(,ncol=3,nrow=0)
for(tg in targenes){
  ncmm = leng(datfgenowp[,tg]==1 & cmm)/nmm *100
  nheal = leng(datfgenowp[,tg]==1 & !cmm)/nh *100
  mat=rbind(mat,c(Type="Pre-MN",Gene=tg,ncmm))
  mat=rbind(mat,c(Type="Control",Gene=tg,nheal))
}
pd=data.frame(Type=mat[,1],Gene=mat[,2],Number=as.numeric(mat[,3]),stringsAsFactors=T)
pd$Number=as.numeric(pd$Number); pd$Gene=factor(pd$Gene,levels=targenes)
pd$Type=factor(pd$Type,levels=c("Pre-MN","Control"))
pd$Number[pd$Type=="Control"]= -pd$Number[pd$Type=="Control"]

ggplot(pd, aes(fill=Type,x=Gene,y=Number)) +
  geom_bar(position="identity",stat="identity", width=.7) + 
  scale_fill_manual(values=c(col1,col2)) + 
  ylab("% of cases") + xlab("")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     axis.text.x=element_text(angle=45,vjust=1,hjust=1),
                     legend.key.size = unit(0.3, 'cm'),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste0(od,"/fig_2a.pdf"),height=3,width=4.5)


#####################
##################### Fig.2c
#####################
pd=data.frame(); pdp=data.frame()
for(tag in names(tars)){
  ctar= tars[[tag]]
  vs=c(); vsp=c()
  for(gene in targenes){
    fish=sf(cbg,ctar,datf[,gene]==1,ifprint=F)
    ovl=fish$exp[4]; exp=fish$exp[5]; pv=fish$pv
    #v=ifelse(exp==0,0,ifelse(ovl<5,0,log10(ovl/exp)))
    v=ifelse(exp==0|ovl==0,0,log10(ovl/exp))
    vs=c(vs,v); vsp=c(vsp,pv)
    
  }
  pd=rbind(pd,data.frame(Type=tag,Gene=targenes,OR=vs))
  pdp=rbind(pdp,data.frame(Type=tag,Gene=targenes,p=vsp))
}
pd$Gene=factor(pd$Gene,levels=targenes)
pd$Type=factor(pd$Type,levels=c("MPN","MDS","AML"))

ggplot(pd, aes(Gene,Type)) + geom_tile(aes(fill=OR)) + 
  scale_fill_gradient(low="white",high=col3) + 
  ylab("")+xlab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  theme(legend.key.size = unit(0.3, 'cm'))
ggsave(paste0(od,"/fig_2c.pdf"),width=3.5,height=1.6)  


######################
###################### Fig.2d
######################
age=tpheno$age
vaf=apply(tgenovaf,1,max)
tars=list(list(age<50,"<50"),
          list(age>=50&age<60,"50-59"),
          list(age>=60&age<65,"60-65"),
          list(age>=65,">65")
)
groups=list(AML=caml1, MDS=cmds1, MPN=cmpn1_)
cbg=T
pd=data.frame()
for(i in 1:length(tars)){
  ctar=tars[[i]][[1]]; tag=tars[[i]][[2]]
  for(gn in names(groups)){
    dft=data.frame(Age=tag,Group=gn,VAF=vaf[ctar&cbg&groups[[gn]]])
    pd=rbind(pd,dft)
  }
}
pd$Age=factor(pd$Age, levels=c(tars[[1]][[2]],tars[[2]][[2]],tars[[3]][[2]],tars[[4]][[2]]))
pd$Group=factor(pd$Group, levels=names(groups))
ns=c()
for(gn in levels(pd$Group)){for(age in levels(pd$Age)){ns=c(ns,nrow(pd[pd$Age==age&pd$Group==gn,]))}}
xx = pd %>% group_by(Group,Age) %>% summarise(y = 0.48, VAF = first(VAF))

ggplot(pd, aes(x=Group,fill=Age,y=VAF)) + ylim(0,0.5)+
  geom_boxplot(outlier.shape=NA) + 
  scale_fill_manual(values=c(col1,col2,col3,col4)) + 
  xlab("")+ ylab("Largest clone VAF")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave(file=paste0(od,"/fig_2d.pdf"), width=4.5, height=3)



###################
################### Fig.2b
###################
tars=list(AML=caml1,MPN=cmpn1, MDS=cmds1)
targenes=c("DNMT3A_R882", "DNMT3A_other","TET2",
           "JAK2","ASXL1","IDH2","SRSF2","SF3B1",
           "CALR","MPL","TP53","U2AF1")
cols = c(col1,col1,"white")

for(tag in names(tars)){
  ctar=tars[[tag]]
  mdat=melt(data.matrix(datfgeno[ctar,])); mdat=mdat[mdat$value!=0,]; colnames(mdat)=c("sample","gene","variant_class");mdat$variant_class="Mutect2"
  datfpile2 = datfpile[ctar,]; colns=colnames(datfpile2); colns[colns=="DNMT3A"]="DNMT3A_R882"; colnames(datfpile2)=colns
  mdatpile=melt(data.matrix(datfpile2)); mdatpile=mdatpile[mdatpile$value!=0,]
  addat=mdatpile[!paste0(mdatpile[,1],"_",mdatpile[,2])%in%paste0(mdat[,1],"_",mdat[,2]),];colnames(addat)=colnames(mdat)=c("sample","gene","variant_class");addat$variant_class="Pileup"
  
  x=rbind(mdat,addat)
  gs=as.character(x$gene); gs[!gs%in%targenes]="Other"
  mdat2=data.frame(sample=x$sample,gene=gs,variant_class=x$variant_class)
  ord=c(targenes,"Other")
  
  tid=mdat2[1,1]
  for(tgn in targenes[!targenes%in%unique(mdat2$gene)]){
    mdat2=rbind(mdat2,c(tid,tgn,"empty"))
  }
  
  pdf(paste0(od,"/fig_2b_",tag,".pdf"),height=9,width=20)
  waterfall(mdat2, fileType = "Custom", variant_class_order=c("Mutect2","Pileup","empty"),
            mainPalette = cols, geneOrder=ord, clinVarCol=cols, 
            main_geneLabSize=18)
  dev.off()
}




