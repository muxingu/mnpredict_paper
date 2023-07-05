#######
####### load required codes
#######
source("./model/rBioLib.r")
source("./model/rBaseLib.r")
source("./model/rStatLibBasic.r")
d=read.table("./directory.txt",row.names = 1)["root",]
setwd(d)

#######
####### read data tables (genotype, genotype with vaf, pileup, phenotype and etc)
#######
tgeno = read.table(paste0(d,"/ids_geno.txt"), head=1, row.names=1, stringsAsFactors=F)
tgenovaf = read.table(paste0(d,"/ids_geno_vaf.txt"), head=1, row.names=1, stringsAsFactors=F)
tpile = read.table(paste0(d,"/ids_pileup.txt"), head=1, row.names=1, stringsAsFactors=F)
tpheno = read.table(paste0(d,"/ids_pheno.txt"), head=1, row.names=1, stringsAsFactors=F, sep="\t")
tdd=read.table(paste0(d,"/ids_dd.txt"),head=1,row.names=1,stringsAsFactors=F)
attach(tdd)

#######
####### process first diagnosis
#######
xdiff=function(vs,tddx){
  xmin = apply(tddx,1,min)
  ret=rep(0,length(vs))
  for(i in 1:length(vs)){
    v=vs[i]; m=xmin[i]
    if(v==Inf | m==Inf){ ret[i]=Inf }else{ret[i]=v-m}
  }
  ret
}
xdaml=xdiff(ddAML,tdd[,c("ddMDS","ddMPN","ddCMML")])
xdmds=xdiff(ddMDS,tdd[,c("ddAML","ddMPN","ddCMML")])
xdmpn=xdiff(ddMPN,tdd[,c("ddMDS","ddAML","ddCMML")])


#######
####### define primary AML, MDS, MPN, CMML
#######
caml=ddAML!=Inf
cmds=ddMDS!=Inf
cmpn=ddMPN!=Inf
ccmml=ddCMML!=Inf
cmm=caml|cmds|cmpn|ccmml
# removing undiagnosed MPN cases
mpnexcl2=c("3067940","5189358","3748062","4598030","4587133","5129345","1646410","5167164","2877547","1376313","3892993","4688598","1506238","4907736","1509907","2938638","5321128","1996377","5052179","3123340","5726749","2021376","3774150","2031877","4933974","5253396","5448790","3447518","2986617","4965434","5844932","3782115","5942923","1851999","4222489","2719916","2980786","4784329","3633452","4733997","2797079","5314773","2381520","3308022","4208194","1243259","3117630","1887035","5390460","5701820","5867242","1473651","1133632","3546379","2719133","3500135","3002779","3927688","3229763","2216622","3398597","3979922","5054709","4584176","3073885","2430544","3651953","5962793","2869300","5339695","4656049","4441601","4597246","6023366","1631206","5264741","5458433","3890683","2886146","4117698","5142292","5051372","4809558","4149178","5973371","2674318","1898467","2401709","3100344","3301128","5549443","3032772","5932637","4336363","3423686","5244895","3844201","3891516","1069482","2039563","4362711","4717527","5038012","5128518","4675192","5033582","5381472","1243448")
cmpnexcl2=row.names(tgeno)%in%mpnexcl2
#
caml1=caml & !cmds&!cmpn&!ccmml & ddAML>0 & !(xdaml>-30&xdaml<=0)
cmpn1=cmpn & ddMPN>0 &(ddMPN<ddAML&ddMPN<ddMDS&ddMPN<ddCMML) & !(xdmpn>-30&xdmpn<=0)  & !cmpnexcl2
cmpn1_=cmpn & ddMPN>0 &(ddMPN<ddAML&ddMPN<ddMDS&ddMPN<ddCMML) & !(xdmpn>-30&xdmpn<=0)
ccmml1=ccmml & ddCMML>0 &(ddCMML<ddAML&ddCMML<ddMDS&ddCMML<ddMPN)
cmds1=cmds & ddMDS>0 &(ddMDS<ddAML&ddMDS<ddMPN&ddMDS<ddCMML) & !(xdmds>-30&xdmds<=0) | ccmml1

tddmin=apply(tdd[,1:4],1,min)
camlu = cmm&tddmin>0 & !(caml1|cmds1|cmpn1|ccmml1)
crem1 = caml&ddAML<=0 | cmds&ddMDS<=0 | cmpn&ddMPN<=0
crem2 = (xdaml>-30&xdaml<=0) | (xdmpn>-30&xdmpn<=0) |(xdmds>-30&xdmds<=0) | camlu
crem=crem1 | crem2
ddall=c();txx=tdd[,1:4];  ddall=apply(txx,1,min);tdd[,"ddall"]= ddall
strmds=ifelse(cmds,"MDS","-");strmpn=ifelse(cmpn,"MPN","-");strcmml=ifelse(ccmml,"CMML ","-")
strs=paste0(strmds,",",strmpn,",",strcmml); strs[caml1]="AML1";strs[strs=="-,-,-"]="-"

#######
####### split DMNT3A_R882 and DMNT3A_other
#######
xx=matrix(0,nrow=nrow(tgeno),ncol=ncol(tgeno)); xx[tgeno!='-']=1; row.names(xx)=row.names(tgeno); colnames(xx)=colnames(tgeno)
datfgeno = data.frame(xx)
xx=matrix(0,nrow=nrow(tgeno),ncol=ncol(tpile)); xx[tpile!='-']=1; row.names(xx)=row.names(tpile); colnames(xx)=colnames(tpile)
datfpile=data.frame(xx)
datfgenovaf = tgenovaf
for(cn in colnames(tpile)){
  cvs= strsplit(tpile[,cn],"_"); ncvs=sapply(cvs,length)
  vafs=as.numeric(ifelse(ncvs==3,sapply(cvs,'[',3),0))
  if(cn=="DNMT3A"){
    mv= data.frame(muts=sapply(cvs,'[',1),vafs)
    nv1=mv$vafs; nv1[!grepl("R882",mv$muts)]=0
    nv11=apply(data.frame(nv1,datfgenovaf$DNMT3A_R882),1,max)
    datfgenovaf$DNMT3A_R882=nv11
    nv2=mv$vafs; nv1[grepl("R882",mv$muts)]=0; nv1[grepl("-",mv$muts)]=0
    nv21=apply(data.frame(nv1,datfgenovaf$DNMT3A_other),1,max)
    datfgenovaf$DNMT3A_other=nv21
  }
  if(cn  %in% colnames(datfgenovaf)){
    nv=apply(data.frame(vafs,datfgenovaf[,cn]),1,max)
    datfgenovaf[,cn]=nv
  }
}
vaf=apply(datfgenovaf,1,max)


#######
####### removing linearly dependent blood/biochemistry parameters
#######
tphenona = apply(tpheno,2,function(x)length(x[is.na(x)]))
cnsdrop=c("yob","daac","genetic_ethnic","ever_smoked","smoking_status","ukb_centre","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","AML","MDS","MPN","CMML","dod","censoring","LY_perc","MO_perc","NE_perc","EO_perc","BA_perc","NRBC_perc","RET_perc","HLR_perc","NRBC",
"PCT","LDLD","APOB","HT","MCH","HLR","RBC","MRV","MSCV"
)
cnsblood=c("WBC","RBC","HGB","HT","MCV","MCH","MCHC","RDW","PLT","PCT","MPV","PDW","LY","MO","NE","EO","BA","RET","MRV","MSCV","IRF","HLR")
cnschem=c("ALB","ALP","ALT","APOA","APOB","AST","BILD","BUN","CA","CHOL","CRE","CRP","CYS","GGT","GLU","HBA1C","HDL","IGF1","LDLD","LPA","OES","PHOS","RF","SHBG","TBIL","TES","TP","TRIG","UA","VITD")
cnspat=c("sex","age","bmi")

#######
####### scaling blood/biochemistry parameters
#######
fillscale=function(datf){
  phemed = apply(datf,2,function(x)median(as.numeric(x[!is.na(x)])))
  phesd = apply(datf,2,function(x)sd(as.numeric(x[!is.na(x)])))
  for(i in 1:ncol(datf)){
    val = phemed[i]
    if(!is.na(val)){
      datf[is.na(datf[,i]),i]=val
    }
  }
  for(i in 1:ncol(datf)){
    datf[,i]=(datf[,i]-phemed[i])/phesd[i]
  }
  datf
}
datfpat=tpheno[,colnames(tpheno)%in%cnspat]
datfpat$sex=ifelse(datfpat$sex=="Male",1,0)
datfpat=data.frame(sex=datfpat$sex,fillscale(datfpat[,2:3]))

#######
####### removing NA values
#######
nablood=apply(tpheno[,colnames(tpheno)%in%cnsblood & tphenona<50000],1,function(x)length(x[is.na(x)]))
nachem=apply(tpheno[,colnames(tpheno)%in%cnschem & tphenona<50000],1,function(x)length(x[is.na(x)]))
naphe=nablood+nachem

datfblood=fillscale(tpheno[,colnames(tpheno)%in%cnsblood  &!colnames(tpheno)%in%cnsdrop & tphenona<50000])
datfchem=fillscale(tpheno[,colnames(tpheno)%in%cnschem  &!colnames(tpheno)%in%cnsdrop & tphenona<50000])
datfpheno = data.frame(datfpat, datfblood, datfchem)

tphesel=tpheno[,(colnames(tpheno)%in%c("age","bmi",cnschem,cnsblood)) & tphenona<50000]
phemed = apply(tphesel,2,function(x)median(as.numeric(x[!is.na(x)])))
phesd = apply(tphesel,2,function(x)sd(as.numeric(x[!is.na(x)])))


#######
####### build dataframe of Mutect2 genotype with pileup results
#######
datfgenowp = datfgeno
for(cn in colnames(datfpile)){
  if(cn=="DNMT3A"){
    nv=ifelse(datfgeno$DNMT3A_R882==1,1,ifelse(datfpile[,cn]==1,1,0))
    datfgenowp$DNMT3A_R882=nv
  }else if(cn  %in% colnames(datfgenowp)){
    nv=ifelse(datfgeno[,cn]==1,1,ifelse(datfpile[,cn]==1,1,0))
    datfgenowp[,cn]=nv
  }
}
nmutcol=colSums(datfgenowp)
nmutrow=rowSums(datfgenowp)
cexcl=cmpnexcl2

#######
####### colours for plots
#######
col1=rgb(78/255,98/255,171/255,1)
col2=rgb(195/255,43/255,74/255,1)
col3=rgb(87/255,187/255,171/255,1)
col4=rgb(253/255,213/255,130/255,1)

