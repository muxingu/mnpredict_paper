######
###### Extract final parameters from stepwise Cox regression
######
args = commandArgs(trailingOnly=TRUE)
tar=args[1]

indir="./cox_stepwise"
init=read.table("./resource/init.txt",stringsAsFactors=F)[,1]
genes=c("ASXL1","BCOR","BCORL1","BRAF","CALR","CBL","CSF1R","DNMT3A_R882","DNMT3A_other","ETV6","EZH2","FLT3","GNAS","GNB1","IDH1","IDH2","JAK2","KDM6A","KIT","KRAS","MPL","MYD88","NPM1","NRAS","PHF6","PIGA","PPM1D","RAD21","RUNX1","SF1","SF3B1","SMC1A","SMC3","SRSF2","STAG2","TET2","TP53","U2AF1","ZRSR2")

if(tar=="aml1"){tard="AML";excl=c("JAK2","CALR","MPL","TP53","U2AF1")}
if(tar=="mds1"){tard="MDS";excl=c("CALR","MPL")}
if(tar=="mpn1"){tard="MPN";excl=c("TP53","U2AF1","MPL")}
t=read.table(paste0(indir,"/",tard,"/summary.txt"),head=1,row.names=1)
vmax=max(t$Conc); subns=c(row.names(t)[1]); vs=t$Conc
for(i in 2:nrow(t)){
incr=(vs[i]-vs[i-1])/vmax *100
if(incr<0.1){break;}else{
 if(!row.names(t)[i]%in%genes){
  subns=c(subns,row.names(t)[i])
 }
}
}
init2=init[!init%in%excl]
ons=c(subns,init2)
write.table(ons,file=paste0(indir,"/final_",tar,".txt"),quote=F,row.names=F,col.names=F)

