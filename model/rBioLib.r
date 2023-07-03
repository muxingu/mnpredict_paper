

file_geneinfo_human="C:/muxingu/genome/human/grch38/ass/gene-name_Homo_sapiens.GRCh38.93.txt"

humanPC=function(){
	PCs(file_geneinfo_human)
}

humanAnnos=function(){
	Annos(file_geneinfo_human)
}

Annos=function(path){
	read.table(path,row.names=1,head=1,stringsAsFactors=F)
}

PCs=function(path){
	t= read.table(path,row.names=1,head=1,stringsAsFactors=F)
	tpc = t[t$Type=="protein_coding",]
	ncou = table(tpc$Name)
	inc = names(ncou)[ncou==1]
	tpc[tpc$Name %in% inc,]
}

tidyDESeq=function(dat){
	res=dat
	res$log2FoldChange[is.nan(res$log2FoldChange)]=0; res$log2FoldChange[is.na(res$log2FoldChange)]=0; 
	res$lfcSE[is.nan(res$lfcSE)]=0; res$lfcSE[is.na(res$lfcSE)]=0; 
	res$stat[is.nan(res$stat)]=0; res$stat[is.na(res$stat)]=0; 
	res$pvalue[is.nan(res$pvalue)]=1; res$pvalue[is.na(res$pvalue)]=1; 
	res$padj[is.nan(res$padj)]=1; res$padj[is.na(res$padj)]=1;
	res$baseMean=round(res$baseMean,digits=1)
	res$log2FoldChange=round(res$log2FoldChange,digits=2)
	res$stat=round(res$stat,digits=1)
	res$pvalue = signif(res$pvalue, digits=3)
	res$padj = signif(res$padj, digits=3)
	res
}

