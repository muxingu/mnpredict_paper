
charAt = function(str, n){
  substr(str,n,n)
}


wtab= function(vs, fn= "", varn= deparse(substitute(vs)), wd=getwd(), rown=T, coln=T){
	setwd(wd)
	varn= gsub("\\[", "_", varn)
	varn= gsub("\\]", "_", varn)
	varn= gsub("\\ ", "", varn)
	varn= gsub("\\|", "_", varn)
	varn= gsub("\\&", "_", varn)
	varn= gsub(">", "more", varn)
	varn= gsub("<", "less", varn)
	if(fn==""){
		fn= paste(varn,".txt", sep="")
	}
	write.table( vs, file=fn, row.names=rown, col.names=coln, quote=F, sep="\t")
}

leng= function(cs){
	cx= !is.na(cs) & !is.nan(cs)
	length(cs[cs&cx])
}

cvalid= function(va){
  !is.na(va) & !is.nan(va) & va!=Inf & va!=-Inf 
}

cnz= function(va){
	cvalid(va) & va!=0
}

cornz= function(va, vb, method= "pearson"){
	cc= cnz(va) & cnz(vb)
	cor(va[cc], vb[cc], method=method)
}


catn= function(vs){
	cat(vs, sep="\n")
}


inrange= function(vec, min=-Inf, max=Inf){
  vec > min & vec< max
}

lt = function(vs, base=10, off=1){
	return( log(off+vs) / log(base) )
}



