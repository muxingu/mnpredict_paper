qplotn= 20


ttest= function(avs, bvs){
  if(var(avs)==0 | var(bvs)==0){
    1
  }else{
	  t.test( avs[cvalid(avs)], bvs[cvalid(bvs) ] )[[3]]
  }
}

ttest_pvs= function(df){
  apply( df,1,function(vs) ttest(vs[1:(length(vs)/2)], vs[(length(vs)/2+1):length(vs)]))
}



sumfisher= function(ctotal, cta, ctb, ifprint=T){
  if(length(ctotal)==1){
    ctot= rep(ctotal, length(cta))
  }else{
    ctot=ctotal
  }
	ctara= cta&ctot
	ctarb= ctb&ctot
	la= leng(ctara)
	lb= leng(ctarb)
	lab= leng(ctara&ctarb)
	ltot= leng(ctot)
	
	mat= matrix( c(la-lab, lab, ltot-lb-(la-lab), lb-lab) , ncol=2)
	fish = fisher.test(mat)
	pv= fish[[1]]
	expec=as.double(la)*as.double(lb)/ltot
	allvs= c(lab, expec, pv)
	if(ifprint){
		print( paste("tot   ",deparse(substitute(cta)),"   ",deparse(substitute(ctb)),"    ovl    exp"),sep="")
		print("")
		print (paste( ltot, la,lb , lab, expec, sep="    ") )
		print("")
		print(pv)
		print(allvs)
		print("")
	}
	list(mat= mat, fish=fish, pv=pv, exp= c(ltot, la,lb , lab, round(expec,digits=3)))
}
sf=sumfisher


bplot= function(tar,lis, cc=T, border="black", sub="",ylab="", col="white", notch=F, main=""){
	par(lwd=2)
	lens= c()
	lvs= list()
	for(i in 1:length(lis)){
		lv= tar[ lis[[i]] & cc ]
		lens= c(lens, length(lv))
		lvv= lv[lv!=Inf & lv!=-Inf & !is.na(lv) & !is.nan(lv)]
		lvs= c(lvs, list(lvv))
	}
	pv= t.test(lvs[[1]], lvs[[2]])[[3]]
	boxplot(lvs, names= lens, outline=F, main=paste0(main,"\n",pv), border=border, sub=sub, ylab=ylab, col=col, notch=notch)
	axis(2, lwd=2)
	#lvs
}

bplot2= function(lis, tar, border="black"){
	par(lwd=2)
	lens= c()
	lvs= list()
	for(i in 1:length(lis)){
		lv= lis[[i]][tar ]
		lvv= lv[lv!=Inf & lv!=-Inf & !is.na(lv) & !is.nan(lv) ]
		lvs= c(lvs, list(lvv))
		lens= c(lens, length(lvv))
	}
	pv= t.test(lvs[[1]], lvs[[2]])[[3]]
	boxplot(lvs, names= lens, outline=F, main=pv, border= border)
	lvs
}
quant= function(vs, n=10){
	nvs= vs[!is.nan(vs) & !is.na(vs)]
	quantile(nvs, probs=seq(0,1, 1/n))
}

lratio= function(vas, vbs, minv){
  van= vas; van[!cnz(vas)] <- 0
  vbn= vbs; vbn[!cnz(vbs)] <- 0
  log2( (vas+minv)/(vbs+minv) )
}


qlist= function(tara, tar, qs){	
	li= list()
	for(i in 1: (length(qs)-1) ){
		ll= qs[[i]]
		hl= qs[[i+1]]
		c= tara>=ll & tara<hl
		li= c(li, list(tar[c]))
	}
	li
}

qplot = function(tara, tar, cc=T, qplotn=10){
	tara= tara[cc]
	tar= tar[cc]
	qs= quantile(tara,  na.rm=T, probs= seq(0,1,1/qplotn))
	qs=c(-Inf, qs)
	li= qlist(tara, tar, qs)
	lens= c()
	for(i in 1:length(li)){
		lens= c(lens, length(li[[i]]))
	}
	lens= qs[1:(length(qs)-1)]
	w=2
	par(lwd=w)
	boxplot(li, outline=F, lwd=w, frame=F, names=lens)
	axis(side=2, label=F, lwd= w*0.8)
}

ftest= function(cbg, ca, cb){
	no= leng(ca&cb)
	mat = matrix( 
		c( no, leng(ca)-no, leng(cb)-no, leng(cbg)-leng(ca)-leng(cb)+no ),
		ncol=2
	)
	fisher.test(mat)[[1]]
}

getexp= function(cbg, cta, ctb){
	na= leng(cta&cbg)
	nb= leng(ctb&cbg)
	nn= leng(cbg)
	c( leng(cta&ctb&cbg), na*nb/nn )
}

mediannz= function(vs){
	median( vs[!is.na(vs)&!is.nan(vs)]) 
}

sdnz= function(vs){
	sd( vs[!is.na(vs)&!is.nan(vs)]) 
}

upperquartile= function(vs){
	qs = quantile(vs, probs=c(0,0.25,0.5,0.75,1))
	qs[[4]]
}

uqnorm = function(vs, ifnz=true){
	vsnz = vs[vs>0]
	vs / upperquartile(vsnz)
}

mediannorm = function(vs, ifnz=true){
	vsnz = vs[vs>0]
	vs / median(vsnz)
}

rnaseqnorm = function(comdat, lens){
	fpkm = sweep( comdat, 2, colSums(comdat)/1000000, "/" )
	fpkm = sweep( fpkm, 1, lens/1000, "/" )
	uqs =c()
	for(i in 1:ncol(fpkm)){
		quant= quantile(fpkm[fpkm[,i]>0,i])
		uqs=c(uqs,quant[[4]])
	}
	uqnorm = sweep( fpkm, 2, uqs/10, "/" )
	list(fpkm=fpkm, uq=uqnorm)
}

vnz = function(vs){
	vs[vs!=0]
}


themen = function(){
	blank() + theme( legend.position="none", axis.text=b(), axis.title = b(),
			 axis.line=element_line(size=.8), axis.ticks=element_line(size=1), axis.ticks.length=unit(0.25,"cm") )
}
blank = function(){
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
}
b = function(){element_blank()}

lm_eqn <- function(df){
	m <- lm(y ~ x, df)
	b = format(coef(m)[1], digits = 2)
	k = format(coef(m)[2], digits = 2)
	r2 = format(summary(m)$r.squared, digits = 3)
	paste0( "y = ",k,"x + ",b,"  (R2=",r2,")")              
}



