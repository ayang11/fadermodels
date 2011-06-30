growth=function(data,len=length(data)) return(data[2:len]-data[1:(len-1)])

#The exponential model is a continuous one parameter timing model. It answers when
setClass('exp',contains='fm')
control.exp=function(model,...){
	xy=getxy(model@raw,class(model))
	pop=git(list(...)$pop,sum(xy$y))
	return(list(x=c(0,xy$x),y=c(xy$y,pop-sum(xy$y)),names=c('lambda'),num=pop,...))
}
ll.exp=function(model,param=NULL,x=model@control$x){
	return(log(c(growth(1-exp(-param*x)),exp(-param*x[length(x)]))))
}
mean.exp=function(model) return(1/model@param$lambda)
vcov.exp=function(model) return(1/model@param$lambda^2)
paramplot.exp=function(model,...) {
	par(mfrow=c(1,1))
	plotSpike(model@param$lambda,model@param$p)
}

#The exponential gamma model is a continuous two parameter timing model. It answers when
setClass('eg',contains='fm')
control.eg=function(model,...){
	xy=getxy(model@raw,class(model))
	pop=git(list(...)$pop,sum(xy$y))
	return(list(x=c(0,xy$x),y=c(xy$y,pop-sum(xy$y)),names=c('r','alpha'),num=pop,...))
}
ll.eg=function(model,param=NULL,x=model@control$x){
	r=param[1];alpha=param[2]
	return(log(c(growth(1-(alpha/(alpha+x))^r),(alpha/(alpha+x[length(x)]))^r)))
}
mean.eg=function(model) return(model@param$alpha/(model@param$r-1))
vcov.eg=function(model) return(model@param$r*model@param$alpha^2/((model@param$r-1)^2*(model@param$r-2)))
paramplot.eg=function(model,...) {
	par(mfrow=c(1,1))
	plotGamma(model@param$r,model@param$alpha,model@param$p)
}


#The weibull gamma model is a continuous three parameter timing model. It answers when
setClass('wg',contains='fm')
control.wg=function(model,...){
	xy=getxy(model@raw,class(model))
	pop=git(list(...)$pop,sum(xy$y))
	return(list(x=c(0,xy$x),y=c(xy$y,pop-sum(xy$y)),names=c('r','alpha','c'),num=pop,...))
}
ll.wg=function(model,param=NULL,x=model@control$x){
	r=param[1];alpha=param[2];k=param[3]
	return(log(c(growth(1-(alpha/(alpha+x^k))^r),(alpha/(alpha+x[length(x)]^k))^r)))
}
mean.wg=function(model) return(exp(log(model@param$alpha^(1/model@param$c))+lgamma(1+1/model@param$c)+lgamma(model@param$r-1/model@param$c)-lgamma(model@param$r)))
vcov.wg=function(model) return(model@param$alpha^(2/model@param$c)/gamma(model@param$r)*(gamma(1+2/model@param$c)*gamma(model@param$r-2/model@param$c)-gamma(1+1/model@param$c)^2*gamma(model@param$r-1/model@param$c)^2/gamma(model@param$r)))
paramplot.wg=function(model,...) {
	par(mfrow=c(1,1))
	plotGamma(model@param$r,model@param$alpha,model@param$p)
}

#The gamma gamma model is a continuous 3 parameter timing model. It answers when. 
setClass('gg',contains='fm')
control.gg=function(model,...){
	xy=getxy(model@raw,class(model))
	pop=git(list(...)$pop,sum(xy$y))
	return(list(x=c(0,xy$x),y=c(xy$y,pop-sum(xy$y)),names=c('r','alpha','s'),num=pop,...))
}
ll.gg=function(model,param=NULL,x=model@control$x){
	r=param[1];alpha=param[2];s=param[3]
	tmp=1/(s*(gamma(r)*gamma(s)/gamma(r+s)))*(x/(alpha+x))^s*hyperg_2F1(1-r,s,s+1,x/(alpha+x))
	return(log(c(growth(tmp),1-tmp[length(tmp)])))
}
mean.gg=function(model) return(model@param$alpha*model@param$s/(model@param$r-1))
vcov.gg=function(model) return(model@param$alpha^2*model@param$s*(model@param$r+model@param$s-1)/((model@param$r-1)^2*(model@param$r-2)))
paramplot.gg=function(model,...) {
	par(mfrow=c(1,1))
	plotGamma(model@param$r,model@param$alpha,model@param$p)
}