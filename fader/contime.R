#The exponential gamma model is a continuous two parameter timing model. It answers when
#data vector should be cumulative. Uses kb.csv
growth=function(data,len=length(data)) return(data[2:len]-data[1:(len-1)])
control.eg=function(model,...){
	xy=getxy(model$raw,class(model))
	pop=git(list(...)$pop,sum(xy$y))
	return(list(x=xy$x,y=c(xy$y,pop-sum(xy$y)),mult=c(xy$y,pop-sum(xy$y)),names=c('r','alpha','p'),num=pop,...))
}
ll.eg=function(model,param=NULL,x=model$control$x){
	r=param[1];alpha=param[2]
	x=c(0,x)
	return(log(c(growth(1-(alpha/(alpha+x))^r),(alpha/(alpha+x[length(x)]))^r)))
}
mean.eg=function(model) return(model$param$alpha/(model$param$r-1))
var.eg=function(model) return(model$param$r*model$param$alpha^2/((model$param$r-1)^2*(model$param$r-2)))


#The weibull gamma model is a continuous three parameter timing model. It answers when
#data vector should be cumulative
control.wg=function(model,...){
	xy=getxy(model$raw,class(model))
	pop=git(list(...)$pop,sum(xy$y))
	return(list(x=xy$x,y=c(xy$y,pop-sum(xy$y)),mult=c(xy$y,pop-sum(xy$y)),names=c('r','alpha','c','p'),num=pop,...))
}
ll.wg=function(model,param=NULL,x=model$control$x){
	r=param[1];alpha=param[2];k=param[3]
	x=c(0,x)
	return(log(c(growth(1-(alpha/(alpha+x^k))^r),(alpha/(alpha+x[length(x)]^k))^r)))
}
mean.wg=function(model) return(exp(log(model$param$alpha^(1/model$param$c))+lgamma(1+1/model$param$c)+lgamma(model$param$r-1/model$param$c)-lgamma(model$param$r)))
var.wg=function(model) return(model$param$alpha^(2/model$param$c)/gamma(model$param$r)*(gamma(1+2/model$param$c)*gamma(model$param$r-2/model$param$c)-gamma(1+1/model$param$c)^2*gamma(model$param$r-1/model$param$c)^2/gamma(model$param$r)))

#The gamma gamma model is a continuous 3 parameter timing model. It answers when. 
library(gsl)
control.gg=function(model,...){
	xy=getxy(model$raw,class(model))
	pop=git(list(...)$pop,sum(xy$y))
	return(list(x=xy$x,y=c(xy$y,pop-sum(xy$y)),mult=c(xy$y,pop-sum(xy$y)),names=c('r','alpha','s','p'),num=pop,...))
}
ll.gg=function(model,param=NULL,x=model$control$x){
	r=param[1];alpha=param[2];s=param[3]
	x=c(0,x)
	tmp=1/(s*(gamma(r)*gamma(s)/gamma(r+s)))*(x/(alpha+x))^s*hyperg_2F1(1-r,s,s+1,x/(alpha+x))
	return(log(c(growth(tmp),1-tmp[length(tmp)])))
}
mean.gg=function(model) return(model$param$alpha*model$param$s/(model$param$r-1))
var.gg=function(model) return(model$param$alpha^2*model$param$s*(model$param$r+model$param$s-1)/((model$param$r-1)^2*(model$param$r-2)))
