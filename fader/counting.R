#The poisson model is a discrete one parameter counting model. It answers how many
setClass('pois',contains='fm')
control.pois=function(model,t=1,...){
	xy=getxy(model@raw,class(model))
	return(list(x=xy$x,y=xy$y,num=sum(xy$y),t=t,names=c('lambda'),allowspike=TRUE,...))
}
ll.pois=function(model,param=NULL,t=model@control$t,x=model@control$x){
	return(log(dpois(x,param*t)))
}
mean.pois=function(model,t=model@control$t) return(sum(model@param$lambda*t*model@param$p))
vcov.pois=function(model,t=model@control$t) return(model@param$lambda*t)
paramplot.pois=function(model,...) {
	par(mfrow=c(1,1))
	plotSpike(c(0,model@param$lambda),c(1-sum(model@param$p),model@param$p))
}


#The nbd model is a discrete two parameter counting model. It answers how many
setClass('nbd',contains='fm')
control.nbd=function(model,t=1,...){
	xy=getxy(model@raw,class(model))
	return(list(x=xy$x,y=xy$y,num=sum(xy$y),t=t,names=c('r','alpha'),allowspike=TRUE,...))
}
ll.nbd=function(model,param=NULL,t=model@control$t,x=model@control$x){
	r=param[1];alpha=param[2]
	return(lgamma(r+x)-lgamma(r)-log(factorial(x))+r*(log(alpha/(alpha+t)))+x*(log(t/(alpha+t))))
}
mean.nbd=function(model,t=model@control$t) return(model@param$r*t/model@param$alpha)
vcov.nbd=function(model,t=model@control$t) return(model@param$r*t/model@param$alpha+model@param$r*t^2/model@param$alpha^2)
paramplot.nbd=function(model,...) {
	par(mfrow=c(1,1))
	plotGamma(model@param$r,model@param$alpha,model@param$p,spikep=1-sum(model@param$p))
}