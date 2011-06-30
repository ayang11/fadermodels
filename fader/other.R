#The was an experimental normal model
setClass('norm',contains='fm')
control.norm=function(model,...){
	return(list(x=model@raw,y=model@raw,names=c('mean','sd'),positives=2,...))
}
ll.norm=function(model,param=NULL,x=model@control$x){
	return(log(dnorm(x,param[1],param[2])))
}
mean.norm=function(model) return(sum(model@param$mean*model@param$p))
vcov.norm=function(model) return(model@param$lambda)