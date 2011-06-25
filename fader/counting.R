#The poisson model is a discrete one parameter counting model. It answers how many
control.pois=function(model,t=1,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,num=sum(xy$y),mult=xy$y,t=t,names=c('lambda','p'),allowspike=TRUE,...))
}
ll.pois=function(model,param=NULL,t=model$control$t,x=model$control$x){
	return(log(dpois(x,param*t)))
}
mean.pois=function(model,t=model$control$t) return(sum(model$param$lambda*t*model$param$p))
var.pois=function(model,t=model$control$t) return(model$param$lambda*t)

#The nbd model is a discrete two parameter counting model. It answers how many
control.nbd=function(model,t=1,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,num=sum(xy$y),mult=xy$y,t=t,names=c('r','alpha','p'),allowspike=TRUE,...))
}
ll.nbd=function(model,param=NULL,t=model$control$t,x=model$control$x){
	r=param[1];alpha=param[2]
	return(lgamma(r+x)-lgamma(r)-log(factorial(x))+r*(log(alpha/(alpha+t)))+x*(log(t/(alpha+t))))
}
mean.nbd=function(model,t=model$control$t) return(model$param$r*t/model$param$alpha)
var.nbd=function(model,t=model$control$t) return(model$param$r*t/model$param$alpha+model$param$r*t^2/model$param$alpha^2)