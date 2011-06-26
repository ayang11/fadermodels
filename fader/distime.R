death=function(data)return(c(data[1:(length(data)-1)]-data[2:length(data)],data[length(data)]))
convert=function(x) rev(cumsum(rev(x)))

#The geometric model is a discrete one parameter timing model. It answers when people will stop ordering. 
control.geom=function(model,...){
	xy=getxy(model$raw,class(model))
	pop=git(list(...)$pop,sum(xy$y))
	return(list(x=c(0,xy$x[-length(xy$x)]),y=c(xy$y,pop-sum(xy$y)),names=c('lambda'),num=pop,...))
}
ll.geom=function(model,param=NULL,x=model$control$x){
	lambda=param[1]
	return(log(c(dgeom(x,lambda),1-pgeom(x[length(x)],lambda))))
}
model.geom=function(model) model.fm(model,zeroone=1)
mean.geom=function(model) return(1/model$param$lambda)
var.geom=function(model) return((1-model$param$lambda)/model$param$lambda^2)


#The Beta Geometric Model is a discrete two parameter timing model. It answers when people will stop ordering. 
control.bg=function(model,...){
	xy=getxy(model$raw,class(model))
	pop=git(list(...)$pop,sum(xy$y))
	return(list(x=c(0,xy$x[-length(xy$x)]),y=c(xy$y,pop-sum(xy$y)),names=c('alpha','beta'),num=pop,...))
}
ll.bg=function(model,param=NULL,x=model$control$x){
	alp=param[1];bet=param[2]
	return(log(c(beta(alp+1,bet+x)/beta(bet,alp),beta(alp,bet+1+x[length(x)])/beta(alp,bet))))
}
mean.bg=function(model) return(model$param$beta/(model$param$alpha-1))
var.bg=function(model){
	alp=model$param$alp;bet=model$param$beta
	return((alp*bet*(alp+bet-1))/((alp-1)^2*(alp-2)))
}


#The discrete weibull is a discrete two parameter timing model. It answers when people will stop ordering.
control.dw=function(model,...){
	xy=getxy(model$raw,class(model))
	pop=git(list(...)$pop,sum(xy$y))
	return(list(x=c(0,xy$x),y=c(xy$y,pop-sum(xy$y)),names=c('lambda','k'),num=pop,...))
}
ll.dw=function(model,param=NULL,x=model$control$x){
	theta=param[1];k=param[2]
	return(log(death((1-theta)^(x^k))))
}
model.dw=function(model) model.fm(model,zeroone=1)


#The Beta discrete weibull is a discrete three parameter timing model. It answers when people will stop ordering. 
control.bdw=function(model,...){
	xy=getxy(model$raw,class(model))
	pop=git(list(...)$pop,sum(xy$y))
	return(list(x=c(0,xy$x),y=c(xy$y,pop-sum(xy$y)),names=c('alpha','beta','k'),num=pop,...))
}
ll.bdw=function(model,param=NULL,x=model$control$x){
	alp=param[1];bet=param[2];k=param[3]
	return(log(death(c(1,beta(alp,bet+(x[-1])^k)/beta(alp,bet)))))
}