death=function(data)return(c(data[1:(length(data)-1)]-data[2:length(data)],data[length(data)]))

#dProb is death probability
#sProb is survival probability
#lProb is living probability

#The geometric model is a discrete one parameter timing model. It answers when people will stop ordering. Uses time.csv
control.geom=function(model,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=max(xy$y),mult=death(xy$y),...))
}
ll.geom=function(model,param=NULL){
	lambda=param[1];x=model$control$x[-length(model$control$x)];
	return(log(c(dgeom(x,lambda),1-pgeom(x[length(x)],lambda))))
}
model.geom=function(model,nseg=1){
	param=findparam(model,2,nseg,func=specialpartition2)
	colnames(param)=c('lambda','p')
	return(param)
}

mean.geom=function(model) return(1/model$param$lambda)
var.geom=function(model) return((1-model$param$lambda)/model$param$lambda^2)
predict.geom=function(model){
	param=model$param
	val=apply(apply(param,1,weightedlik,model=model),1,sum)
	return(rev(cumsum(rev(model$control$num*val))))
}
residuals.geom=function(model) predict(model)-model$control$y
print.geom=function(x){print(x$param)}
plot.geom=function(model)fmhistoplot(model,model$control$y)

data.geom=read.csv('time.csv')
mod=fm(data.geom,'geom',1);rmse(mod);mod;plot(mod);mean(mod);var(mod)



#The Beta Geometric Model is a discrete two parameter timing model. It answers when people will stop ordering. 
control.bg=function(model,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=max(xy$y),mult=death(xy$y),...))
}
ll.bg=function(model,param=NULL){
	alp=param[1];bet=param[2];x=model$control$x[-length(model$control$x)]
	return(log(c(beta(alp+1,bet+x)/beta(bet,alp),beta(alp,bet+1+x[length(x)])/beta(alp,bet))))
}
model.bg=function(model,nseg=1){
	param=findparam(model,3,nseg)
	colnames(param)=c('alpha','beta','p')
	return(param)
}

mean.bg=function(model) return(1/model$param$lambda)
var.bg=function(model) return((1-model$param$lambda)/model$param$lambda^2)
predict.bg=function(model){
	param=model$param
	val=apply(apply(param,1,weightedlik,model=model),1,sum)
	return(rev(cumsum(rev(model$control$num*val))))
}
residuals.bg=function(model) predict(model)-model$control$y
print.bg=function(x){print(x$param)}
plot.bg=function(model)fmhistoplot(model,model$control$y)

data.bg=read.csv('time.csv')
mod=fm(data.bg,'bg',1);rmse(mod);mod;plot(mod);mean(mod);var(mod)



#The discrete weibull is a discrete two parameter timing model. It answers when people will stop ordering.
control.dw=function(model,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=max(xy$y),mult=death(xy$y),...))
}
ll.dw=function(model,param=NULL){
	theta=param[1];k=param[2];x=model$control$x
	return(log(death((1-theta)^(x^k))))
}
model.dw=function(model,nseg=1){
	param=findparam(model,3,nseg,func=specialpartition3)
	colnames(param)=c('lambda','k','p')
	return(param)
}

mean.dw=function(model) return(NA)
var.dw=function(model) return(NA)
predict.dw=function(model){
	param=model$param
	val=apply(apply(param,1,weightedlik,model=model),1,sum)
	return(rev(cumsum(rev(model$control$num*val))))
}
residuals.dw=function(model) predict(model)-model$control$y
print.dw=function(x){print(x$param)}
plot.dw=function(model)fmhistoplot(model,model$control$y)

data.dw=read.csv('time.csv')
mod=fm(data.dw,'dw',2);rmse(mod);mod;plot(mod);mean(mod);var(mod)

#The Beta discrete weibull is a discrete three parameter timing model. It answers when people will stop ordering. 
control.bdw=function(model,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=max(xy$y),mult=death(xy$y),...))
}
ll.bdw=function(model,param=NULL){
	alp=param[1];bet=param[2];k=param[3];x=model$control$x
	return(log(death(c(1,beta(alp,bet+(x[-1])^k)/beta(alp,bet)))))
}
model.bdw=function(model,nseg=1){
	param=findparam(model,4,nseg)
	colnames(param)=c('alpha','beta','k','p')
	return(param)
}

mean.bdw=function(model) return(NA)
var.bdw=function(model) return(NA)
predict.bdw=function(model){
	param=model$param
	val=apply(apply(param,1,weightedlik,model=model),1,sum)
	return(rev(cumsum(rev(model$control$num*val))))
}
residuals.bdw=function(model) predict(model)-model$control$y
print.bdw=function(x){print(x$param)}
plot.bdw=function(model)fmhistoplot(model,model$control$y)

data.bdw=read.csv('time.csv')
mod=fm(data.bdw,'bdw',1);rmse(mod);mod;plot(mod);mean(mod);var(mod)