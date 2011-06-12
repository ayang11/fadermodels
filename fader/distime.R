death=function(data)return(c(data[1:(length(data)-1)]-data[2:length(data)],data[length(data)]))

#The geometric model is a discrete one parameter timing model. It answers when people will stop ordering. Uses time.csv
control.geom=function(model,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=max(xy$y),mult=death(xy$y),...))
}
ll.geom=function(model,param=NULL,x=model$control$x[-length(model$control$x)]){
	lambda=param[1]
	return(log(c(dgeom(x,lambda),1-pgeom(x[length(x)],lambda))))
}
predict.geom=function(model,...) rev(cumsum(rev(standardpredict(model,...))))
model.geom=function(model) standardmodel(model,c('lambda','p'),zeroone=1)
mean.geom=function(model) return(1/model$param$lambda)
var.geom=function(model) return((1-model$param$lambda)/model$param$lambda^2)
residuals.geom=function(model) standardresid(model)
print.geom=function(model) standardprint(model)
myplot.geom=function(model,...) standardplot(model,...)


#The Beta Geometric Model is a discrete two parameter timing model. It answers when people will stop ordering. 
control.bg=function(model,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=max(xy$y),mult=death(xy$y),...))
}
ll.bg=function(model,param=NULL,x=model$control$x[-length(model$control$x)]){
	alp=param[1];bet=param[2]
	return(log(c(beta(alp+1,bet+x)/beta(bet,alp),beta(alp,bet+1+x[length(x)])/beta(alp,bet))))
}
predict.bg=function(model,...) rev(cumsum(rev(standardpredict(model,...))))
model.bg=function(model) standardmodel(model,c('alpha','beta','p'))
mean.bg=function(model) return(model$param$beta/(model$param$alpha-1))
var.bg=function(model){
	alp=model$param$alp;bet=model$param$beta
	return((alp*bet*(alp+bet-1))/((alp-1)^2*(alp-2)))
}
residuals.bg=function(model) standardresid(model)
print.bg=function(model) standardprint(model)
myplot.bg=function(model,...) standardplot(model,...)


#The discrete weibull is a discrete two parameter timing model. It answers when people will stop ordering.
control.dw=function(model,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=max(xy$y),mult=death(xy$y),...))
}
ll.dw=function(model,param=NULL,x=model$control$x){
	theta=param[1];k=param[2]
	return(log(death((1-theta)^(x^k))))
}
predict.dw=function(model,...) rev(cumsum(rev(standardpredict(model,...))))
model.dw=function(model) standardmodel(model,c('lambda','k','p'),zeroone=1)
mean.dw=function(model) return(NA)
var.dw=function(model) return(NA)
residuals.dw=function(model) standardresid(model)
print.dw=function(model) standardprint(model)
myplot.dw=function(model,...) standardplot(model,...)


#The Beta discrete weibull is a discrete three parameter timing model. It answers when people will stop ordering. 
control.bdw=function(model,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=max(xy$y),mult=death(xy$y),...))
}
ll.bdw=function(model,param=NULL,x=model$control$x){
	alp=param[1];bet=param[2];k=param[3]
	return(log(death(c(1,beta(alp,bet+(x[-1])^k)/beta(alp,bet)))))
}
predict.bdw=function(model,...) rev(cumsum(rev(standardpredict(model,...))))
model.bdw=function(model) standardmodel(model,c('alpha','beta','k','p'))
mean.bdw=function(model) return(NA)
var.bdw=function(model) return(NA)
residuals.bdw=function(model) standardresid(model)
print.bdw=function(model) standardprint(model)
myplot.bdw=function(model,...) standardplot(model,...)