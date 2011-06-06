#The poisson model is a discrete one parameter counting model. It answers how many
control.pois=function(model,t=1,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=sum(xy$y),mult=xy$y,t=t,...))
}
ll.pois=function(model,param=NULL,t=model$control$t,x=model$control$x){
	return(log(dpois(x,param*t)))
}
predict.pois=function(model,t=model$control$t,...) standardpredict(model,t=t,...)
model.pois=function(model,nseg=1) standardmodel(model,c('lambda','p'),nseg)
mean.pois=function(model,t=model$control$t) return(sum(model$param$lambda*t*model$param$p))
var.pois=function(model,t=model$control$t) return(model$param$lambda*t)
residuals.pois=function(model,t=model$control$t) predict(model,t=t)-model$control$y
print.pois=function(model) standardprint(model)
myplot.pois=function(model,...) standardplot(model,...)


#The nbd model is a discrete two parameter counting model. It answers how many
control.nbd=function(model,t=1,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=sum(xy$y),mult=xy$y,t=t,...))
}
ll.nbd=function(model,param=NULL,t=model$control$t,x=model$control$x){
	r=param[1];alpha=param[2]
	return(lgamma(r+x)-lgamma(r)-log(factorial(x))+r*(log(alpha/(alpha+t)))+x*(log(t/(alpha+t))))
}
predict.nbd=function(model,t=model$control$t,...) standardpredict(model,t=t,...)
model.nbd=function(model,nseg=1) standardmodel(model,c('r','alpha','p'),nseg)
mean.nbd=function(model,t=model$control$t) return(model$param$r*t/model$param$alpha)
var.nbd=function(model,t=model$control$t) return(model$param$r*t/model$param$alpha+model$param$r*t^2/model$param$alpha^2)
residuals.nbd=function(model,t=model$control$t) predict(model,t=t)-model$control$y
print.nbd=function(model) standardprint(model)
myplot.nbd=function(model,...)standardplot(model,...)


#The exponential gamma model is a continuous two parameter timing model. It answers when
#data vector should be cumulative. Uses kb.csv
growth=function(data,len=length(data)) return(c(data[1],data[2:len]-data[1:(len-1)]))
control.eg=function(model,num=max(xy$y),...){
	xy=getxy(model$raw,class(model))
	y=c(xy$y,num)
	return(list(x=xy$x+1,y=y,len=length(xy$y),num=num,mult=growth(y),...))
}
ll.eg=function(model,param=NULL,x=model$control$x){
	r=param[1];alpha=param[2]
	return(log(c(growth(1-(alpha/(alpha+x))^r),(alpha/(alpha+x[length(x)]))^r)))
}
predict.eg=function(model,...) cumsum(standardpredict(model,...))
model.eg=function(model,nseg=1) standardmodel(model,c('r','alpha','p'),nseg)
mean.eg=function(model) return(model$param$alpha/(model$param$r-1))
var.eg=function(model) return(model$param$r*model$param$alpha^2/((model$param$r-1)^2*(model$param$r-2)))
residuals.eg=function(model) standardresid(model)
print.eg=function(model) standardprint(model)
myplot.eg=function(model,...) standardplot(model,...)


#The weibull gamma model is a continuous three parameter timing model. It answers when
#data vector should be cumulative
control.wg=function(model,num=max(xy$y),...){
	xy=getxy(model$raw,class(model))
	y=c(xy$y,num)
	return(list(x=xy$x+1,y=y,len=length(xy$y),num=num,mult=growth(y),...))
}
ll.wg=function(model,param=NULL,x=model$control$x){
	r=param[1];alpha=param[2];k=param[3]
	return(log(c(growth(1-(alpha/(alpha+x^k))^r),(alpha/(alpha+x[length(x)]^k))^r)))
}
predict.wg=function(model,...) cumsum(standardpredict(model,...))
model.wg=function(model,nseg=1) standardmodel(model,c('r','alpha','c','p'),nseg)
mean.wg=function(model) return(exp(log(model$param$alpha^(1/model$param$c))+lgamma(1+1/model$param$c)+lgamma(model$param$r-1/model$param$c)-lgamma(model$param$r)))
var.wg=function(model) return(model$param$alpha^(2/model$param$c)/gamma(model$param$r)*(gamma(1+2/model$param$c)*gamma(model$param$r-2/model$param$c)-gamma(1+1/model$param$c)^2*gamma(model$param$r-1/model$param$c)^2/gamma(model$param$r)))
residuals.wg=function(model) standardresid(model)
print.wg=function(model) standardprint(model)
myplot.wg=function(model,...) standardplot(model,...)


#The gamma gamma model is a continuous 3 parameter timing model. It answers when. 
#install.packages('gsl')
library(gsl)
control.gg=function(model,num=max(xy$y),...){
	xy=getxy(model$raw,class(model))
	y=c(xy$y,num)
	return(list(x=xy$x+1,y=y,len=length(xy$y),num=num,mult=growth(y),...))
}
ll.gg=function(model,param=NULL,x=model$control$x){
	r=param[1];alpha=param[2];s=param[3]
	tmp=1/(s*(gamma(r)*gamma(s)/gamma(r+s)))*(x/(alpha+x))^s*hyperg_2F1(1-r,s,s+1,x/(alpha+x))
	return(log(c(growth(tmp),1-tmp[length(tmp)])))
}
predict.gg=function(model,...) cumsum(standardpredict(model,...))
model.gg=function(model,nseg=1) standardmodel(model,c('r','alpha','s','p'),nseg)
mean.gg=function(model) return(model$param$alpha*model$param$s/(model$param$r-1))
var.gg=function(model) return(model$param$alpha^2*model$param$s*(model$param$r+model$param$s-1)/((model$param$r-1)^2*(model$param$r-2)))
residuals.gg=function(model) standardresid(model)
print.gg=function(model) standardprint(model)
myplot.gg=function(model,...) standardplot(model,...)