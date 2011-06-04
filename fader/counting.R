control.pois=function(model,t=1,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=sum(xy$y),mult=xy$y,t=t,...))
}
ll.pois=function(model,param=NULL,t=model$control$t){
	return(log(dpois(model$control$x,param*t)))
}
model.pois=function(model,nseg=1){
	param=findparam(model,2,nseg)
	colnames(param)=c('lambda','p')
	return(param)
}

mean.pois=function(model,t=model$control$t) return(sum(model$param$lambda*t*model$param$p))
var.pois=function(model,t=model$control$t) return(model$param$lambda*t)
predict.pois=function(model,num=model$control$x,t=model$control$t){
	param=model$param
	val=apply(apply(param,1,weightedlik,model=model,t=t),1,sum)
	return(model$control$num*val)
}
residuals.pois=function(model,t=model$control$t) predict(model,t=t)-model$control$y
print.pois=function(x){print(x$param)}
plot.pois=function(model,t=model$control$t)fmhistoplot(model,model$control$y,t=t)

data.pois=read.csv('count.csv')
mod=fm(data.pois,'pois',6,t=3);rmse(mod);mod;plot(mod,t=1);mean(mod);var(mod)


#The poisson model is a discrete two parameter counting model. It answers how many
control.nbd=function(model,t=1,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=sum(xy$y),mult=xy$y,t=t,...))
}
ll.nbd=function(model,param=NULL,t=model$control$t){
	r=param[1];alpha=param[2];x=model$control$x
	return(lgamma(r+x)-lgamma(r)-log(factorial(x))+r*(log(alpha/(alpha+t)))+x*(log(t/(alpha+t))))
}
model.nbd=function(model,nseg=1){
	param=findparam(model,3,nseg)
	colnames(param)=c('r','alpha','p')
	return(param)
}

mean.nbd=function(model,t=model$control$t) return(model$param$r*t/model$param$alpha)
var.nbd=function(model,t=model$control$t) return(model$param$r*t/model$param$alpha+model$param$r*t^2/model$param$alpha^2)
predict.nbd=function(model,t=model$control$t){
	param=model$param
	val=apply(apply(param,1,weightedlik,model=model,t=t),1,sum)
	return(model$control$num*val)
}
residuals.nbd=function(model,t=model$control$t) predict(model,t=t)-model$control$y
print.nbd=function(x){print(x$param)}
plot.nbd=function(model,t=model$control$t)fmhistoplot(model,model$control$y,t=t)

data.nbd=read.csv('count.csv')
mod=fm(data.nbd,'nbd',2);rmse(mod);mod;plot(mod,t=0.5);mean(mod);var(mod)




#The exponential gamma model is a continuous two parameter timing model. It answers when
#data vector should be cumulative. Uses kb.csv
growth=function(data,len=length(data)) return(c(data[1],data[2:len]-data[1:(len-1)]))
control.eg=function(model,num=max(xy$y),...){
	xy=getxy(model$raw,class(model))
	y=c(xy$y,num)
	return(list(x=xy$x+1,y=y,len=length(xy$y),num=num,mult=growth(y),...))
}
ll.eg=function(model,param=NULL){
	r=param[1];alpha=param[2];x=model$control$x;
	return(log(c(growth(1-(alpha/(alpha+x))^r),(alpha/(alpha+x[length(x)]))^r)))
}
model.eg=function(model,nseg=1){
	param=findparam(model,3,nseg)
	colnames(param)=c('r','alpha','p')
	return(param)
}

mean.eg=function(model) return(model$param$alpha/(model$param$r-1))
var.eg=function(model) return(model$param$r*model$param$alpha^2/((model$param$r-1)^2*(model$param$r-2)))
predict.eg=function(model){
	param=model$param
	val=apply(apply(param,1,weightedlik,model=model),1,sum)
	return(cumsum(model$control$num*val))
}
residuals.eg=function(model) predict(model)-model$control$y
print.eg=function(x){print(x$param)}
plot.eg=function(model)fmhistoplot(model,model$control$y)

data.eg=read.csv('kb.csv',header=FALSE)
mod=fm(data.eg,'eg',1,num=1499);rmse(mod);mod;plot(mod);mean(mod);var(mod)


#The weibull gamma model is a continuous three parameter timing model. It answers when
#data vector should be cumulative
control.wg=function(model,num=max(xy$y),...){
	xy=getxy(model$raw,class(model))
	y=c(xy$y,num)
	return(list(x=xy$x+1,y=y,len=length(xy$y),num=num,mult=growth(y),...))
}
ll.wg=function(model,param=NULL){
	r=param[1];alpha=param[2];k=param[3];x=model$control$x;
	return(log(c(growth(1-(alpha/(alpha+x^k))^r),(alpha/(alpha+x[length(x)]^k))^r)))
}
model.wg=function(model,nseg=1){
	param=findparam(model,4,nseg)
	colnames(param)=c('r','alpha','c','p')
	return(param)
}

mean.wg=function(model) return(exp(log(model$param$alpha^(1/model$param$c))+lgamma(1+1/model$param$c)+lgamma(model$param$r-1/model$param$c)-lgamma(model$param$r)))
var.wg=function(model) return(model$param$alpha^(2/model$param$c)/gamma(model$param$r)*(gamma(1+2/model$param$c)*gamma(model$param$r-2/model$param$c)-gamma(1+1/model$param$c)^2*gamma(model$param$r-1/model$param$c)^2/gamma(model$param$r)))
predict.wg=function(model){
	param=model$param
	val=apply(apply(param,1,weightedlik,model=model),1,sum)
	return(cumsum(model$control$num*val))
}
residuals.wg=function(model) predict(model)-model$control$y
print.wg=function(x){print(x$param)}
plot.wg=function(model)fmhistoplot(model,model$control$y)

data.wg=read.csv('kb.csv',header=FALSE)
mod=fm(data.wg,'wg',1,num=1499);rmse(mod);mod;plot(mod);mean(mod);var(mod)


#The gamma gamma model is a continuous 3 parameter timing model. It answers when. 
#install.packages('gsl')
library(gsl)
control.gg=function(model,num=max(xy$y),...){
	xy=getxy(model$raw,class(model))
	y=c(xy$y,num)
	return(list(x=xy$x+1,y=y,len=length(xy$y),num=num,mult=growth(y),...))
}
ll.gg=function(model,param=NULL){
	r=param[1];alpha=param[2];s=param[3];t=model$control$x;
	tmp=1/(s*(gamma(r)*gamma(s)/gamma(r+s)))*(t/(alpha+t))^s*hyperg_2F1(1-r,s,s+1,t/(alpha+t))
	return(log(c(growth(tmp),1-tmp[length(tmp)])))
}
model.gg=function(model,nseg=1){
	param=findparam(model,4,nseg)
	colnames(param)=c('r','alpha','s','p')
	return(param)
}

mean.gg=function(model) return(model$param$alpha*model$param$s/(model$param$r-1))
var.gg=function(model) return(model$param$alpha^2*model$param$s*(model$param$r+model$param$s-1)/((model$param$r-1)^2*(model$param$r-2)))
predict.gg=function(model){
	param=model$param
	val=apply(apply(param,1,weightedlik,model=model),1,sum)
	return(cumsum(model$control$num*val))
}
residuals.gg=function(model) predict(model)-model$control$y
print.gg=function(x){print(x$param)}
plot.gg=function(model)fmhistoplot(model,model$control$y)

data.gg=read.csv('kb.csv',header=FALSE)
mod=fm(data.wg,'gg',1,num=1499);rmse(mod);mod;plot(mod);mean(mod);var(mod)