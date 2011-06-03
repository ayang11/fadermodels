growth=function(data,len=length(data)){
	return(data[2:len]-data[1:(len-1)])
}


control.pois=function(model,t=1,...){
	xy=getxy(model$raw,'Pois')
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=sum(xy$y),mult=xy$y,t=t,...))
}
ll.pois=function(model,param=NULL){
	return(log(dpois(model$control$x,param*model$control$t)))
}
model.pois=function(model,nseg=1){
	param=findparam(model,2,nseg)
	colnames(param)=c('lambda','p')
	return(param)
}

mean.pois=function(model,t=1) return(sum(model$param$lambda*t*model$param$p))
var.pois=function(model,t=1) return(model$param$lambda*t)
predict.pois=function(model,num=model$control$x,t=1){
	param=model$param
	val=apply(apply(param,1,weightedlik,model=model),1,sum)
	return(model$control$num*val)
}
residuals.pois=function(model) predict(model)-model$control$y
print.pois=function(x){print(x$param)}
plot.pois=function(model)fmhistoplot(model,model$control$y)

data.pois=read.csv('count.csv')
mod=fm(data.pois,'pois',6);rmse(mod);mod;plot(mod);mean(mod);var(mod)



#The poisson model is a discrete two parameter counting model. It answers how many
control.nbd=function(model,t=1,...){
	xy=getxy(model$raw,'NBD')
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=sum(xy$y),mult=xy$y,t=t,...))
}
ll.nbd=function(model,param=NULL){
	r=param[1];alpha=param[2];x=model$control$x;t=model$control$t
	return(lgamma(r+x)-lgamma(r)-log(factorial(x))+r*(log(alpha/(alpha+t)))+x*(log(t/(alpha+t))))
}
model.nbd=function(model,nseg=1){
	param=findparam(model,3,nseg)
	colnames(param)=c('r','alpha','p')
	return(param)
}

mean.nbd=function(model,t=1) return(model$param$r*t/model$param$alpha)
var.nbd=function(model,t=1) return(model$param$r*t/model$param$alpha+model$param$r*t^2/model$param$alpha^2)
predict.nbd=function(model,num=model$control$x,t=1){
	param=model$param
	val=apply(apply(param,1,weightedlik,model=model),1,sum)
	return(model$control$num*val)
}
residuals.nbd=function(model) predict(model)-model$control$y
print.nbd=function(x){print(x$param)}
plot.nbd=function(model)fmhistoplot(model,model$control$y)

data.nbd=read.csv('count.csv')
mod=fm(data.nbd,'nbd',1);rmse(mod);mod;plot(mod);mean(mod);var(mod)



#####################I haven't finished anything below
#The exponential gamma model is a continuous two parameter timing model. It answers when
#data vector should be cumulative
ll.eg=function(param,data=data){
	r=exp(param[1])
	alpha=exp(param[2])
	len=length(data)
	new=growth(data)
	F=1-(alpha/(alpha+(1:len)))^r
	f=growth(F,len)
	return(sum(new*log(f)))
}
model.eg=function(data){
	param=exp(optim(runif(2),egLL,data=data,control=list(fnscale=-1))$par)
	return(param)
}
predict.eg=function(param,data=data){
	n=max(data)
	r=param[1]
	alpha=param[2]
	return(data.frame(act=data,model=n*(1-(alpha/(alpha+(1:length(data))))^r)))
}
mean.eg=function(param){
	r=param[1]
	alpha=param[2]
	return(alpha/(r-1))
}	
var.eg=function(param){
	r=param[1]
	alpha=param[2]
	return(r*alpha^2/((r-1)^2*(r-2)))
}

#The weibull gamma model is a continuous three parameter timing model. It answers when
#data vector should be cumulative
ll.wg=function(param,data=data){
	r=exp(param[1])
	alpha=exp(param[2])
	k=exp(param[3])
	len=length(data)
	new=growth(data)
	F=1-(alpha/(alpha+(1:len)^k))^r
	f=growth(F,len)
	return(sum(new*log(f)))
}
model.wg=function(data){
	param=exp(optim(runif(3),wgLL,data=data,control=list(fnscale=-1))$par)
	return(param)
}
predict.wg=function(param,data=data){
	n=max(data)
	r=param[1]
	alpha=param[2]
	k=param[3]
	return(data.frame(act=data,model=n*(1-(alpha/(alpha+(1:length(data))^k))^r)))
}
mean.wg=function(param){
	r=param[1]
	alpha=param[2]
	k=param[3]
	return(exp(log(alpha^(1/k))+lgamma(1+1/k)+lgamma(r-1/k)-lgamma(r)))
}	
var.wg=function(param){
	r=param[1]
	alpha=param[2]
	k=param[3]
	return(alpha^(2/k)/gamma(r)*(gamma(1+2/c)*gamma(r-2/c)-gamma(1+1/k)^2*gamma(r-1/k)^2/gamma(r)))
}



ll.gg=function(param,data=data){
	r=exp(param[1])
	alpha=exp(param[2])
	s=exp(param[3])
	len=length(data)
	t=1:len
	F=1/(s*(gamma(r)*gamma(s)/gamma(r+s)))*(t/(alpha+t))^s*hyperg_2F1(1-r,s,s+1,t/(alpha+t))
	new=growth(data)
	f=growth(F,len)
	return(sum(new*log(f)))
}
model.gg=function(data){
	param=exp(optim(runif(3),ggLL,data=data,control=list(fnscale=-1))$par)
	return(param)
}
predict.gg=function(param,data=data){
	n=max(data)
	r=param[1]
	alpha=param[2]
	k=param[3]
	t=1:length(data)
	return(data.frame(act=data,model=n*(1/(s*(gamma(r)*gamma(s)/gamma(r+s)))*(t/(alpha+t))^s*hyperg_2F1(1-r,s,s+1,t/(alpha+t)))))
}
mean.gg=function(param){
	r=param[1]
	alpha=param[2]
	s=param[3]
	return(alpha*s/(r-1))
}	
var.gg=function(param){
	r=param[1]
	alpha=param[2]
	s=param[3]
	return(alpha^2*s*(r+s-1)/((r-1)^2*(r-2)))
}