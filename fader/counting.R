growth=function(data,len=length(data)){
	return(data[2:len]-data[1:(len-1)])
}

#The poisson model is a discrete one parameter counting model. It answers how many
pois=function(data=data,classes=1,t=1){
	a=list(raw=data,t=t,classes=classes)
	class(a)='pois'
	a$param=model(a)
	
	return(a)
}
ll.pois=function(lambda,num,t){
	return(log(dpois(num,lambda*t)))
}
train.pois=function(param,num,frequency,t=1){
	param=exp(param)
	return(frequency*ll.pois(param,num,t))
}
model.pois=function(model){
	frequency=model$raw
	num=0:(length(model$raw)-1)
	classes=model$classes
	param=optim(runif(2*classes),function(param){
				one=param[(1:classes*2)-1]
				two=param[(1:classes*2)]
				probs=exp(two)/sum(exp(two))
				print(sapply(one,function(x) (train.pois(x,num,frequency,1))))
				return(sum(sapply(1:length(one),function(x) sum(log(probs[x])+train.pois(one[x],num,frequency,model$t)))))
			},control=list(fnscale=-1))$par
	param=data.frame(matrix(exp(param),ncol=2,byrow=TRUE))
	colnames(param)=c('lambda','p')
	param$p=param$p/sum(param$p)
	return(param)
}
mean.pois=function(model,t=1){
	return(sum(model$param$lambda*t*model$param$p))
}
var.pois=mean.pois
predict.pois=function(model,num=0:length(model$raw)-1,t=1){
	n=sum(model$raw)
	param=model$param
	lambda=param$lambda
	p=param$p
	val=apply(sapply(1:length(lambda),function(x){
						p[x]*exp(ll.pois(lambda[x],num,model$t))
					}),1,sum)
	return(n*val)
}
residuals.pois=function(model) predict(model)-model$raw
print.pois=function(x){print(x$param)}
plot.pois=function(model){
	breaks=-0.5:(length(model$raw)-1)
	myhist=list(breaks=breaks,counts=model$raw)
	myhist2=list(breaks=breaks,counts=predict(model))
	class(myhist)='histogram'
	class(myhist2)='histogram'
	plot(myhist,col='black')
	plot(myhist2,add=TRUE,col=rgb(1,0,0,.5))
}
count=read.csv('count.csv')
data.pois=count$num
pois(data.pois)

tmp=runif(1000)<.3
raw=tmp*rpois(1000,2)+(1-tmp)*rpois(1000,50)
pois(sapply(0:max(raw),function(x)length(which(raw==x))),2)
tmpmodel=pois(sapply(0:max(raw),function(x)length(which(raw==x))),2)
plot(tmpmodel)


param=c(log(2),0,log(50),0)
param=c(log(35),0,log(35),0)
one=param[(1:classes*2)-1]
two=param[(1:classes*2)]
probs=exp(two)/sum(exp(two))
print(sum(sapply(1:length(one),function(x) sum(log(probs[x])+train.pois(one[x],num,frequency,1)))))

return(sum(probs*sapply(one,function(x) sum(train.pois(x,num,frequency,model$t)))))

#The poisson model is a discrete two parameter counting model. It answers how many
ll.nbd=function(param,data=data,count=data$count,num=data$num,t=1){
	r=exp(param[1])
	alpha=exp(param[2])
	l=lgamma(r+count)-lgamma(r)-log(factorial(count))+r*(log(alpha/(alpha+t)))+count*(log(t/(alpha+t)))
	return(sum(num*l))
}
model.nbd=function(data,t=1){
	param=exp(optim(runif(2),nbdLL,data=data,t=t,control=list(fnscale=-1))$par)
	return(param)
}
predict.nbd=function(param,data=data,count=data$count,num=data$num,t=1){
	n=sum(num)
	r=param[1]
	alpha=param[2]
	l=lgamma(r+count)-lgamma(r)-log(factorial(count))+r*(log(alpha/(alpha+t)))+count*(log(t/(alpha+t)))
	return(data.frame(count=data$count,num=data$num,model=n*exp(l)))
}
mean.nbd=function(param,t=1){
	r=param[1]
	alpha=param[2]
	return(r*t/alpha)
}
var.nbd=function(param,t=1){
	r=param[1]
	alpha=param[2]
	return(r*t/alpha+r*t^2/alpha^2)
}

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