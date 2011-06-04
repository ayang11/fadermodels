#The beta binomial model is a discrete choice model. It answers which. Use count.csv
control.bb=function(model,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=sum(xy$y),mult=xy$y,...))
}
#helper functions function for the model
ll.bb=function(model,param=NULL){
	len=max(model$control$x); alp=param[1]; bet=param[2]; x=model$control$x
	return(lchoose(len,x)+lgamma(alp+x)+lgamma(bet+len-x)+lgamma(alp+bet)-lgamma(alp+bet+len)-lgamma(alp)-lgamma(bet))
}
model.bb=function(model,nseg=1){
	param=findparam(model,3,nseg)
	colnames(param)=c('alpha','beta','p')
	return(param)
}

predict.bb=function(model){
	param=model$param
	val=apply(apply(param,1,weightedlik,model=model),1,sum)
	return(model$control$num*val)
}
mean.bb=function(model,len=max(model$control$x)){
	param=model$param; alps=param$alpha; bets=param$beta
	return(len*alps/(alps+bets))
}
var.bb=function(model,len=max(model$control$x)){
	param=model$param; alps=param$alpha; bets=param$beta; ps=param$p
	return(len*alps*bets*(alps+bets+len)/((alps+bets)^2*(alps+bets+1)))
}
residuals.bb=function(model) predict(model)-model$control$y
print.bb=function(x){print(x$param)}
plot.bb=function(model)	fmhistoplot(model,model$control$y)


data.bb=read.csv('bb.csv')
mod=fm(data.bb,nseg=2);rmse(mod);mod;plot(mod);mean(mod);var(mod)



#####################Below is not done.


#The dirichlet model is a multinomial discrete choice model. It answers which
ll.dir=function(param,data=data){
	alp=exp(param)
	ll=apply(data,1,function(x){
				n=sum(x)
				a=sum(alp)
				lfactorial(n)-sum(lfactorial(x))+lgamma(a)+sum(lgamma(alp+x))-sum(lgamma(alp))-lgamma(a+n)
			})
	return(sum(ll))
}
model.dir=function(data,tries=10){
	attempts=sapply(1:tries,function(x){
				tmp=optim(runif(ncol(data),-5,5),dirLL,data=data,control=list(fnscale=-1,reltol=10^-11))
				return(c(tmp$par,tmp$value))
			})
	print(attempts)
	n=nrow(attempts)
	loc=which.max(attempts[n,])
	param=exp(attempts[-n,loc])
	return(param)
}
#expected number of purchases next period
dirExp=function(param,data=data){
	nrow=nrow(data)
	ncol=ncol(data)
	alp=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
	n=matrix(apply(data,1,sum),nrow=nrow,ncol=ncol,byrow=FALSE)
	exp=(alp+data)/(n+sum(param))
	return(exp)
}
#expected number of purchases this period
dirFreq=function(param,data=data){
	return(apply(dirPredI(param,data)/(dirPenI(param,data))/100,2,sum))
}
dirPredI=function(param,data=data){
	nrow=nrow(data)
	ncol=ncol(data)
	alp=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
	n=matrix(apply(data,1,sum),nrow=nrow,ncol=ncol,byrow=FALSE)
	pred=n*alp/sum(param)
	return(pred)
}
#penetration level
dirPen=function(param,data=data){
	return(apply(dirPenI(param,data),2,mean))
}
dirPenI=function(param,data=data){
	return(1-dirNo(param,data))
}
#probability of no purchases
dirNo=function(param,data=data){
	s=sum(param)
	nrow=nrow(data)
	ncol=ncol(data)
	alp=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
	n=matrix(apply(data,1,sum),nrow=nrow,ncol=ncol,byrow=FALSE)
	p0=lgamma(s)+lgamma(s-alp+n)-lgamma(s+n)-lgamma(s-alp)
	return(exp(p0))
}
#probability of 100% loyal customer
dirLoyal=function(param,data=data){
	return(apply(dirLoyalI(param,data),2,sum)/apply(dirPenI(param,data),2,sum))
}
dirLoyalI=function(param,data=data){
	s=sum(param)
	nrow=nrow(data)
	ncol=ncol(data)
	alp=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
	n=matrix(apply(data,1,sum),nrow=nrow,ncol=ncol,byrow=FALSE)
	loyal=lgamma(alp+n)+lgamma(s)-lgamma(alp)-lgamma(s+n)
	return(exp(loyal))
}
#probability of once only customer
dirOnce=function(param,data=data){
	return(apply(dirOnceI(param,data),2,sum)/apply(dirPenI(param,data),2,sum))
}
dirOnceI=function(param,data=data){
	s=sum(param)
	nrow=nrow(data)
	ncol=ncol(data)
	alp=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
	n=matrix(apply(data,1,sum),nrow=nrow,ncol=ncol,byrow=FALSE)
	once=log(n)+log(alp)+lgamma(s)+lgamma(s-alp+n-1)-lgamma(s-alp)-lgamma(s+n)
	return(exp(once))
}
#Market share of each segment
dirMkt=function(param,data=data){
	return(param/sum(param))
}
#share of category requirements
dirScr=function(param,data=data){
	nrow=nrow(data)
	ncol=ncol(data)
	mkt=matrix(dirMkt(param,data),nrow=nrow,ncol=ncol,byrow=TRUE)
	pen=dirPenI(param,data)
	n=matrix(apply(data,1,sum),nrow=nrow,ncol=ncol,byrow=FALSE)
	return(apply(mkt*n,2,sum)/apply(pen*n,2,sum))
}
var.dir=function(param,n=1){
	alp1=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
	alp2=matrix(param,nrow=nrow,ncol=ncol,byrow=FALSE)
	s=sum(param)
	return(alp1*alp2*(s+n)/(s^2*(s+1)))
}
#mod=dirModel(data)
mod=c(0.195,0.054,0.060,0.116,0.362,0.139,0.151,0.175)

dirPen(mod,data)
dirFreq(mod,data)
dirOnce(mod,data)
dirLoyal(mod,data)
dirMkt(mod,data)
dirScr(mod,data)