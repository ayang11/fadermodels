#The beta binomial model is a discrete choice model. It answers which. Use count.csv

#num is number of segments
bb=function(data=data,num=1){
	a=list(param=model.bb(data,num),raw=data)
	class(a)='bb'
	return(a)
}
#helper functions function for the model
ll.bb=function(alp,bet,len){
	x=0:len
	ll=lchoose(len,x)+lgamma(alp+x)+lgamma(bet+len-x)+lgamma(alp+bet)-lgamma(alp+bet+len)-lgamma(alp)-lgamma(bet)
	return(ll)
}
train.bb=function(param,data=data){
	alp=exp(param[1])
	bet=exp(param[2])
	len=length(data)-1
	return(data*ll.bb(alp,bet,len))
}
model.bb=function(data,num){
	param=optim(runif(3*num),function(param,data=data){
						one=param[(1:num)*3-2]
						two=param[(1:num)*3-1]
						three=param[(1:num)*3]
						probs=exp(three)/sum(exp(three))
						return(sum(probs*sapply(1:num,function(x) sum(train.bb(c(one[x],two[x]),data)))))
					},data=data,control=list(fnscale=-1))$par
	param=data.frame(matrix(exp(param),ncol=3,byrow=TRUE))
	colnames(param)=c('alpha','beta','p')
	param$p=param$p/sum(param$p)
	return(param)
}

#functions to use
predict.bb=function(model,len=length(model$raw)-1,num=sum(model$raw)){
	param=model$param
	alps=param$alpha
	bets=param$beta
	ps=param$p
	val=apply(sapply(1:length(alps),function(x) ps[x]*exp(ll.bb(alps[x],bets[x],len))),1,sum)
	return(num*val)
}
mean.bb=function(model,len=length(model$raw)-1){
	param=model$param
	alps=param$alpha
	bets=param$beta
	len*alps/(alps+bets)
}
var.bb=function(model,len=length(model$raw)-1){
	param=model$param
	alps=param$alpha
	bets=param$beta
	ps=param$p
	len*alps*bets*(alps+bets+len)/((alps+bets)^2*(alps+bets+1))
}
residuals.bb=function(model) predict(model)-model$raw
print.bb=function(x){print(x$param)}
plot.bb=function(model){
	breaks=-0.5:length(model$raw)
	myhist=list(breaks=breaks,counts=model$raw)
	myhist2=list(breaks=breaks,counts=predict(model))
	class(myhist)='histogram'
	class(myhist2)='histogram'
	plot(myhist,col='black')
	plot(myhist2,add=TRUE,col=rgb(1,0,0,.5))
}
bb.data=read.csv('bb.csv')$n
model=bb(bb.data)
model2=bb(bb.data,2)
resid(model)
plot(model)


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