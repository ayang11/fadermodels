#The beta binomial model is a discrete choice model. It answers which. Use count.csv
control.bb=function(model,...){
	xy=getxy(model$raw,class(model))
	return(list(x=xy$x,y=xy$y,len=length(xy$y),num=sum(xy$y),mult=xy$y,...))
}
ll.bb=function(model,param=NULL,x=model$control$x){
	len=max(model$control$x); alp=param[1]; bet=param[2] 
	return(lchoose(len,x)+lgamma(alp+x)+lgamma(bet+len-x)+lgamma(alp+bet)-lgamma(alp+bet+len)-lgamma(alp)-lgamma(bet))
}
predict.bb=function(model,...) standardpredict(model,...)
model.bb=function(model,nseg=1) standardmodel(model,c('alpha','beta','p'),nseg)
mean.bb=function(model,len=max(model$control$x)){
	param=model$param; alps=param$alpha; bets=param$beta
	return(len*alps/(alps+bets))
}
var.bb=function(model,len=max(model$control$x)){
	param=model$param; alps=param$alpha; bets=param$beta; ps=param$p
	return(len*alps*bets*(alps+bets+len)/((alps+bets)^2*(alps+bets+1)))
}
residuals.bb=function(model) standardresid(model)
print.bb=function(model) standardprint(model)
myplot.bb=function(model,...) standardplot(model,...)


#####################Below is not done.

#The dirichlet model is a multinomial discrete choice model. It answers which
control.dir=function(model,...){
	return(list(num=1,mult=1,...))
}
ll.dir=function(model,param=NULL){
	ll=apply(model$raw,1,function(x){
				n=sum(x)
				a=sum(param)
				lfactorial(n)-sum(lfactorial(x))+lgamma(a)+sum(lgamma(param+x))-sum(lgamma(param))-lgamma(a+n)
			})
	return(ll)
}
predict.dir=function(model,...) standardpredict(model,...)
model.dir=function(model,nseg=1) standardmodel(model,c(1:ncol(model$raw),'p'),nseg)
mean.dir=function(model){
	a=dirPredI(model)
	b=dirPenI(model)
	return(lapply(1:length(a),function(i) apply(a[[i]]/b[[i]]/100,2,sum)))
}
var.dir=function(model,n=1){
	params=model$param[-ncol(model$param)]
	nrow=nrow(model$raw)
	ncol=ncol(model$raw)
	res=list()
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		alp1=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
		alp2=matrix(param,nrow=nrow,ncol=ncol,byrow=FALSE)
		s=sum(param)
		res[[i]]=alp1*alp2*(s+n)/(s^2*(s+1))
	}
	return(res)
}
residuals.dir=function(model) standardresid(model)
print.dir=function(model) standardprint(model)
myplot.dir=function(model,...) standardplot(model,...)

#expected number of purchases next period
dirExp=function(model){
	params=model$param[-ncol(model$param)]
	data=model$raw
	nrow=nrow(data)
	ncol=ncol(data)
	res=list()
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		alp=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
		n=matrix(apply(data,1,sum),nrow=nrow,ncol=ncol,byrow=FALSE)
		res[[i]]=(alp+data)/(n+sum(param))
	}
	return(res)
}
#penetration level
dirPen=function(model){
	return(lapply(dirPenI(model),function(x) apply(x,2,mean)))
}
#probability of no purchases
dirNo=function(model){
	params=model$param[-ncol(model$param)]
	data=model$raw
	nrow=nrow(data)
	ncol=ncol(data)
	res=list()
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		s=sum(param)
		alp=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
		n=matrix(apply(data,1,sum),nrow=nrow,ncol=ncol,byrow=FALSE)
		res[[i]]=exp(lgamma(s)+lgamma(s-alp+n)-lgamma(s+n)-lgamma(s-alp))
	}
	return(res)
}
#probability of 100% loyal customer
dirLoyal=function(model){
	a=dirLoyalI(model)
	b=dirPenI(model)
	return(lapply(1:length(a),function(i) apply(a[[i]],2,sum)/apply(b[[i]],2,sum)))
}
#probability of once only customer
dirOnce=function(model){
	a=dirOnceI(model)
	b=dirPenI(model)
	return(lapply(1:length(a),function(i) (apply(a[[i]],2,sum)/apply(b[[i]],2,sum))))
}
#Market share of each segment
dirMkt=function(model){
	params=model$param[-ncol(model$param)]
	res=list()
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		res[[i]]=(param/sum(param))
	}
	return(res)
}
#share of category requirements
dirScr=function(model){
	params=model$param[-ncol(model$param)]
	data=model$raw
	nrow=nrow(data)
	ncol=ncol(data)
	res=list()
	pen=dirPenI(model)
	dmkt=dirMkt(model)
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		mkt=matrix(dmkt[[i]],nrow=nrow,ncol=ncol,byrow=TRUE)
		n=matrix(apply(data,1,sum),nrow=nrow,ncol=ncol,byrow=FALSE)
		res[[i]]=(apply(mkt*n,2,sum)/apply(pen[[i]]*n,2,sum))
	}
	return(res)
}
dirPredI=function(model){
	params=model$param[-ncol(model$param)]
	data=model$raw
	nrow=nrow(data)
	ncol=ncol(data)
	res=list()
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		alp=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
		n=matrix(apply(data,1,sum),nrow=nrow,ncol=ncol,byrow=FALSE)
		res[[i]]=n*alp/sum(param)
	}
	return(res)
}
dirPenI=function(model){
	return(lapply(dirNo(model),function(x) 1-x))
}
dirLoyalI=function(model){
	params=model$param[-ncol(model$param)]
	data=model$raw
	nrow=nrow(data)
	ncol=ncol(data)
	res=list()
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		s=sum(param)
		alp=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
		n=matrix(apply(data,1,sum),nrow=nrow,ncol=ncol,byrow=FALSE)
		res[[i]]=exp(lgamma(alp+n)+lgamma(s)-lgamma(alp)-lgamma(s+n))
	}
	return(res)
}
dirOnceI=function(model){
	params=model$param[-ncol(model$param)]
	data=model$raw
	nrow=nrow(data)
	ncol=ncol(data)
	res=list()
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		s=sum(param)
		alp=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
		n=matrix(apply(data,1,sum),nrow=nrow,ncol=ncol,byrow=FALSE)
		res[[i]]=exp(log(n)+log(alp)+lgamma(s)+lgamma(s-alp+n-1)-lgamma(s-alp)-lgamma(s+n))
	}
	return(res)
}