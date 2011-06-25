#The beta binomial model is a discrete choice model. It answers which. Use count.csv
control.bb=function(model,...){
	x=git(model$raw$x)	
	y=git(model$raw$y)	
	return(list(x=x,y=y,m=git(model$raw$m,max(x)),num=sum(y),mult=y,names=c('alpha','beta'),allowspike=TRUE,...))
}
ll.bb=function(model,param=NULL,x=model$control$x){
	m=model$control$m; alp=param[1]; bet=param[2]
	return(lchoose(m,x)+lgamma(alp+x)+lgamma(bet+m-x)+lgamma(alp+bet)-lgamma(alp+bet+m)-lgamma(alp)-lgamma(bet))
}
mean.bb=function(model,len=max(model$control$x)){
	param=model$param; alps=param$alpha; bets=param$beta
	return(len*alps/(alps+bets))
}
var.bb=function(model,len=max(model$control$x)){
	param=model$param; alps=param$alpha; bets=param$beta; ps=param$p
	return(len*alps*bets*(alps+bets+len)/((alps+bets)^2*(alps+bets+1)))
}


#####################Below is not done.
#The dirichlet model is a multinomial discrete choice model. It answers which
control.dir=function(model,...){
	return(list(num=1,mult=1,names=c(1:ncol(model$raw)),...))
}
ll.dir=function(model,param=NULL){
	ll=apply(model$raw,1,function(x){
				n=sum(x)
				a=sum(param)
				lfactorial(n)-sum(lfactorial(x))+lgamma(a)+sum(lgamma(param+x))-sum(lgamma(param))-lgamma(a+n)
			})
	return(ll)
}
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