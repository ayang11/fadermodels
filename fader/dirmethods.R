#penetration level. Average penetration for each consumer
pen.dir=function(model){
	return(lapply(pen.ind.dir(model),function(x) colMeans(x)))
}
#probability that a customer will be 100% loyal in your segment
loyal.dir=function(model){
	a=loyal.ind.dir(model)
	b=pen.ind.dir(model)
	return(lapply(1:length(a),function(i) colSums(a[[i]])/colSums(b[[i]])))
}
#probability that a customer will be a one time customer in your segment
once.dir=function(model){
	a=once.ind.dir(model)
	b=pen.ind.dir(model)
	return(lapply(1:length(a),function(i) (colSums(a[[i]])/colSums(b[[i]]))))
}
#Market share of each segment
mkt.dir=function(model){
	params=model@param[-ncol(model@param)]
	res=list()
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		res[[i]]=(param/sum(param))
	}
	return(res)
}
#share of category requirements. share of brands that are bought at least once
scr.dir=function(model){
	params=model@param[-ncol(model@param)]
	data=model@raw
	nrow=nrow(data)
	ncol=ncol(data)
	res=list()
	pen=pen.ind.dir(model)
	dmkt=mkt.dir(model)
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		mkt=matrix(dmkt[[i]],nrow=nrow,ncol=ncol,byrow=TRUE)
		n=matrix(rowSums(data),nrow=nrow,ncol=ncol,byrow=FALSE)
		res[[i]]=(colSums(mkt*n)/colSums(pen[[i]]*n))
	}
	return(res)
}
#individual level probability of no purchases. It is P(0 purchases in a segment)
no.ind.dir=function(model){
	params=model@param[-ncol(model@param)]
	data=model@raw
	nrow=nrow(data)
	ncol=ncol(data)
	res=list()
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		s=sum(param)
		alp=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
		n=matrix(rowSums(data),nrow=nrow,ncol=ncol,byrow=FALSE)
		res[[i]]=exp(lgamma(s)+lgamma(s-alp+n)-lgamma(s+n)-lgamma(s-alp))
	}
	return(res)
}
#individual level expected number of purchases this period. Expected frequency per buyer
freq.ind.dir=function(model){
	params=model@param[-ncol(model@param)]
	data=model@raw
	nrow=nrow(data)
	ncol=ncol(data)
	res=list()
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		alp=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
		n=matrix(rowSums(data),nrow=nrow,ncol=ncol,byrow=FALSE)
		res[[i]]=n*alp/sum(param)
	}
	return(res)
}
#individual level penetration probability. It is 1-P(zero purchases)
pen.ind.dir=function(model){
	return(lapply(no.ind.dir(model),function(x) 1-x))
}
#individual level 100% loyal probability. It is P(n purchases from that segment)
loyal.ind.dir=function(model){
	params=model@param[-ncol(model@param)]
	data=model@raw
	nrow=nrow(data)
	ncol=ncol(data)
	res=list()
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		s=sum(param)
		alp=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
		n=matrix(rowSums(data),nrow=nrow,ncol=ncol,byrow=FALSE)
		res[[i]]=exp(lgamma(alp+n)+lgamma(s)-lgamma(alp)-lgamma(s+n))
	}
	return(res)
}
#individual level 100% loyal probability. It is P(1 purchase in segment)
once.ind.dir=function(model){
	params=model@param[-ncol(model@param)]
	data=model@raw
	nrow=nrow(data)
	ncol=ncol(data)
	res=list()
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		s=sum(param)
		alp=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
		n=matrix(rowSums(data),nrow=nrow,ncol=ncol,byrow=FALSE)
		res[[i]]=exp(log(n)+log(alp)+lgamma(s)+lgamma(s-alp+n-1)-lgamma(s-alp)-lgamma(s+n))
	}
	return(res)
}
#expected number of purchases next period if given one more transaction opportunity
exp.ind.dir=function(model){
	params=model@param[-ncol(model@param)]
	data=model@raw
	nrow=nrow(data)
	ncol=ncol(data)
	res=list()
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		alp=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
		n=matrix(rowSums(data),nrow=nrow,ncol=ncol,byrow=FALSE)
		res[[i]]=(alp+data)/(n+sum(param))
	}
	return(res)
}