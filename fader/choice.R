#The beta binomial model is a discrete choice model. It answers which. Use count.csv
setClass('bb',contains='fm')
control.bb=function(model,...){
	checkdata(model,condition=c('x','y'))
	return(list(x=model@raw$x,y=model@raw$y,m=git(model@raw$m,max(model@raw$x)),num=sum(model@raw$y),names=c('alpha','beta'),allowspike=TRUE,...))
}
ll.bb=function(model,param=NULL,x=model@control$x){
	m=model@control$m; alp=param[1]; bet=param[2]
	return(lchoose(m,x)+lgamma(alp+x)+lgamma(bet+m-x)+lgamma(alp+bet)-lgamma(alp+bet+m)-lgamma(alp)-lgamma(bet))
}
mean.bb=function(model,len=max(model@control$x)){
	param=model@param; alps=param$alpha; bets=param$beta
	warning("This hasn't been checked")
	return(len*alps/(alps+bets))
}
vcov.bb=function(model,len=max(model@control$x)){
	param=model@param; alps=param$alpha; bets=param$beta; ps=param$p
	warning("This hasn't been checked")
	return(len*alps*bets*(alps+bets+len)/((alps+bets)^2*(alps+bets+1)))
}
paramplot.bb=function(model,...) {
	par(mfrow=c(1,1))
	plotBeta(model@param$alpha,model@param$beta,model@param$p,spikep=1-sum(model@param$p),...)
}


#The dirichlet model is a multinomial discrete choice model. It answers which
setClass('dir',contains='fm')
control.dir=function(model,...) {
	return(list(x=colnames(model@raw),y=rep(1,nrow(model@raw)),num=nrow(model@raw),plot.y=colSums(model@raw),names=paste('x',c(1:ncol(model@raw)),sep=''),...))
}
ll.dir=function(model,param=NULL){
	ll=apply(model@raw,1,function(x){
				n=sum(x)
				a=sum(param)
				lfactorial(n)-sum(lfactorial(x))+lgamma(a)+sum(lgamma(param+x))-sum(lgamma(param))-lgamma(a+n)
			})
	return(ll)
}
predict.dir=function(model,...) colSums(model@param$p*t(sapply(freq.ind.dir(model),function(x) colSums(x))))
mean.dir=function(model){
	a=freq.ind.dir(model)
	b=pen.ind.dir(model)
	warning("This hasn't been checked")
	return(lapply(1:length(a),function(i) colSums(a[[i]]/b[[i]]/100)))
}
vcov.dir=function(model,n=1){
	params=model@param[-ncol(model@param)]
	nrow=nrow(model@raw)
	ncol=ncol(model@raw)
	res=list()
	for(i in 1:nrow(params)){
		param=unlist(params[i,])
		alp1=matrix(param,nrow=nrow,ncol=ncol,byrow=TRUE)
		alp2=matrix(param,nrow=nrow,ncol=ncol,byrow=FALSE)
		s=sum(param)
		res[[i]]=alp1*alp2*(s+n)/(s^2*(s+1))
	}
	warning("This hasn't been checked")
	return(res)
}
chitest.dir=function(model,...) chitest.fm(model,act=model@control$plot.y)
barplot.dir=function(model,...) barplot.fm(model,x=model@control$x,act=model@control$plot.y,...)
residuals.dir=function(model,...) residuals.fm(model,y=model@control$plot.y,...)