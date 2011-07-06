#Standard Functions
sse=function(model)sum(resid(model)^2)
rmse=function(model)return(sqrt(sse(model)/length(model$control$y)))
vcov.fm=function(model) return(NA)
mean.fm=function(model) return(NA)
spike.fm=function(model,param,x=model$control$x) return((x==0)*(1-sum(param$p)))
residuals.fm=function(model,y=model$control$y,...) predict(model,...)-y
print.fm=function(model) {
	cat(toupper(class(model)[1]),'Model\nLL =',model$ll,'\n\nParameters:\n')
	if(hasSpike(model))
		cat('Spike Probability =',1-sum(model$param$p),'\n')
	print(model$param)
}
barplot.fm=function(model,x=model$control$x,act=model$control$y,legend.text=TRUE,mod=predict(model,x=x),...) {
	mat=matrix(c(act,mod),nrow=2,byrow=TRUE)
	rownames(mat)=c('Act','Model')
	colnames(mat)=x[1:ncol(mat)]
	barplot(mat,beside=TRUE,legend.text=legend.text,...)
}
predict.fm=function(model,num=git(model$control$num),...){
	return(num*exp(loglikfunc(model,model$param,...)))
}
update.fm=function(model,...){
	param=data.frame(...)
	for(i in colnames(param))
		if(i %in% colnames(model$param))
			model$param[i]=param[i]
	if(!isTRUE(all.equal(sum(model$param$p),1))){
		if(sum(model$param$p)<1&&sum(model$param$p)>0)
			if(git(model$control$allowspike,FALSE)) model$control$spike=TRUE else stop('Probabilities must sum to 1')
		stop('Probabilities must sum to 1')
	}
	model$ll=sum(git(model$control$y)*loglikfunc(model,model$param))
	model$bic=-2*model$ll+model$numparam*log(git(model$control$num))
	return(model)
}
paramplot.fm=function(model,...){plot(1,1)}
summary.fm=function(model)print.fm(model)
logLik.fm=function(model)return(model$ll)


#reserved variables 
#names is the name of the parameters. This is required.
#x is what will be plotted on the x axis
#y is what will be plotted on the y axis and will be used in multiplication in the ll function
#num is the number of people (used for prediction from probabilities).
#allowspike is whether or not the model allows spikes.
#positives,zerooone,and zeroonesum are used to tell the model that the domain of the parameter at that index is restricted 
#other stuff is used by some models
strip=function(a,...) {
	tmp=list(a,...)
	tmp$names=NULL
	tmp$x=NULL
	tmp$y=NULL
	tmp$num=NULL
	tmp$allowspike=NULL
	tmp$positives=NULL
	tmp$zerooone=NULL
	tmp$zeroonesum=NULL
	tmp$tx=NULL
	tmp$T=NULL
	tmp$plot.y=NULL
	return(tmp)
}
fm=function(data,modname='bb',nseg=1,...){
	mod=list(raw=data,nseg=nseg)
	if(modname %in% c('pnbd','pexp','bgnbd','bgbb'))
		class(mod)=c(modname,'integ','fm')
	else
		class(mod)=c(modname,'fm')
	mod$control=do.call(control,strip(mod,...))
	obj=findMax(mod)
	mod$param=obj$param
	mod$ll=obj$value
	mod$numparam=obj$num
	mod$bic=-2*mod$ll+mod$numparam*log(git(mod$control$num))
	return(mod)
}
findMax=function(model){
	names=c(model$control$names,'p')
	len=length(names)
	obj=findparam(model,len,model$nseg)
	colnames(obj$param)=names
	return(obj)
}
findparam=function(model,dim=1,nseg=1){
	obj=list(value=-Inf)
	spi=hasSpike(model)
	num=dim*nseg-1*!spi
	positives=git(model$control$positives,1:dim)
	zeroone=git(model$control$zeroone,c())
	zeroonesum=git(model$control$zeroonesum,dim)
	for(i in 1:git(model$control$tries)){
		if(i>1) print(paste('Try Number',i,'Best LL =',obj$value))
		errs=1
		repeat{
			if(errs>=5+i) stop('Problem finding a good starting point')
			curr=try(optim(runif(num,-errs*i,errs*i),function(param,model=model){
								param=partition(if(spi) param else c(param,1),cols=dim,zospad=spi,positives=positives,zeroone=zeroone,zeroonesum=zeroonesum)
								colnames(param)=c(model$control$names,'p')
								sum(git(model$control$y)*loglikfunc(model,param))
							},control=list(fnscale=-1),model=model),silent=TRUE)
			if(class(curr)=='try-error') {
				errs=errs+1
				next
			}
			break
		}
		if(obj$value<curr$value)
			obj=curr
	}
	param=obj$par
	value=obj$value
	param=partition(if(spi) param else c(param,1),cols=dim,zospad=spi,positives=positives,zeroone=zeroone,zeroonesum=zeroonesum)
	return(list(param=param,value=value,num=num))
}
partition=function(params,cols=1,positives=1:cols,zeroone=c(),zeroonesum=cols,zospad=0){
	tmp=data.frame(matrix(params,ncol=cols))
	colnames(tmp)=paste('x',1:cols,sep='')
	if(length(union(union(positives,zeroone),zeroonesum))>0) tmp[union(union(positives,zeroone),zeroonesum)]=exp(tmp[union(union(positives,zeroone),zeroonesum)])
	if(length(zeroone)>0) tmp[zeroone]=tmp[zeroone]/(1+tmp[zeroone])
	if(length(zeroonesum)>0) tmp[zeroonesum]=tmp[zeroonesum]/(sum(tmp[zeroonesum])+zospad)
	return(tmp)
}
loglikfunc=function(model,param,...){
	spi=hasSpike(model)
	if(!spi&&nrow(param)==1)
		return(ll(model,param=unlist(param[-length(param)]),...))
	val=(if(spi) spike(model,param) else 0) +rowSums(apply(param,1,weightedlik,model=model,...))
	return(log(val))
}
weightedlik=function(param,model=NULL,...){
	logl=ll(model,param=param[-length(param)],...)
	return(param[length(param)]*exp(logl))
}



#Class Functions
ll=function(model,...) UseMethod('ll')
control=function(model,...) UseMethod('control')
paramplot=function(model,...) UseMethod('paramplot')
condexp=function(model,...) UseMethod('condexp')
palive=function(model,...) UseMethod('palive')
spike=function(model,...) UseMethod('spike')
chitest=function(model) UseMethod('chitest')
indiv=function(model,...) UseMethod('indiv')