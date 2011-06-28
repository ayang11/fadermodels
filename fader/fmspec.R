#options(warn=-1)

#Standard Functions
vcov.fm=function(model) return(NA)
mean.fm=function(model) return(NA)
spike.fm=function(model,param,x=model@control$x) return((x==0)*(1-sum(param$p)))
residuals.fm=function(model,...) predict(model,...)-model@control$y
print.fm=function(model) {
	print(paste(toupper(class(model)[1]),'Model LL =',model@ll))
	if(git(model@control$allowspike,FALSE)&&git(model@control$spike,FALSE))
		print(paste('Spike P =',1-sum(model@param$p)))
	print(model@param)
}
barplot.fm=function(model,x=model@control$x,mod=model@control$y,legend.text=TRUE,act=predict(model,x=x),...) {
	mat=matrix(c(mod,act),nrow=2,byrow=TRUE)
	rownames(mat)=c('Act','Model')
	colnames(mat)=x[1:ncol(mat)]
	barplot(mat,beside=TRUE,legend.text=legend.text,...)
}
model.fm=function(model,...){
	names=c(model@control$names,'p')
	len=length(names)
	obj=findparam(model,len,model@nseg,...)
	colnames(obj$param)=names
	return(obj)
}
likfunc=function(model,param,...){
	spi=hasSpike(model)
	val=(if(spi) spike(model,param) else 0) +apply(apply(param,1,weightedlik,model=model,...),1,sum)
	return(val)
}
predict.fm=function(model,num=git(model@control$num),...){
	return(num*likfunc(model,model@param,...))
}
update.fm=function(model,...){
	param=data.frame(...)
	for(n in colnames(model@param))
		if(!(n %in% colnames(param)))
			param=cbind(param,model@param[n])
	if(!isTRUE(all.equal(sum(param$p),1)))
		if(sum(param$p)<1&&sum(param$p)>0)
			if(git(model@control$allowspike,FALSE)) model@control$spike=TRUE else stop('Probabilities must sum to 1')
	if(nrow(param)==nrow(model@param)&&ncol(param)==ncol(model@param)&&colnames(param)==colnames(model@param)){
		model@param=param
		model@ll=sum(git(model@control$y)*log(likfunc(model,model@param)))
		model@bic=-2*model@ll+model@numparam*log(git(model@control$num))
		return(model)
	}
	stop("Please keep number of segments constant and don't add extra parameters")
}
post=function(model,...){
	lik=apply(model@param,1,weightedlik,model=model,...)
	colnames(lik)=c(1:ncol(lik))
	if(hasSpike(model)){
		lik=cbind(spike(model,model@param),lik)
		colnames(lik)=c('Spike',1:(ncol(lik)-1))
	}
	total=apply(lik,1,sum)
	return(lik/total)
}

#reserved variables 
#names is the name of the parameters. This is required.
#x is what will be plotted on the x axis
#y is what will be plotted on the y axis
#num is the number of people (used for prediction from probabilities).
#allowspike is whether or not the model allows spikes.
strip=function(a,names=NULL,x=NULL,y=NULL,num=NULL,allowspike=NULL,...) return(list(a,...))
fm=function(data,modname='bb',nseg=1,...){
	mod=new(modname,raw=data,nseg=nseg)
	mod@control=do.call(control,strip(mod,...))
	obj=model(mod)
	mod@param=obj$param
	mod@ll=obj$value
	mod@numparam=obj$num
	mod@bic=-2*mod@ll+mod@numparam*log(git(mod@control$num))
	return(mod)
}
partition=function(params,by=1,pos=1:by,zeroone=c(),zeroonesum=by,zospad=0){
	tmp=data.frame(matrix(params,ncol=by))
	colnames(tmp)=paste('x',1:by,sep='')
	if(length(union(union(pos,zeroone),zeroonesum))>0) tmp[union(union(pos,zeroone),zeroonesum)]=exp(tmp[union(union(pos,zeroone),zeroonesum)])
	if(length(zeroone)>0) tmp[zeroone]=tmp[zeroone]/(1+tmp[zeroone])
	if(length(zeroonesum)>0) tmp[zeroonesum]=tmp[zeroonesum]/(sum(tmp[zeroonesum])+zospad)
	return(tmp)
}
weightedlik=function(param,model=NULL,...){
	logl=ll(model,param=param[-length(param)],...)
	return(param[length(param)]*exp(logl))
}
findparam=function(model,dim=1,nseg=1,...){
	obj=list(value=-Inf)
	spi=hasSpike(model)
	num=dim*nseg-1*!spi
	for(i in 1:git(model@control$tries)){
		if(i>1) print(paste('Try Number',i,'Best LL =',obj$value))
		errs=1
		while(errs<5){
			curr=try(optim(runif(num,-errs,errs),function(param,model=model){
								param=partition(if(spi) param else c(param,1),dim,zospad=spi,...)
								colnames(param)=c(model@control$names,'p')
								sum(git(model@control$y)*log(likfunc(model,param)))
							},control=list(fnscale=-1),model=model),silent=TRUE)
			if(class(curr)=='try-error') 
				errs=errs+1
			else {
				if(obj$value<curr$value)
					obj=curr
				break
			}
		}
		if(errs>=5) stop('Problem finding a good starting point')
	}
	param=obj$par
	value=obj$value
	param=partition(if(spi) param else c(param,1),dim,zospad=spi,...)
	return(list(param=param,value=value,num=num))
}



#Class Functions
ll=function(model,...) UseMethod('ll')
model=function(model,...) UseMethod('model')
control=function(model,...) UseMethod('control')
paramplot=function(model,...) UseMethod('paramplot')
sse=function(model)sum(resid(model)^2)
rmse=function(model)return(sqrt(sse(model)/length(model@control$y)))
condexp=function(model,...) UseMethod('condexp')
palive=function(model,...) UseMethod('palive')
spike=function(model,...) UseMethod('spike')
chitest=function(model) UseMethod('chitest')
setMethod('show','fm',print.fm)
setClass('fm',representation(control='list',raw='data.frame',nseg='numeric',param='data.frame',ll='numeric',numparam='numeric',bic='numeric'))