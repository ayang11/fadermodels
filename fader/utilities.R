sse=function(model){
	sum(resid(model)^2)
}
rmse=function(model){
	return(sqrt(sse(model)/model$control$len))
}

standardresid=function(model) predict(model)-model$control$y
standardprint=function(model) print(model$param)
standardplot=function(model,...) fmhistoplot(model,model$control$y,...)
standardmodel=function(model,names,nseg,len=length(names),...){
	param=findparam(model,len,nseg,...)
	colnames(param)=names
	return(param)
}
standardpredict=function(model,num=model$control$num,...){
	param=model$param
	val=apply(apply(param,1,weightedlik,model=model,...),1,sum)
	return(num*val)
}

getxy=function(raw,title){
	cols=ncol(raw)
	if(is.null(cols)||cols==1){
		y=unlist(raw)
		x=0:(length(y)-1)
	}else if(cols==2){
		y=raw[,2]
		x=raw[,1]
	}else throw(paste('Error in the data format for the',title))
	return(list(x=x,y=y))
}
partition=function(params,by=1,pos=1:by,bound=c(),prob=by){
	tmp=data.frame(matrix(params,ncol=by))
	colnames(tmp)=paste('x',1:by,sep='')
	if(length(union(union(pos,bound),prob))>0) tmp[union(union(pos,bound),prob)]=exp(tmp[union(union(pos,bound),prob)])
	if(length(bound)>0) tmp[bound]=tmp[bound]/(1+tmp[bound])
	if(length(prob)>0) tmp[prob]=tmp[prob]/sum(tmp[prob])
	return(tmp)
}
fmhistoplot=function(model,counts,x=model$control$x,...){
	len=length(x)
	breaks=c(min(x)-0.5,(x[-1]+x[-length(x)])/2,max(x)+0.5)
	add=FALSE
	if(all(x==model$control$x)){
		add=TRUE
		myhist=list(breaks=breaks,counts=counts[1:len])
		class(myhist)='histogram'
		plot(myhist,col='black')
	}
	myhist2=list(breaks=breaks,counts=predict(model,x=x,...)[1:len])
	class(myhist2)='histogram'
	plot(myhist2,add=add,col=rgb(1,0,0,.5))
}

fm=function(data,modname='bb',nseg=1,...){
	a=list(raw=data,nseg=nseg)
	class(a)=modname
	a$control=control(a,...)
	a$param=model(a,nseg)
	return(a)
}
weightedlik=function(param,model=NULL,...){
	logl=ll(model,param=param[-length(param)],...)
	return(param[length(param)]*exp(logl))
}
findparam=function(model,dim=1,nseg=1,...){
	param=optim(runif(dim*nseg),function(param,model=model){
				param=partition(param,dim,...)
				sum(model$control$mult*log(apply(apply(param,1,weightedlik,model=model),1,sum)))
			},control=list(fnscale=-1),model=model)$par
	param=partition(param,dim,...)
	return(param)
}

var=function(model,...) UseMethod('var')
ll=function(model,...) UseMethod('ll')
model=function(model,...) UseMethod('model')
control=function(model,...) UseMethod('control')
myplot=function(model,...) UseMethod('myplot')