setwd('C:\\Users\\Andrew\\Desktop\\Workspace\\models\\fader')
sse=function(model){
	sum(resid(model)^2)
}
rmse=function(model){
	return(sqrt(sse(model)/model$control$len))
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
partition=function(params,by=1){
	tmp=data.frame(matrix(params,ncol=by))
	colnames(tmp)=paste('x',1:by,sep='')
	return(tmp)
}
specialpartition=function(params,by=2){
	param=partition(params,by)
	param=exp(param)
	param[ncol(param)]=param[ncol(param)]/sum(param[ncol(param)])
	return(param)
}
specialpartition2=function(params,by=2){
	param=partition(params,by)
	param=exp(param)
	param[-ncol(param)]=param[-ncol(param)]/(1+param[-ncol(param)])
	param[ncol(param)]=param[ncol(param)]/sum(param[ncol(param)])
	return(param)
}
specialpartition3=function(params,by=2){
	param=partition(params,by)
	param=exp(param)
	param[1]=param[1]/(1+param[1])
	param[ncol(param)]=param[ncol(param)]/sum(param[ncol(param)])
	return(param)
}

fmhistoplot=function(model,counts,...){
	x=model$control$x
	len=length(x)
	breaks=c(min(x)-0.5,(x[-1]+x[-length(x)])/2,max(x)+0.5)
	myhist=list(breaks=breaks,counts=counts[1:len])
	myhist2=list(breaks=breaks,counts=predict(model,...)[1:len])
	class(myhist)='histogram'
	class(myhist2)='histogram'
	plot(myhist,col='black')
	plot(myhist2,add=TRUE,col=rgb(1,0,0,.5))
}


fm=function(data,modname='bb',nseg=1,...){
	a=list(raw=data,nseg=nseg)
	class(a)=modname
	a$control=control(a,...)
	a$param=model(a,nseg)
	return(a)
}
weightedlik=function(param,model,...){
	logl=ll(model,param=param[-length(param)],...)
	return(param[length(param)]*exp(logl))
}
findparam=function(model,dim=1,nseg=1,func=specialpartition){
	param=optim(runif(dim*nseg),function(param,model=model){
				param=func(param,dim)
				sum(model$control$mult*log(apply(apply(param,1,weightedlik,model=model),1,sum)))
			},control=list(fnscale=-1),model=model)$par
	param=func(param,dim)
	return(param)
}

cov=function(x,y,...) UseMethod("cov")
var=function(x,...) UseMethod('var')
ll=function(x,...) UseMethod('ll')
model=function(x,...) UseMethod('model')
control=function(x,...) UseMethod('control')
#var,cov,mean,sd,resid,plot,predict,print

#bg,wg,geom go with time.csv
#pois,nbd go with count.csv
