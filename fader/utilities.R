#Standard Functions
standardresid=function(model) predict(model)-model$control$y
standardprint=function(model) {print(paste(toupper(class(model)),'Model LL =',model$ll));print(model$param)}
standardplot=function(model,...) {
	fmhistoplot(model,model$control$y,...)
}
standardmodel=function(model,names,len=length(names),...){
	obj=findparam(model,len,model$nseg,...)
	colnames(obj$param)=names
	return(obj)
}
standardpredict=function(model,num=model$control$num,...){
	param=model$param
	val=apply(apply(param,1,weightedlik,model=model,...),1,sum)
	return(num*val)
}
standardtest=function(filename,model,nseg=1,...){
	data=read.csv(r(filename))
	print('Model')
	mod=fm(data,model,nseg,...)
	print(mod)
	print('Prediction/RMSE')
	print(data.frame(act=mod$control$y,mod=predict(mod)))
	print(rmse(mod))
	myplot(mod)
	print('Mean/Var')
	print(mean(mod))
	print(var(mod))
	return(mod)
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

#Helper Functions
fm=function(data,modname='bb',nseg=1,...){
	a=list(raw=data,nseg=nseg)
	class(a)=modname
	a$control=control(a,...)
	obj=model(a)
	a$param=obj$param
	a$ll=obj$value
	return(a)
}
partition=function(params,by=1,pos=1:by,zeroone=c(),zeroonesum=by){
	tmp=data.frame(matrix(params,ncol=by))
	colnames(tmp)=paste('x',1:by,sep='')
	if(length(union(union(pos,zeroone),zeroonesum))>0) tmp[union(union(pos,zeroone),zeroonesum)]=exp(tmp[union(union(pos,zeroone),zeroonesum)])
	if(length(zeroone)>0) tmp[zeroone]=tmp[zeroone]/(1+tmp[zeroone])
	if(length(zeroonesum)>0) tmp[zeroonesum]=tmp[zeroonesum]/sum(tmp[zeroonesum])
	return(tmp)
}
weightedlik=function(param,model=NULL,...){
	logl=ll(model,param=param[-length(param)],...)
	return(param[length(param)]*exp(logl))
}
findparam=function(model,dim=1,nseg=1,...){
	obj=list(value=-Inf)
	tries=if(is.null(model$control$tries)) 1 else model$control$tries
	for(i in 1:tries){
		curr=optim(runif(dim*nseg-1,-1,1),function(param,model=model){
				param=partition(c(param,1),dim,...)
				sum(model$control$mult*log(apply(apply(param,1,weightedlik,model=model),1,sum)))
			},control=list(fnscale=-1),model=model)
		if(obj$value<curr$value)
			obj=curr
	}
	param=obj$par
	value=obj$value
	param=partition(c(param,1),dim,...)
	return(list(param=param,value=value))
}
fmhistoplot=function(model,counts,x=model$control$x,...){
	len=length(x)
	breaks=c(min(x)-0.5,(x[-1]+x[-length(x)])/2,max(x)+0.5)
	add=FALSE
	#if(length(model$control$x)==len && all(x==model$control$x)){
		add=TRUE
		myhist=list(breaks=breaks,counts=counts[1:len])
		class(myhist)='histogram'
		plot(myhist,col='black')
	#}else
	myhist2=list(breaks=breaks,counts=predict(model,x=x,...)[1:len])
	class(myhist2)='histogram'
	plot(myhist2,add=add,col=rgb(1,0,0,.5))
}

#Class Functions
var=function(model,...) UseMethod('var')
ll=function(model,...) UseMethod('ll')
model=function(model,...) UseMethod('model')
control=function(model,...) UseMethod('control')
myplot=function(model,...) UseMethod('myplot')
sse=function(model)sum(resid(model)^2)
rmse=function(model)return(sqrt(sse(model)/length(model$control$y)))