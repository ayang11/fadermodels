#Standard Functions
var.fm=function(model) return(NA)
mean.fm=function(model) return(NA)
spike.fm=function(model,param,x=model$control$x) return((x==0)*(1-sum(param$p)))
residuals.fm=function(model,...) predict(model,...)-model$control$y
print.fm=function(model) {
	print(paste(toupper(class(model)[1]),'Model LL =',model$ll))
	if(git(model$control$allowspike,FALSE)&&git(model$control$spike,FALSE))
		print(paste('Spike P =',1-sum(model$param$p)))
	print(model$param)
}
barplot.fm=function(model,x=model$control$x,y=model$control$y,legend.text=TRUE,...) {
	mat=matrix(c(y,predict(model,x=x)),nrow=2,byrow=TRUE)
	rownames(mat)=c('Act','Model')
	colnames(mat)=x[1:ncol(mat)]
	barplot(mat,beside=TRUE,legend.text=legend.text,...)
}
model.fm=function(model,...){
	names=c(model$control$names,'p')
	len=length(names)
	obj=findparam(model,len,model$nseg,...)
	colnames(obj$param)=names
	return(obj)
}
predict.fm=function(model,num=git(model$control$num),...){
	param=model$param
	spi=git(model$control$allowspike,FALSE)&&git(model$control$spike,FALSE)
	val=(if(spi) spike(model,param) else 0) +apply(apply(param,1,weightedlik,model=model,...),1,sum)
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
	barplot(mod)
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
git=function(val,default=1) return(if(is.null(val)) default else val)
#reserved variables 
#names is the name of the parameters. This is required.
#x is what will be plotted on the x axis
#y is what will be plotted on the y axis
#num is the number of people (used for prediction from probabilities).
#mult is the number of people in each segment of an observation. It is multiplied by likelihood.
#allowspike is whether or not the model allows spikes.
strip=function(a,names=NULL,x=NULL,y=NULL,num=NULL,mult=NULL,allowspike=NULL,...) return(list(a,...))
fm=function(data,modname='bb',nseg=1,...){
	a=list(raw=data,nseg=nseg)
	class(a)=c(modname,'fm')
	a$control=do.call(control,strip(a,...))
	obj=model(a)
	a$param=obj$param
	a$ll=obj$value
	return(a)
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
	spi=git(model$control$allowspike,FALSE)&&git(model$control$spike,FALSE)
	for(i in 1:git(model$control$tries)){
		if(i>1) print(paste('Try Number',i,'Best LL =',obj$value))
		errs=1
		while(errs<5){
			curr=try(optim(runif(dim*nseg-1*!spi,-errs,errs),function(param,model=model){
								param=partition(if(spi) param else c(param,1),dim,zospad=spi,...)
								colnames(param)=c(model$control$names,'p')
								sum(git(model$control$mult)*log((if(spi) spike(model,param) else 0) +apply(apply(param,1,weightedlik,model=model),1,sum)))
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
	return(list(param=param,value=value))
}


#Class Functions
var=function(model,...) UseMethod('var')
ll=function(model,...) UseMethod('ll')
model=function(model,...) UseMethod('model')
control=function(model,...) UseMethod('control')
myplot=function(model,...) UseMethod('myplot')
sse=function(model)sum(resid(model)^2)
rmse=function(model)return(sqrt(sse(model)/length(model$control$y)))
condexp=function(model,...) UseMethod('condexp')
palive=function(model,...) UseMethod('palive')
spike=function(model,...) UseMethod('spike')