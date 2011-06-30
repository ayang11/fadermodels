#Bayesian posterior probabilities
post=function(model,...){
	lik=apply(model@param,1,weightedlik,model=model,...)
	colnames(lik)=c(1:ncol(lik))
	if(hasSpike(model)){
		lik=cbind(spike(model,model@param),lik)
		colnames(lik)=c('Spike',1:(ncol(lik)-1))
	}
	total=rowSums(lik)
	return(lik/total)
}

#Model Validation Tests

#This performs a LL ratio test. Be sure model 1 is nested within model 2
#The null hypothesis is that the two models are identical
lltest=function(nullmodel,altmodel) {
	if(all(class(nullmodel)==class(altmodel))&&isTRUE(all.equal(nullmodel@raw==altmodel@raw)))
		return(1-pchisq(-2*nullmodel@ll+2*altmodel@ll,altmodel@numparam-nullmodel@numparam))
	stop('The models are not nested')
	
}
#This performs a chi squared goodness of fit test
#The null hypothesis is that the two count data sets are from the same distribution.
chitest.fm=function(model,mod=predict(model),act=model@control$y,...) {
	return(chisq.test(mod,p=prop.table(act),...))
}
bictest=function(...) {
	loc=which.min(vapply(list(...),function(x) x@bic,0))
	print(paste('Model',loc,'has the best BIC value'))
	return(list(...)[[loc]])
}

#Helper Plot functions
plotGamma=function(shape,rate,p,bounds=c(0,10),points=100,type='l',spikep=0,...){
	x=ppoints(points)*(bounds[2]-bounds[1])+bounds[1]
	y=rowSums(vapply(1:length(shape),function(i){p[i]*dgamma(x,shape[i],rate[i])},x))
	y[x==0]=y[x==0]+spikep
	plot(x,y,type=type,...)
}
plotBeta=function(sh1,sh2,p,bounds=c(0,1),points=100,type='l',spikep=0,...){
	x=ppoints(points)*(bounds[2]-bounds[1])+bounds[1]
	y=rowSums(vapply(1:length(sh1),function(i){p[i]*dbeta(x,sh1[i],sh2[i])},x))
	y[x==0]=y[x==0]+spikep
	plot(x,y,type=type,...)
}
plotSpike=function(val,prob,lwd=5,xlim=c(0,max(val)),ylim=c(0,max(c(1,prob))),...){
	plot(0,0,type='n',xlim=xlim,ylim=ylim,...)
	segments(val,0,val,prob,lwd=lwd,...)
}

#Other Helper Functions
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
checkdata=function(model,condition=NULL){
	if(!(all(condition %in% colnames(model@raw))))
		stop(paste(toupper(class(model)),'must have column headers:',paste(condition,collapse=' '),'in the data file'))
}
hasSpike=function(model) return(git(model@control$allowspike,FALSE)&&git(model@control$spike,FALSE))
git=function(val,default=1) return(if(is.null(val)) default else val)