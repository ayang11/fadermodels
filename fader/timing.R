getChurn=function(data,len=length(data)){
	return(data[1:(len-1)]-data[2:len])
}

#dProb is death probability
#sProb is survival probability
#lProb is living probability

#The geometric model is a discrete one parameter timing model. It answers when people will stop ordering. 
ll.geom=function(param,data=data){
	len=length(data)
	loss=c(getChurn(data,len),data[len])
	dProb=pgeom(0:(len-2),param)
	sProb=1-pgeom(len-2,param)
	return(sum(loss*log(c(dProb,sProb))))
}
model.geom=function(data){
	param=optimize(geomLL,c(0,1),data=data,maximum=TRUE)
	return(param$maximum)
}
predict.geom=function(param,data=data,n=data[1]){
	len=length(data)
	return(data.frame(act=data,model=n*c(1,1-pgeom(0:(len-2),param))))
}
mean.geom=function(param){return(1/param)}
var.geom=function(param){return((1-param)/param^2)}


#The Beta Geometric Model is a discrete two parameter timing model. It answers when people will stop ordering. 
ll.bg=function(param,data=data){
	alp=exp(param[1])
	bet=exp(param[2])
	len=length(data)
	dProb=beta(alp+1,bet+(1:(len-1))-1)/beta(bet,alp)
	sProb=beta(alp,bet+len-1)/beta(alp,bet)
	loss=c(getChurn(data,len),data[len])
	return(sum(loss*log(c(dProb,sProb))))
}
model.bg=function(data){
	param=optim(runif(2),bgLL,data=data,control=list(fnscale=-1))
	return(exp(param$par))
}
predict.bg=function(param,data=data,n=data[1]){
	alp=param[1]
	bet=param[2]
	len=length(data)
	lProb=1-cumsum(beta(alp+1,bet+(1:(len-1))-1)/beta(bet,alp))
	return(data.frame(act=data,model=n*c(1,lProb)))
}
mean.bg=function(param){
	alp=param[1]
	bet=param[2]
	return(bet/(alp-1))
}
var.bg=function(param){
	alp=param[1]
	bet=param[2]
	alp*bet*(alp+bet-1)/((alp-1)^2*(alp-2))
}


#The discrete weibull is a discrete two parameter timing model. It answers when people will stop ordering. 
ll.dw=function(param,data=data){
	theta=exp(param[1])/(1+exp(param[1]))
	k=exp(param[2])
	len=length(data)
	sProb=(1-theta)^((0:len)^k)
	dProb=sProb[1:(len-1)]-sProb[2:(len)]
	loss=c(getChurn(data,len),data[len])
	return(sum(loss*log(c(dProb,sProb[len]))))
}
model.dw=function(data){
	param=optim(runif(2),dwLL,data=data,control=list(fnscale=-1))$par
	param=exp(param)
	param[1]=param[1]/(1+param[1])
	return(param)
}
predict.dw=function(param,data=data,n=data[1]){
	theta=param[1]
	k=param[2]
	len=length(data)-1
	sProb=(1-theta)^((0:len)^k)
	return(data.frame(act=data,model=n*sProb))
}

#The Beta discrete weibull is a discrete three parameter timing model. It answers when people will stop ordering. 
ll.bdw=function(param,data=data){
	alp=exp(param[1])
	bet=exp(param[2])
	k=exp(param[3])
	len=length(data)
	sProb=c(1,beta(alp,bet+(1:len)^k)/beta(alp,bet))
	dProb=sProb[1:(len-1)]-sProb[2:(len)]
	loss=c(getChurn(data,len),data[len])
	return(sum(loss*log(c(dProb,sProb[len]))))
}
model.bdw=function(data){
	param=optim(runif(3),bdwLL,data=data,control=list(fnscale=-1))
	return(exp(param$par))
}
predict.bdw=function(param,data=data,n=data[1]){
	alp=param[1]
	bet=param[2]
	k=param[3]
	len=length(data)-1
	sProb=c(1,beta(alp,bet+(1:len)^k)/beta(alp,bet))
	return(data.frame(act=data,model=n*sProb))
}