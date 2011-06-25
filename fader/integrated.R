count=function(raw){
	tmp=aggregate(raw$num,by=list(raw$x),sum)
	cou=rep(0,1+max(tmp[1]))
	cou[unlist(tmp[1])+1]=unlist(tmp[2])
	names(cou)=0:max(tmp[1])
	return(cou)
}
cumin=function(model){
	T=model$control$T
	maxT=max(T)
	return(vapply(0:maxT,function(t){
				sum(model$control$mult*mean(model,((t+T)>maxT)*(t+T-maxT)))
			},0))
}
cumout=function(model,T=max(model$control$T)){
	return(vapply(0:max(T),function(t){
				sum(model$control$mult*condexp(model,t))
			},0))
}
cumtracking=function(model,T=max(model$control$T)){
	first=cumin(model)
	second=cumout(model,T)
	return(c(first,second[-1]+first[length(first)]))
}

#The poisson exponential model
control.pexp=function(model,...){
	cou=count(model$raw)
	return(list(num=as.numeric(names(cou)),y=cou,x=model$raw$x,tx=model$raw$tx,T=model$raw$T,mult=model$raw$num,names=c('lambda','mu'),...))
}
ll.pexp=function(model,param=NULL,x=model$control$x){
	lambda=param[1]; mu=param[2]; T=model$control$T;tx=model$control$tx
	return(log(lambda^x*mu*exp(-(lambda+mu)*tx)/(lambda+mu)+lambda^(x+1)*exp(-(lambda+mu)*T)/(lambda+mu)))
}
indiv.pexp=function(param,model=NULL){
	lambda=param['lambda']; mu=param['mu']
	return(param['p']*vapply(model$control$num,function(x){
				sum(model$control$mult*vapply(model$control$T,function(t){
									(lambda*t)^x * exp(-(lambda+mu)*t)/factorial(x)+
											lambda^x*mu/((lambda+mu)^(x+1))*(1-exp(-(lambda+mu)*t)*sum(((lambda+mu)*t)^(1:x)/(factorial(1:x))))
								},0))
			},0))
}
predict.pexp=function(model,...){
	return(apply(apply(model$param,1,indiv.pexp,model=model),1,sum))
}
mean.pexp=function(model,t=mean(model$control$T)) return(model$param$lambda/model$param$mu-model$param$lambda/model$param$mu*exp(-model$param$mu*t))
var.pexp=function(model,t=mean(model$control$T)){
	lambda=model$param$lambda; mu=model$param$mu;
	return((lambda*(1/mu-1/mu*exp(-mu*t))+2*lambda^2*(1/mu^2-exp(-mu*t)/mu^2-t*exp(-mu*t)/mu))-mean(model,t)^2)
}
myplot.pexp=function(model,...) myplot.fm(model,x=model$control$num,...)
palive.pexp=function(model){
	lambda=model$param$lambda; mu=model$param$mu;tx=model$control$tx;T=model$control$T
	return(1/(1+(mu/(lambda+mu))*(exp((lambda+mu)*(T-tx))-1)))
}
condexp.pexp=function(model,t){
	lambda=model$param$lambda; mu=model$param$mu
	return((lambda/mu-lambda/mu*exp(-mu*t))*palive(model))
}


#The pareto NBD model
control.pnbd=function(model,...){
	cou=count(model$raw)
	return(list(num=as.numeric(names(cou)),y=cou,x=model$raw$x,tx=model$raw$tx,T=model$raw$T,mult=model$raw$num,names=c('r','alpha','s','beta'),...))
}
ll.pnbd=function(model,param=NULL,x=model$control$x){
	r=(param[1]); alp=(param[2]); s=(param[3]); bet=(param[4])
	tx=model$control$tx; T=model$control$T
	maxab = max(alp,bet)
	absab = abs(alp-bet)
	param2 = s+1
	if (alp < bet) param2 = r + x
	part1 = (alp^r*bet^s/gamma(r))*gamma(r+x)
	part2 = 1/((alp+T)^(r+x)*(bet+T)^s)
	if (absab == 0){
		F1 = 1/((maxab+tx)^(r+s+x))
		F2 = 1/((maxab+T)^(r+s+x))
	}
	else{
		F1 = hyperg_2F1(r+s+x,param2,r+s+x+1,absab/(maxab+tx))/((maxab+tx)^(r+s+x))
		F2 = hyperg_2F1(r+s+x,param2,r+s+x+1,absab/(maxab+T))/((maxab+T)^(r+s+x))
	}
	return(log(part1*(part2+(s/(r+s+x))*(F1-F2))))
}
indiv.pnbd=function(param,model=NULL){
	r=param['r']; alp=param['alpha']; s=param['s']; bet=param['beta']; 
	return(param['p']*vapply(model$control$num,function(x){
		sum(model$control$mult*vapply(model$control$T,function(t){
			A1=gamma(r+x)/(gamma(r)*factorial(x))*(alp/(alp+t))^r*(t/(alp+t))^x*(bet/(bet+t))^s
			if(alp>=bet)
				A2=alp^r*bet^s*beta(r+x,s+1)/(alp^(r+s)*beta(r,s))*hyperg_2F1(r+s,s+1,r+s+x+1,(alp-bet)/alp)
			else
				A2=alp^r*bet^s*beta(r+x,s+1)/(bet^(r+s)*beta(r,s))*hyperg_2F1(r+s,r+x,r+s+x+1,(bet-alp)/bet)
			return(A1+A2-sum(vapply(0:x,function(i){
				if(alp>=bet)
					A3=alp^r*bet^s*beta(r+x,s+1)*gamma(r+s+i)/((alp+t)^(r+s+i)*beta(r,s)*gamma(r+s))*hyperg_2F1(r+s+i,s+1,r+s+x+1,(alp-bet)/(alp+t))
				else
					A3=alp^r*bet^s*beta(r+x,s+1)*gamma(r+s+i)/((bet+t)^(r+s+i)*beta(r,s)*gamma(r+s))*hyperg_2F1(r+s+i,r+x,r+s+x+1,(bet-alp)/(bet+t))
				return(A3*t^i/factorial(i))
			},0)))
},0))},0))
}
predict.pnbd=function(model,...){
	return(apply(apply(model$param,1,indiv.pnbd,model=model),1,sum))
}
mean.pnbd=function(model,t=mean(model$control$T)){
	r=model$param$r; alp=model$param$alpha; s=model$param$s; bet=model$param$beta; 
	return(r*bet/(alp*(s-1))*(1-(bet/(bet+t))^(s-1)))
}
var.pnbd=function(model,t=mean(model$control$T)){
	r=model$param$r; alp=model$param$alpha; s=model$param$s; bet=model$param$beta; 
	return((mean(model,t)+2*r*(r+1)*bet/(alp^2*(s-1))*(bet/(s-2)-bet/(s-2)*(bet/(bet+t))^(s-2)-t*(bet/(bet+t))^(s-1)))-mean(model,t)^2)
}
myplot.pnbd=function(model,...) myplot.fm(model,x=model$control$num,...)
palive.pnbd=function(model){
	r=model$param$r; alp=model$param$alpha; s=model$param$s; bet=model$param$beta; 
	x=model$control$x;tx=model$control$tx;T=model$control$T
	
	ans=lgamma(r+x)+r*log(alp)+s*log(bet)-ll(model,unlist(model$param[1,-ncol(model$param)]))-lgamma(r)-(r+x)*log(alp+T)-s*log(bet+T)
	return(exp(ans))
	
#	
#	if(alp>bet)
#		ans=(1+s/(r+x+s)*(((alp+T)/(alp+tx))^(r+x)*((bet+T)/(alp+tx))^(s)*hyperg_2F1(r+s+x,s+1,r+s+x+1,(alp-bet)/(alp+tx))-((bet+T)/(alp+T))^s*hyperg_2F1(r+s+x,s+1,r+s+x+1,(alp-bet)/(alp+T))))^-1
#	else if(alp<bet)
#		ans=(1+s/(r+x+s)*(((alp+T)/(bet+tx))^(r+x)*((bet+T)/(bet+tx))^(s)*hyperg_2F1(r+s+x,r+x,r+s+x+1,(bet-alp)/(bet+tx))-((alp+T)/(bet+T))^(r+x)*hyperg_2F1(r+s+x,r+x,r+s+x+1,(bet-alp)/(bet+T))))^-1
#	else
#		ans=(1+s/(r+x+s)*(((alp+T)/(alp+tx)^(r+s+x))-1))^-1
	
#	maxab = max(alp,bet)
#	absab = abs(alp-bet)
#	param2 = if (alp < bet) r + x else s+1
#	if (absab == 0){
#		F1 = 1/((maxab+tx)^(r+s+x))
#		F2 = 1/((maxab+T)^(r+s+x))
#	}
#	else{
#		F1 = hyperg_2F1(r+s+x,param2,r+s+x+1,absab/(maxab+tx))/((maxab+tx)^(r+s+x))
#		F2 = hyperg_2F1(r+s+x,param2,r+s+x+1,absab/(maxab+T))/((maxab+T)^(r+s+x))
#	}
#	return((1+(s/(r+s+x))*(alp+T)^(r+x)*(bet+T)^s*(F1-F2))^(-1))
}
condexp.pnbd=function(model,t){
	r=model$param$r; alp=model$param$alpha; s=model$param$s; bet=model$param$beta; 
	T=model$control$T; x=model$control$x
	p=palive(model)
	first=((r+x)*(bet+T))/((alp+T)*(s-1))
	second=1-((bet+T)/(bet+T+t))^(s-1)
	return(p*first*second)
}


#The BG/BB model
control.bgbb=function(model,...){
	cou=count(model$raw)
	return(list(num=as.numeric(names(cou)),y=cou,x=model$raw$x,tx=model$raw$tx,T=model$raw$T,mult=model$raw$num,n=mean(model$raw$T),names=c('alpha','beta','gamma','delta'),...))
}
ll.bgbb=function(model,param=NULL,x=model$control$x){
	a = param[1]; b = param[2]; g = param[3]; d = param[4]; 
	T = model$control$T; x = model$control$x; tx = model$control$tx;
	denom_ab = lbeta(a,b); denom_gd = lbeta(g,d);
	lik = exp(lbeta(a+x, b+T-x) - denom_ab + lbeta(g,d+T) - denom_gd);
	count = T - tx - 1;
	lik=lik+vapply(1:length(lik),function(j){
				ifelse(count[j]>=0,sum(vapply(0:count[j],function(i){exp(lbeta(a+x[j],b+tx[j]-x[j]+i) - denom_ab + lbeta(g+1, d+tx[j]+i) - denom_gd)},0)),0)
			},0)
	return (log(lik))
}
indiv.bgbb=function(param,model=NULL){
	a=param['alpha']; b=param['beta']; g=param['gamma']; d=param['delta'] 
	return(param['p']*vapply(model$control$num,function(x){
						sum(model$control$mult*vapply(model$control$T,function(n){
											i=x:(n-1)
											return(choose(n,x)*beta(a+x,b+n-x)*beta(g,d+n)/(beta(a,b)*beta(g,d))+ifelse(x>n-1,0,sum(choose(i,x)*beta(a+x,b+i-x)*beta(g+1,d+i)/(beta(a,b)*beta(g,d)))))
											
										},0))
			},0))
}
predict.bgbb=function(model,...){
	return(apply(apply(model$param,1,indiv.bgbb,model=model),1,sum))
}
mean.bgbb=function(model,n = model$control$n){
	a=model$param$a; b=model$param$b; g=model$param$g; d=model$param$d 
	return((a/(a+b))*(d/(g-1))*(1-gamma(d+g)/gamma(d+g+n)*gamma(1+d+n)/gamma(1+d)))
}
myplot.bgbb=function(model,...) myplot.fm(model,x=model$control$num,...)
condexp.bgbb = function(model,n2) {
	a=model$param$a; b=model$param$b; g=model$param$g; d=model$param$d 
	n = model$control$n; x = model$control$x
	logsum=-ll(model,unlist(model$param[1,-ncol(model$param)]))+lbeta(a+x+1,b+n-x) - lbeta(a,b)+ lgamma(g+d) - lgamma(1+d)
	return(exp(logsum)*d/(g-1)*(gamma(1+d+n)/gamma(g+d+n)-gamma(1+d+n+n2)/gamma(g+d+n+n2)))
}


#The BG/NBD model
control.bgnbd=function(model,...){
	cou=count(model$raw)
	return(list(num=as.numeric(names(cou)),y=cou,x=model$raw$x,tx=model$raw$tx,T=model$raw$T,mult=model$raw$num,names=c('r','alpha','a','b'),...))
}
ll.bgnbd=function(model,param=NULL,x=model$control$x){
	tx=model$control$tx; T=model$control$T
	r=(param[1]); alp=(param[2]); a=(param[3]); b=(param[4])
	lA1=lgamma(r+x)+r*log(alp)-lgamma(r)
	lA2=lgamma(a+b)+lgamma(b+x)-lgamma(b)-lgamma(a+b+x)
	lA3=(r+x)*log(1/(alp+T))
	lA4=log(a/(b+x-1))+(r+x)*log(1/(alp+tx))
	return(lA1+lA2+log(exp(lA3)+(x>0)*exp(lA4)))
}
indiv.bgnbd=function(param,model=NULL){
	r=param['r']; alp=param['alpha']; a=param['a']; b=param['b']; 
	return(param['p']*vapply(model$control$num,function(x){
						sum(model$control$mult*vapply(model$control$T,function(t){
											exp(log(beta(a,b+x))+lgamma(r+x)+r*log(alp/(alp+t))+x*log(t/(alp+t))-log(beta(a,b))-lgamma(r)-lfactorial(x))+
													(x>0)*beta(a+1,b+x-1)/(beta(a,b))*
													(1-(alp/(alp+t))^r*ifelse(x-1>=0,sum(vapply(0:(x-1),function(j){
																	exp(lgamma(r+j)+j*log(t/(alp+t))-lgamma(r)-lfactorial(j))
																},0)),0))
										}))},0))
}
predict.bgnbd=function(model,...){
	return(apply(apply(model$param,1,indiv.bgnbd,model=model),1,sum))
}
mean.bgnbd=function(model,t=mean(model$control$T)){
	r=model$param$r; alp=model$param$alpha; a=model$param$a; b=model$param$b; 
	return((a+b-1)/(a-1)*(1-(alp/(alp+t))^r*hyperg_2F1(r,b,a+b-1,t/(alp+t))))
}
myplot.bgnbd=function(model,...) myplot.fm(model,x=model$control$num,...)
palive.bgnbd=function(model,t=mean(model$control$T)){
	r=model$param$r; alp=model$param$alpha; a=model$param$a; b=model$param$b; 
	x=model$control$x; tx=model$control$tx; T=model$control$T
	return(1/(1+(x>0)*(a/(b+x-1))*((alp+T)/(alp+tx))^(r+x)))
}
condexp.bgnbd=function(model,t){
	r=model$param$r; alp=model$param$alpha; a=model$param$a; b=model$param$b;
	x=model$control$x; T=model$control$T; tx=model$control$tx
	return((a+b+x-1)/(a-1)*(1-((alp+T)/(alp+T+t))^(r+x)*hyperg_2F1(r+x,b+x,a+b+x-1,t/(alp+T+t)))/(1+(x>0)*a/(b+x-1)*((alp+T)/(alp+tx))^(r+x)))
}