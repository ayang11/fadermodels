#mean/variance and dont work for current models. Fix that
cumin=function(model,T=model@control$T,points=40,x=(0:points)/points*max(T)){
	maxT=max(T)
	val=(vapply(x,function(t){
							return(sum(model@control$y*vapply(((t+T)>maxT)*(t+T-maxT),function(i){
														return(sum(model@param$p*mean(model,i)))
													},0)))
						},0))
	names(val)=x
	return(val)
}
cumout=function(model,T=max(model@control$T),points=40,x=(0:points)/points*max(T)){
	val=(vapply(x,function(t){
							sum(model@control$y*condexp(model,t))
						},0))
	names(val)=x+T
	return(val)
}
cumtracking=function(model,T.out=max(model@control$T)){
	first=cumin(model)
	second=cumout(model,T.out)
	val=(c(first,second[-1]+first[length(first)]))
	return(val)
}


#Integrated models class
setClass('integ',contains='fm')
control.integ=function(model,names=names,...){
	cou=tapply(model@raw$num,model@raw$x,sum)
	return(list(num=sum(cou),x.total=as.numeric(names(cou)),plot.y=cou,x=model@raw$x,tx=model@raw$tx,T=model@raw$T,y=model@raw$num,names=names,...))
}
predict.integ=function(model,...){
	return(apply(apply(model@param,1,indiv,model=model,...),1,sum))
}
chitest.integ=function(model,...) chitest.fm(model,act=model@control$plot.y)
barplot.integ=function(model,...) barplot.fm(model,x=model@control$x.total,act=model@control$plot.y,...)
residuals.integ=function(model,...) residuals.fm(model,y=model@control$plot.y,...)


#The poisson exponential model
setClass('pexp',contains='integ')
control.pexp=function(model,...) control.integ(model,names=c('lambda','mu'),...)
ll.pexp=function(model,param=NULL,x=model@control$x){
	lambda=param[1]; mu=param[2]; T=model@control$T;tx=model@control$tx
	return(log(lambda^x*mu*exp(-(lambda+mu)*tx)/(lambda+mu)+lambda^(x+1)*exp(-(lambda+mu)*T)/(lambda+mu)))
}
indiv.pexp=function(param,model=NULL,x.total=model@control$x.total,T=model@control$T){
	lambda=param['lambda']; mu=param['mu']
	return(param['p']*vapply(x.total,function(x.ind){
						sum(model@control$y*vapply(T,function(t){
											(lambda*t)^x.ind * exp(-(lambda+mu)*t)/factorial(x.ind)+
													lambda^x.ind*mu/((lambda+mu)^(x.ind+1))*(1-exp(-(lambda+mu)*t)*sum(((lambda+mu)*t)^(1:x.ind)/(factorial(1:x.ind))))
										},0))
					},0))
}
mean.pexp=function(model,t=mean(model@control$T)) {
	return(model@param$lambda/model@param$mu-model@param$lambda/model@param$mu*exp(-model@param$mu*t))
}
vcov.pexp=function(model,t=mean(model@control$T)){
	lambda=model@param$lambda; mu=model@param$mu;
	return((lambda*(1/mu-1/mu*exp(-mu*t))+2*lambda^2*(1/mu^2-exp(-mu*t)/mu^2-t*exp(-mu*t)/mu))-mean(model,t)^2)
}
palive.pexp=function(model,all.lambda=model@param$lambda, all.mu=model@param$mu,all.p=model@param$p){
	x=model@control$x; tx=model@control$tx; T=model@control$T
	return(apply(vapply(1:length(all.lambda),function(i){
								lambda=all.lambda[i];mu=all.mu[i];p=all.p[i]
								return(p*(1/(1+(mu/(lambda+mu))*(exp((lambda+mu)*(T-tx))-1))))
							},rep(0,length(x))),1,sum))
}
condexp.pexp=function(model,t){
	all.lambda=model@param$lambda; all.mu=model@param$mu; all.p=model@param$p
	x=model@control$x; tx=model@control$tx; T=model@control$T
	return(apply(vapply(1:length(all.lambda),function(i){
								lambda=all.lambda[i];mu=all.mu[i];p=all.p[i]
								return(p*(lambda/mu-lambda/mu*exp(-mu*t))*palive(model,lambda,mu,1))
							},rep(0,length(x))),1,sum))
}
paramplot.pexp=function(model,...) {
	par(mfrow=c(1,2))
	plotSpike(model@param$lambda,model@param$p,main='Purchase Process')
	plotSpike(model@param$mu,model@param$p,main='Death Process')
	par(mfrow=c(1,1))
}


#The pareto NBD model
setClass('pnbd',contains='integ')
control.pnbd=function(model,...) control.integ(model,names=c('r','alpha','s','beta'),...)
ll.pnbd=function(model,param=NULL,x=model@control$x){
	r=(param[1]); alp=(param[2]); s=(param[3]); bet=(param[4])
	tx=model@control$tx; T=model@control$T
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
indiv.pnbd=function(param,model=NULL,x.total=model@control$x.total,T=model@control$T){
	r=param['r']; alp=param['alpha']; s=param['s']; bet=param['beta']; 
	return(param['p']*vapply(x.total,function(x){
						sum(model@control$y*vapply(T,function(t){
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
mean.pnbd=function(model,t=mean(model@control$T)){
	r=model@param$r; alp=model@param$alpha; s=model@param$s; bet=model@param$beta; 
	return(r*bet/(alp*(s-1))*(1-(bet/(bet+t))^(s-1)))
}
vcov.pnbd=function(model,t=mean(model@control$T)){
	r=model@param$r; alp=model@param$alpha; s=model@param$s; bet=model@param$beta; 
	return((mean(model,t)+2*r*(r+1)*bet/(alp^2*(s-1))*(bet/(s-2)-bet/(s-2)*(bet/(bet+t))^(s-2)-t*(bet/(bet+t))^(s-1)))-mean(model,t)^2)
}
palive.pnbd=function(model,all.r=model@param$r, all.alp=model@param$alpha, all.s=model@param$s, all.bet=model@param$beta,all.p=model@param$p){
	x=model@control$x; tx=model@control$tx; T=model@control$T
	return(apply(vapply(1:length(all.r),function(i){
								r=all.r[i];alp=all.alp[i];s=all.s[i];bet=all.bet[i];p=all.p[i]
								maxab = max(alp,bet)
								absab = abs(alp-bet)
								param2 = if (alp < bet) r + x else s+1
								if (absab == 0){
									F1 = 1/((maxab+tx)^(r+s+x))
									F2 = 1/((maxab+T)^(r+s+x))
								}
								else{
									F1 = hyperg_2F1(r+s+x,param2,r+s+x+1,absab/(maxab+tx))/((maxab+tx)^(r+s+x))
									F2 = hyperg_2F1(r+s+x,param2,r+s+x+1,absab/(maxab+T))/((maxab+T)^(r+s+x))
								}
								return(p*((1+(s/(r+s+x))*(alp+T)^(r+x)*(bet+T)^s*(F1-F2))^(-1)))
							},rep(0,length(x))),1,sum))
}

condexp.pnbd=function(model,t){
	all.r=model@param$r; all.alp=model@param$alpha; all.s=model@param$s; all.bet=model@param$beta; all.p=model@param$p
	x=model@control$x; tx=model@control$tx; T=model@control$T
	return(apply(vapply(1:length(all.r),function(i){
								r=all.r[i];alp=all.alp[i];s=all.s[i];bet=all.bet[i];p=all.p[i]
								first=((r+x)*(bet+T))/((alp+T)*(s-1))
								second=1-((bet+T)/(bet+T+t))^(s-1)
								return(p*palive(model,r,alp,s,bet,1)*first*second)
							},rep(0,length(x))),1,sum))
}
paramplot.pnbd=function(model,...) {
	par(mfrow=c(1,2))
	plotGamma(model@param$r,model@param$alpha,model@param$p,main='Purchase Process')
	plotGamma(model@param$s,model@param$beta,model@param$p,main='Death Process')
	par(mfrow=c(1,1))
}


#The BG/BB model
setClass('bgbb',contains='integ')
control.bgbb=function(model,...) control.integ(model,names=c('alpha','beta','gamma','delta'),n=mean(model@raw$T),...)
ll.bgbb=function(model,param=NULL,x=model@control$x){
	a = param[1]; b = param[2]; g = param[3]; d = param[4]; 
	T = model@control$T; x = model@control$x; tx = model@control$tx;
	denom_ab = lbeta(a,b); denom_gd = lbeta(g,d);
	lik = exp(lbeta(a+x, b+T-x) - denom_ab + lbeta(g,d+T) - denom_gd);
	count = T - tx - 1;
	lik=lik+vapply(1:length(lik),function(j){
				(if(count[j]>=0)sum(vapply(0:count[j],function(i){exp(lbeta(a+x[j],b+tx[j]-x[j]+i) - denom_ab + lbeta(g+1, d+tx[j]+i) - denom_gd)},0)) else 0)
			},0)
	return (log(lik))
}
indiv.bgbb=function(param,model=NULL,x.total=model@control$x.total,T=model@control$T){
	a=param['alpha']; b=param['beta']; g=param['gamma']; d=param['delta'] 
	return(param['p']*vapply(x.total,function(x){
						sum(model@control$y*vapply(T,function(n){
											i=x:(n-1)
											return(choose(n,x)*beta(a+x,b+n-x)*beta(g,d+n)/(beta(a,b)*beta(g,d))+(if(x>n-1) 0 else sum(choose(i,x)*beta(a+x,b+i-x)*beta(g+1,d+i)/(beta(a,b)*beta(g,d)))))
											
										},0))
					},0))
}
mean.bgbb=function(model,n = model@control$n){
	a=model@param$a; b=model@param$b; g=model@param$g; d=model@param$d 
	return((a/(a+b))*(d/(g-1))*(1-gamma(d+g)/gamma(d+g+n)*gamma(1+d+n)/gamma(1+d)))
}
condexp.bgbb = function(model,n2) {
	all.a=model@param$a; all.b=model@param$b; all.g=model@param$g; all.d=model@param$d; all.p=model@param$p
	n = model@control$n; x=model@control$x
	return(apply(vapply(1:length(all.a),function(i){
								a=all.a[i];b=all.b[i];g=all.g[i];d=all.d[i];p=all.p[i]
								logsum=-ll(model,unlist(model@param[1,-ncol(model@param)]))+lbeta(a+x+1,b+n-x) - lbeta(a,b)+ lgamma(g+d) - lgamma(1+d)
								return(p*(exp(logsum)*d/(g-1)*(gamma(1+d+n)/gamma(g+d+n)-gamma(1+d+n+n2)/gamma(g+d+n+n2))))
							},rep(0,length(x))),1,sum))
}
paramplot.bgbb=function(model,...) {
	par(mfrow=c(1,2))
	plotBeta(model@param$alpha,model@param$beta,model@param$p,main='Purchase Process')
	plotBeta(model@param$gamma,model@param$delta,model@param$p,main='Death Process')
	par(mfrow=c(1,1))
}


#The BG/NBD model
setClass('bgnbd',contains='integ')
control.bgnbd=function(model,...) control.integ(model,names=c('r','alpha','a','b'),...)
ll.bgnbd=function(model,param=NULL,x=model@control$x){
	tx=model@control$tx; T=model@control$T
	r=(param[1]); alp=(param[2]); a=(param[3]); b=(param[4])
	lA1=lgamma(r+x)+r*log(alp)-lgamma(r)
	lA2=lgamma(a+b)+lgamma(b+x)-lgamma(b)-lgamma(a+b+x)
	lA3=(r+x)*log(1/(alp+T))
	lA4=log(a/(b+x-1))+(r+x)*log(1/(alp+tx))
	return(lA1+lA2+log(exp(lA3)+(x>0)*exp(lA4)))
}
indiv.bgnbd=function(param,model=NULL,x.total=model@control$x.total,T=model@control$T){
	r=param['r']; alp=param['alpha']; a=param['a']; b=param['b']; 
	return(param['p']*vapply(x.total,function(x){
						sum(model@control$y*vapply(T,function(t){
											exp(log(beta(a,b+x))+lgamma(r+x)+r*log(alp/(alp+t))+x*log(t/(alp+t))-log(beta(a,b))-lgamma(r)-lfactorial(x))+
													(x>0)*beta(a+1,b+x-1)/(beta(a,b))*
													(1-(alp/(alp+t))^r*(if(x-1>=0)sum(vapply(0:(x-1),function(j){
																					exp(lgamma(r+j)+j*log(t/(alp+t))-lgamma(r)-lfactorial(j))
																				},0)) else 0))
										},0))},0))
}
mean.bgnbd=function(model,t=mean(model@control$T)){
	r=model@param$r; alp=model@param$alpha; a=model@param$a; b=model@param$b; 
	return((a+b-1)/(a-1)*(1-(alp/(alp+t))^r*hyperg_2F1(r,b,a+b-1,t/(alp+t))))
}
palive.bgnbd=function(model,t=mean(model@control$T)){
	all.r=model@param$r; all.alp=model@param$alpha; all.a=model@param$a; all.b=model@param$b; all.p=model@param$p
	x=model@control$x; tx=model@control$tx; T=model@control$T
	return(apply(vapply(1:length(all.r),function(i){
								r=all.r[i];alp=all.alp[i];a=all.a[i];b=all.b[i];p=all.p[i]
								return(p*(1/(1+(x>0)*(a/(b+x-1))*((alp+T)/(alp+tx))^(r+x))))
							},rep(0,length(x))),1,sum))
}
condexp.bgnbd=function(model,t){
	all.r=model@param$r; all.alp=model@param$alpha; all.a=model@param$a; all.b=model@param$b; all.p=model@param$p
	x=model@control$x; tx=model@control$tx; T=model@control$T
	return(apply(vapply(1:length(all.r),function(i){
								r=all.r[i];alp=all.alp[i];a=all.a[i];b=all.b[i];p=all.p[i]
								return(p*((a+b+x-1)/(a-1)*(1-((alp+T)/(alp+T+t))^(r+x)*hyperg_2F1(r+x,b+x,a+b+x-1,t/(alp+T+t)))/(1+(x>0)*a/(b+x-1)*((alp+T)/(alp+tx))^(r+x))))
							},rep(0,length(x))),1,sum))
}
paramplot.bgnbd=function(model,...) {
	par(mfrow=c(1,2))
	plotGamma(model@param$r,model@param$alpha,model@param$p,main='Purchase Process')
	plotBeta(model@param$a,model@param$b,model@param$p,main='Death Process')
	par(mfrow=c(1,1))
}