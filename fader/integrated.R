#The pareto NBD model
control.pnbd=function(model,...){
	return(list(x=model$raw$x,tx=model$raw$tx,T=model$raw$T,mult=1))
}
ll.pnbd=function(model,param=NULL,x=model$control$x[-length(model$control$x)]){
	r=(param[1]); alp=(param[2]); s=(param[3]); bet=(param[4])
	x=model$raw$x; t=model$raw$tx; T=model$raw$T
	maxab = max(alp,bet)
	absab = abs(alp-bet)
	param2 = s+1
	if (alp < bet) param2 = r + x
	part1 = (alp^r*bet^s/gamma(r))*gamma(r+x)
	part2 = 1/((alp+T)^(r+x)*(bet+T)^s)
	if (absab == 0){
		F1 = 1/((maxab+t)^(r+s+x))
		F2 = 1/((maxab+T)^(r+s+x))
	}
	else{
		F1 = hyperg_2F1(r+s+x,param2,r+s+x+1,absab/(maxab+t))/((maxab+t)^(r+s+x))
		F2 = hyperg_2F1(r+s+x,param2,r+s+x+1,absab/(maxab+T))/((maxab+T)^(r+s+x))
	}
	return(log(part1*(part2+(s/(r+s+x))*(F1-F2))))
}
predict.pnbd=function(model,...) rev(cumsum(rev(standardpredict(model,...))))
model.pnbd=function(model,nseg=1) standardmodel(model,c('r','alpha','s','beta','p'),nseg)
mean.pnbd=function(model) return(1/model$param$lambda)
var.pnbd=function(model) return((1-model$param$lambda)/model$param$lambda^2)
residuals.pnbd=function(model) standardresid(model)
print.pnbd=function(model) standardprint(model)
myplot.pnbd=function(model,...) standardplot(model,...)


#The BG/BB model
control.bgbb=function(model,...){
	return(list(x=model$raw$x,tx=model$raw$tx,T=model$raw$T,mult=model$raw$num,num=1))
}
ll.bgbb=function(model,param=NULL,x=model$control$x){
	a = param[1]; b = param[2]; g = param[3]; d = param[4]; 
	T = model$control$T; x = model$control$x; tx = model$control$tx;
	denom_ab = lbeta(a,b); denom_gd = lbeta(g,d);
	lik = exp(lbeta(a+x, b+T-x) - denom_ab + lbeta(g,d+T) - denom_gd);
	count = T - tx - 1;
	lik=lik+sapply(1:length(lik),function(j){
				ifelse(count[j]>=0,sum(sapply(0:count[j],function(i){exp(lbeta(a+x[j],b+tx[j]-x[j]+i) - denom_ab + lbeta(g+1, d+tx[j]+i) - denom_gd)})),0)
			})
	return (log(lik))
}
predict.bgbb=function(model,...) standardpredict(model,...)
model.bgbb=function(model,nseg=1) standardmodel(model,c('alpha','beta','gamma','delta','p'),nseg)
mean.bgbb=function(model) return(1/model$param$lambda)
var.bgbb=function(model) return((1-model$param$lambda)/model$param$lambda^2)
residuals.bgbb=function(model) standardresid(model)
print.bgbb=function(model) standardprint(model)
myplot.bgbb=function(model,...) standardplot(model,...)



# Functions:
# 1) BGBB.convertRFM: takes in data of repeat transactions and returns Recency-Frequency table
# 2) BGBB.indiv.likelihood: takes in model params, recency (tx) and frequency (x), number of time periods (nDays) and returns probability of this occurrence
# 3) BGBB.log.likelihood: takes in model params, RF table and number of periods (nDays) and returns log likelihood of entire data represented in table
# 4) BGBB.pmf: takes in model params, num periods (nDays) and number of transactions (nTx) and returns probability of nTx in nDays periods 
# 5) BGBB.PredictionPMF: takes in model params, num periods (nDays) and returns probability mass function across all possible num transactions (Equation 7)
# 6) BGBB.PredictionTrans: Equation 8 
# 7) BGBB.FuturePMF: Equation 9
# 8) BGBB.PAlive: Equation 11 (added 2/8/11)
# 9) BGBB.CondExpectation: Equation 13
# 10) BGBB.PlotData: plot histogram, expected transactions (within training period)
# 11) BGBB.PlotDataFuture: plot histogram, expected transactions (for validation period)
# 12) BGBB.ReadFile: read CSV file return data within the begin and end columns provided as function arguments
# 13) BGBB.Simulate: takes in model params, num subjects and num time periods (nDays), returns simulation data reflecting these params


# Uses BGBB log likelihood on multiple cohorts.
MultiCohort.log.likelihood = function(par, tableList, nCohorts, nTotal) {
	llsum = 0
	for (jj in 1:nCohorts) {
		llsum = llsum + BGBB.log.likelihood(par, tableList[[ jj ]], nTotal-jj+1)
	}
	return (llsum)
}

#result <- optim(par_init, MultiCohort.log.likelihood, tableList=rf.tables, 
#		nCohorts=kNCohorts, nTotal = nTotal, control=list(fnscale=-1))

# Modifying to use lchoose instead of combn (can't allocate 31c12).
# Change made Feb 11. Checked with WBEZ test, seems ok.   
BGBB.pmf = function(par, nDays, nTx) {
	a = par[1]; b = par[2]; g = par[3]; d = par[4]; n = nDays; x = nTx;
	denom_ab = lbeta(a,b); denom_gd = lbeta(g,d);
	
	# this can handle x = 0, combn cannot
	sum = exp(lchoose(n,x) + (lbeta(a+x, b+n-x) - denom_ab + lbeta(g,d+n) - denom_gd));
	
	if (x < n) {
		for (i in x:(n-1) ) {
			sum = sum + exp(lchoose(i,x)+(lbeta(a+x,b+i-x) - denom_ab + lbeta(g+1, d+i) - denom_gd));		
		}
	}
	return (sum);
}	

#Equation 7 of BGBB paper
BGBB.PredictionPMF = function(par, nDays) {
	# Calculate pmf of BGBB
	pmf = vector("numeric",nDays+1);
	for (i in 0:nDays) {
		pmf[i+1] = BGBB.pmf(par,nDays,i);
	}
	return (pmf);
}

#Equation 8 of BGBB paper
BGBB.PredictionTrans = function(par, nDays) {
	a = par[1]; b = par[2]; g = par[3]; d = par[4]; n = nDays;
	trans = matrix(0,2,nDays);
	
	for (i in 1:n) {
		trans[1,i] = (a/(a+b))*(d/(g-1))*(1 - (gamma(g+d)/gamma(g+d+i))*(gamma(1+d+i)/gamma(1+d)));
		if (i == 1) {
			trans[2,i] = trans[1,i];
		}
		else {
			trans[2,i] = trans[1,i] - trans[1,i-1];
		}
	}
	return (trans);
}

#Equation 9 of BGBB paper, assume nEnd > nDays, no check for this.
# I think it works now - needed protection against x = nInt  
BGBB.FuturePMF = function(par, nDays, nInt) {
	a = par[1]; b = par[2]; g = par[3]; d = par[4]; n = nDays; 
	pmf = vector("numeric",nInt+1);
	denom_ab = lbeta(a,b); denom_gd = lbeta(g,d);
	
	for (x in 0:nInt) {
		sum = 0;
		if (x == 0) {
			sum = 1 - exp(lbeta(g,d+n) - denom_gd);
		}
		
		sum = sum + exp(lchoose(nInt,x) + lbeta(a+x,b+nInt-x) - denom_ab + lbeta(g,d+n+nInt) - denom_gd);
		
		if (x < nInt) {
			
			for (i in x:(nInt-1)) {
				sum = sum + exp(lchoose(i,x) + lbeta(a+x,b+i-x) - denom_ab + lbeta(g+1,d+n+i) - denom_gd);
			}
			
		}
		
		pmf[x+1]=sum;
	}
	return(pmf);
}

# P(Alive) at n+1 (equation 11 of paper) #added function on 2/8/11
BGBB.PAlive = function(par, table, nDays) {
	a = par[1]; b = par[2]; g = par[3]; d = par[4]; n = nDays;
	denom_ab = lbeta(a,b); denom_gd = lbeta(g,d);
	
	nEntries = dim(table)[1];
	prob = vector("numeric",nEntries);
	
	for (j in 1:nEntries) {
		logsum = -log(BGBB.indiv.likelihood(par,table[j,1],table[j,2],n));
		logsum = logsum + lbeta(a+table[j,1],b+n-table[j,1]) - denom_ab + lbeta(g,d+n+1) - denom_gd;
		prob[j] = exp(logsum); 
	}
	return(prob);	
}	

# Conditional Expectations (equation 13 of paper)
BGBB.CondExpectation = function(par, table, nDays, nInt) {
	a = par[1]; b = par[2]; g = par[3]; d = par[4]; n = nDays;
	
	nEntries = dim(table)[1];
	condExp = vector("numeric",nEntries);
	
	for (j in 1:nEntries) {
		logsum = -log(BGBB.indiv.likelihood(par,table[j,1],table[j,2],n));
		logsum = logsum + lbeta(a+table[j,1]+1,b+n-table[j,1]) - lbeta(a,b)+ lgamma(g+d) - lgamma(1+d);
		sum = exp(logsum)*d/(g-1);
		condExp[j] = sum*(gamma(1+d+n)/gamma(g+d+n)-gamma(1+d+n+nInt)/gamma(g+d+n+nInt));
	}
	return(condExp);
}

# function to generate histograms (for only the training data period)
BGBB.PlotData = function(data,nDays,predPMF,predTrans) {
	windows(h=20,w=20)
	par( mfrow=c(2,2) ) # can change if want 1 or 2 plots only
	
	freq.by.donor = apply(data,1,sum);
	bks = seq(-0.5,max(freq.by.donor)+1,by=1);
	
	hist(freq.by.donor,breaks=bks,labels=FALSE,ylim=c(0, 10000)); # can change labels to TRUE   
	lines(0:nDays,predPMF,col="red");
	
	tracking.incr = apply(data,2,sum);
	print(tracking.incr);
	print(predTrans[2,]);
	
	plot(1:nDays,tracking.incr,type="l", lty=1, col=1, lwd=2,
			main="Incremental Tracking Plot", xlab="Year", ylab="Yearly Transactions" );
	lines(1:nDays, predTrans[2,],col="red");
	
	tracking.cumu = cumsum(tracking.incr); 
	plot(1:nDays,tracking.cumu,type="l", lty=1, col=1, lwd=2,
			main="Cumulative Tracking Plot", xlab="Year", ylab="Cumulative Transactions");	
	lines(1:nDays,predTrans[1,],col="red");
}

# function to generate future prediction plots
BGBB.PlotDataFuture = function(data,nDays,nEnd,predPMF,predTrans) {
	windows(h=20,w=20)
	par( mfrow=c(2,2) ) # can change if want 1 or 2 plots only
	
	futuredata = data[,(nDays+1):nEnd]; # cannot handle beyond N
	freq.by.donor = apply(futuredata,1,sum);
	bks = seq(-0.5,max(freq.by.donor)+1,by=1);
	
	hist(freq.by.donor,breaks=bks,labels=TRUE,ylim=c(0, 16000));
	lines(0:(nEnd-nDays),predPMF,col="red");
	
	tracking.incr = apply(data,2,sum);
	print(tracking.incr);
	print(predTrans[2,]);
	
	plot(1:nEnd,tracking.incr,type="l", lty=1, col=1, lwd=2,
			main="Incremental Tracking Plot", xlab="Year", ylab="Yearly Transactions" );
	lines(1:nEnd, predTrans[2,],col="red");
	
	tracking.cumu = cumsum(tracking.incr); 
	plot(1:nEnd,tracking.cumu,type="l", lty=1, col=1, lwd=2,
			main="Cumulative Tracking Plot", xlab="Year", ylab="Cumulative Transactions");	
	lines(1:nEnd,predTrans[1,],col="red");
}

BGBB.ReadFile = function(file, bcol, ecol) {
	data = read.csv(file,header=FALSE);
	nSubjects = dim(data)[1];
	return(data[1:nSubjects,bcol:ecol]);
}

BGBB.simulate = function(par, nSubjects, nDays) {
	a = par[1]; # alpha
	b = par[2]; # beta
	g = par[3]; # gamma
	d = par[4]; # delta
	
	data = matrix(0,nrow = nSubjects, ncol = nDays);
	
	p = rbeta(nSubjects, a, b);
	theta = rbeta(nSubjects, g, d);
	
	for (i in 1:nSubjects) {
		dDay = rgeom(1,theta[i]);
		if (dDay > 0) {
			data[i,1:min(nDays, dDay)] = rbinom(min(nDays, dDay), 1, p[i]);
			# 1 in rbinom is to indicate that each value is from one trial
		}
	}
	
	return (data);
}
