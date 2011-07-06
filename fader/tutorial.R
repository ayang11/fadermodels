#Installation
setwd('<INSERT DIRECTORY>')
install.packages('fadermodels_1.0.zip',repos=NULL)
library('fadermodels')
install.packages('gsl')
library('gsl')

#Introduction
data(data.count,package='fadermodels')
mod=fm(data.count,'pois')
mod
mean(mod) #Mean of the Poisson fit
vcov(mod) #Variance of the Poisson fit
predict(mod) #Predicts the number of people at each occurance 
resid(mod) #finds the residuals of the prediction
sse(mod) #finds the sse of the prediction
rmse(mod) #finds the rmse of the prediction
barplot(mod)
chitest(mod)

#Latent Class models
mod2=fm(data.count,'pois',nseg=2)
mod3=fm(data.count,'pois',nseg=3)
mod4=fm(data.count,'pois',nseg=4,tries=5)
mod5=fm(data.count,'pois',nseg=5,tries=5)
mod6=fm(data.count,'pois',nseg=6,tries=5)
lltest(mod,mod2)
lltest(mod2,mod3)
lltest(mod3,mod4)
lltest(mod4,mod5)
lltest(mod5,mod6)
best=bictest(mod,mod2,mod3,mod4,mod5,mod6)
barplot(best)
chitest(best)
rmse(best) 

#Spikes and posterior
mod2s=fm(data.count,'pois',nseg=2,spike=T)
lltest(mod2,mod2s)
barplot(mod2s)
post(mod2s)

#Mixing distribution
mixmod=fm(data.count,'nbd')
mixmod
barplot(mixmod)
mean(mixmod)
chitest(mixmod)

#Analyzing the models
mixmod2=fm(data.count,'nbd',nseg=2,spike=T)
mixmod2
best=bictest(mod,mod2,mod3,mod4,mod5,mod6,mod2s,mixmod,mixmod2)
best
paramplot(mixmod)
paramplot(mod2s)

#Changing parameters
newmod=update(mixmod,r=1,alpha=0.1)
newmod
barplot(newmod)
paramplot(newmod)

#Predictions through time
predict(best,t=2)
predict(best,t=1/2)
predict(best,t=3,x=0:50)


#Dirichlet
data(data.dir,package='fadermodels')
dir=fm(data.dir,'dir',tries=10)
barplot(dir)
mean(dir)
pen.dir(dir) #Penetration level
loyal.dir(dir) #Probability that a customer will be 100% loyal in your segment
once.dir(dir) #Probability that a customer will be a one time customer in your segment
mkt.dir(dir) #Market share of each segment
scr.dir(dir) #Share of category requirements
dir=update(dir,x1=.195,x2=.054,x3=0.06,x4=0.116,x5=.362,x6=.139,x7=.151,x8=.175)

#Integrated Models

data(data.conint,package='fadermodels')
data(data.disint,package='fadermodels')
a=fm(data.conint,'pexp',1)
b=fm(data.conint,'pexp',2,tries=2)
c=fm(data.conint,'bgnbd',1)
d=fm(data.conint,'bgnbd',2,tries=2)
e=fm(data.conint,'pnbd',1)
f=fm(data.conint,'pnbd',2,tries=2)
g=fm(data.disint,'bgbb',1)
h=fm(data.disint,'bgbb',2,tries=2)

par(mfrow=c(2,4))
c1=cumtracking(a)
c2=cumtracking(b)
c3=cumtracking(c)
c4=cumtracking(d)
c5=cumtracking(e)
c6=cumtracking(f)
c7=cumtracking(g)
c8=cumtracking(h)

#cumulative tracking plots
plot(c1)
plot(c2)
plot(c3)
plot(c4)
plot(c5)
plot(c6)
plot(c7)
plot(c8)

#weekly tracking plots
plot(diff(c1))
plot(diff(c2))
plot(diff(c3))
plot(diff(c4))
plot(diff(c5))
plot(diff(c6))
plot(diff(c7))
plot(diff(c8))

#Model and actual purchases
barplot(a)
barplot(b)
barplot(c)
barplot(d)
barplot(e)
barplot(f)
barplot(g)
barplot(h)


#Discrete time Models
data(data.distime,package='fadermodels')
geommod=fm(data.distime,'geom',pop=1000)
geommod2=fm(data.distime,'geom',nseg=2,pop=1000)
bgmod=fm(data.distime,'bg',pop=1000)
dwmod=fm(data.distime,'dw',pop=1000)
bdwmod=fm(data.distime,'bdw',pop=1000)
post(geommod2)
best=bictest(geommod,geommod2,bgmod,dwmod,bdwmod)
best
chitest(best)
barplot(best)
rmse(best)
paramplot(best)
predict(geommod,x=c(0,3:15))


#Continuous time models
data(data.contime,package='fadermodels')
expmod=fm(data.contime,'exp',pop=1499)
egmod=fm(data.contime,'eg',pop=1499)
wgmod=fm(data.contime,'wg',pop=1499)
ggmod=fm(data.contime,'gg',pop=1499)
best=bictest(expmod,egmod,wgmod,ggmod)
best
data.frame(data.contime,mod=predict(best)[-length(predict(best))])
paramplot(best)
paramplot(best,bounds=c(0,1),main='Distribution of Mixing')
predict(egmod,x=c(0,3:15))


#Models for choice
data(data.choice,package='fadermodels')
data(data.choice2,package='fadermodels')
bbmod=fm(data.choice,'bb')
bbmods=fm(data.choice,'bb',spike=T)
bbmod2=fm(data.choice,'bb',nseg=2)
lltest(bbmod,bbmods)
lltest(bbmod,bbmod2)
best=bictest(bbmod,bbmods,bbmod2)
paramplot(bbmod2)
