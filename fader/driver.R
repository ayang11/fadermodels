#install.packages('gsl')
rm(list=ls())
setwd('C:\\Users\\Andrew\\Desktop\\Workspace\\models\\fader')
source('utilities.R')
source('counting.R')
source('contime.R')
source('distime.R')
source('integrated.R')
source('other.R')
source('choice.R')
r=function(str) return(paste('raw',str,sep='\\'))

standardtest('count.csv','pois')
standardtest('count.csv','nbd',nseg=2)
standardtest('contime.csv','eg',num=1499)
standardtest('contime.csv','wg',num=1499)
standardtest('contime.csv','gg',num=1499)
standardtest('distime.csv','geom')
standardtest('distime.csv','bg')
standardtest('distime.csv','dw')
standardtest('distime.csv','bdw')
standardtest('choice.csv','bb')

data.dir=read.csv(r('dirich.csv'))
mod=fm(data.dir,'dir',1);mod
predict(mod)
#mod=c(0.195,0.054,0.060,0.116,0.362,0.139,0.151,0.175)
mean(mod)
dirExp(mod)
dirPen(mod)
dirNo(mod)
dirLoyal(mod)
dirOnce(mod)
dirMkt(mod)
dirScr(mod)


data.bgbb=read.csv(r('bgbb.csv'))
mod=fm(data.bgbb,'bgbb',2,tries=4);mod;rmse(mod);mean(mod)
myplot(mod)

data.pnbd=read.csv(r('pnbd.csv'))
mod=fm(data.pnbd,'pnbd',2);mod
myplot(mod)


data.bgnbd=data.pnbd
mod=fm(data.bgnbd,'bgnbd',2);mod
myplot(mod)

data.pexp=data.pnbd
mod=fm(data.pexp,'pexp',4);mod;rmse(mod)
mean(mod)
myplot(mod)

x=c(rnorm(80,15,2),rnorm(20,-10,4))
mod=fm(x,'norm',nseg=2,tries=20,b=2);mod

# TODO check need tracking plots
# TODO Figure out how to add spikes to each model
# TODO Double check all code
# TODO clean up dirichlet