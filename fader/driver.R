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
source('fmspec.R')
r=function(str) return(paste('raw',str,sep='\\'))


standardtest('count.csv','pois',nseg=2,spike=T)
standardtest('count.csv','pois',nseg=1,spike=T)
standardtest('count.csv','nbd',nseg=1,spike=T)
standardtest('contime.csv','exp',pop=1499)
standardtest('contime.csv','eg',pop=1499)
standardtest('contime.csv','wg',pop=1499)
standardtest('contime.csv','gg',pop=1499)
standardtest('distime.csv','geom',pop=1000)
standardtest('distime.csv','bg',pop=1000)
standardtest('distime.csv','dw',pop=1000)
standardtest('distime.csv','bdw',pop=1000)
standardtest('choice.csv','bb')

a=standardtest('pnbd.csv','pexp')
a=standardtest('pnbd.csv','bgnbd')
a=standardtest('bgbb.csv','bgbb')
cumtracking(a)


data=read.csv(r('pnbd.csv'))
a=fm(data,'pexp',1)
cumtracking(a)
b=fm(data,'pexp',2)
cumtracking(b)

a=fm(data,'bgnbd',1)
cumout(a)
b=fm(data,'bgnbd',2)
cumout(b)

data=read.csv(r('bgbb.csv'))
a=fm(data,'bgbb',1)
cumout(a)
b=fm(data,'bgbb',2)
cumout(b)

standardtest('pnbd.csv','pnbd',2)
standardtest('dirich.csv','dir')


#mod=c(0.195,0.054,0.060,0.116,0.362,0.139,0.151,0.175)
data.dir=read.csv(r('dirich.csv'))
mod=fm(data.dir,'dir',1,tries=10);mod
predict(mod)
mean(mod)
dirExp(mod)
dirPen(mod)
dirNo(mod)
dirLoyal(mod)
dirOnce(mod)
dirMkt(mod)
dirScr(mod)


x=c(rnorm(800,15,2),rnorm(200,-10,4))
mod=fm(x,'norm',nseg=2,tries=50,b=2);mod

# TODO tracking plot for pnbd is still different
# TODO clean up dirichlet, and check it and integrated models if they work for multiple segments


a=fm(data,'pnbd')
b=fm(data,'pnbd',2)
post(b)