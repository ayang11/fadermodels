#install.packages('gsl')
rm(list=ls())
setwd('C:\\Users\\Andrew\\Desktop\\Workspace\\models\\fader')
source('utilities.R')
source('counting.R')
source('contime.R')
source('distime.R')
source('choice.R')
r=function(str) return(paste('raw',str,sep='\\'))

standardtest('count.csv','pois')
standardtest('count.csv','nbd')
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
mod=c(0.195,0.054,0.060,0.116,0.362,0.139,0.151,0.175)
mean(mod)
dirExp(mod)
dirPen(mod)
dirNo(mod)
dirLoyal(mod)
dirOnce(mod)
dirMkt(mod)
dirScr(mod)


data.bgbb=read.csv(r('bgbb.csv'))
mod=fm(data.bgbb,'bgbb')

data.pnbd=read.csv(r('pnbd.csv'))
mod=fm(data.pnbd,'pnbd')

# TODO Figure out how to graph all the cool stats for the integrated models
# TODO Figure out how to add spikes to each model
# TODO Double check all code
# TODO clean up dirichlet