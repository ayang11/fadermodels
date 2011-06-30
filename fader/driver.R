#install.packages('gsl')
rm(list=ls())
setwd('C:\\Users\\Andrew\\Desktop\\Workspace\\models\\fader')
source('fmspec.R')
source('utilities.R')
source('counting.R')
source('contime.R')
source('distime.R')
source('integrated.R')
source('other.R')
source('choice.R')
source('dirmethods.R')
r=function(str) return(paste('raw',str,sep='\\'))

# TODO tracking plot for pnbd is still different 
# TODO write instructions

data.conint=read.csv(r('pnbd.csv'))
data.disint=read.csv(r('bgbb.csv'))
data.count=read.csv(r('count.csv'))
data.contime=read.csv(r('contime.csv'))
data.distime=read.csv(r('distime.csv'))
data.choice=read.csv(r('choice.csv'))
data.dir=read.csv(r('dirich.csv'))
