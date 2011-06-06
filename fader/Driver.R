# TODO: Add comment
# 
# Author: Andrew
###############################################################################


setwd('C:\\Users\\Andrew\\Desktop\\Workspace\\models\\fader')
source('utilities.R')
source('counting.R')
source('timing.R')
source('choice.R')

data.pois=read.csv('count.csv')
mod=fm(data.pois,'pois',6,t=1);rmse(mod);mod;myplot(mod,t=1);mean(mod);var(mod)

data.nbd=read.csv('count.csv')
mod=fm(data.nbd,'nbd',2);rmse(mod);mod;myplot(mod,t=1);mean(mod);var(mod)
predict(mod,x=c(30:50),t=2)

data.eg=read.csv('kb.csv',header=FALSE)
mod=fm(data.eg,'eg',1,num=1499);rmse(mod);mod;myplot(mod);mean(mod);var(mod)

data.wg=read.csv('kb.csv',header=FALSE)
mod=fm(data.wg,'wg',1,num=1499);rmse(mod);mod;myplot(mod);mean(mod);var(mod)

data.gg=read.csv('kb.csv',header=FALSE)
mod=fm(data.wg,'gg',1,num=1499);rmse(mod);mod;myplot(mod);mean(mod);var(mod)

data.geom=read.csv('time.csv')
mod=fm(data.geom,'geom',1);rmse(mod);mod;myplot(mod);mean(mod);var(mod)

#Fix this
data.bg=read.csv('time.csv')
mod=fm(data.bg,'bg',1);rmse(mod);mod;myplot(mod);mean(mod);var(mod)

data.dw=read.csv('time.csv')
mod=fm(data.dw,'dw',2);rmse(mod);mod;myplot(mod);mean(mod);var(mod)

data.bdw=read.csv('time.csv')
mod=fm(data.bdw,'bdw',1);rmse(mod);mod;myplot(mod);mean(mod);var(mod)

data.bb=read.csv('bb.csv')
mod=fm(data.bb,nseg=1);rmse(mod);mod;myplot(mod);mean(mod);var(mod)

data.dir=read.csv('dirich.csv')
mod=c(0.195,0.054,0.060,0.116,0.362,0.139,0.151,0.175)
dirPen(mod,data.dir)
dirFreq(mod,data.dir)
dirOnce(mod,data.dir)
dirLoyal(mod,data.dir)
dirMkt(mod,data.dir)
dirScr(mod,data.dir)