#install.packages('gsl')
rm(list=ls())
setwd('Y:\\Portable Apps\\Workspace\\models\\fader')
source('utilities.R')
source('counting.R')
source('contime.R')
source('distime.R')
source('integrated.R')
source('other.R')
source('choice.R')
r=function(str) return(paste('raw',str,sep='\\'))
#pois LL = -929.04, lambda = 4.456
#pois + spike LL = -831, lambda = 5.492, p=0.811
#NBD LL = -649.7, r = 0.969, alpha = 0.218.
#exp with spike: p=0.085, theta=0.066
#EG LL = -681.4, r = 0.050, alpha = 7.973.
#WG LL = -681, r = 0.031, alpha = 6.199, c=1.241.
#geom LL = -1637.09, theta = 0.103.
#BG LL = -1611.16, alpha = 0.668, beta = 3.806.
#dw, theta = 0.15, c=0.77, regular: theta=0.38, c=0.57
#bdw, alpha=0.21, beta=1.43, c=1.72. Regular: alpha=0.46, beta=0.78, c=1.28
#Ben's Knick Knack
#BB LL=-200.5, alpha=0.439, beta=95.411, p=0.0046
#BB alpha = 0.487, beta= 0.826,LL= -35,516.1

#standardtest('count.csv','pois',nseg=1,spike=T)
#standardtest('count.csv','nbd',nseg=1)
#standardtest('contime2.csv','eg',pop=1499)
#standardtest('contime2.csv','wg',pop=1499)
#standardtest('contime2.csv','gg',pop=1499)
#standardtest('distime2.csv','geom',pop=1000)
#standardtest('distime2.csv','bg',pop=1000)
#standardtest('distime2.csv','dw',pop=1000)
#standardtest('distime2.csv','bdw',pop=1000)
#standardtest('choice.csv','bb')


#LL:	-508.0080
data.dir=read.csv(r('dirich.csv'))
mod=fm(data.dir,'dir',1,tries=10);mod
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


#BG/BB alpha = 1.204, beta =0.750,gamma= 0.657, delta= 2.783, LL= -33,225.6
data.bgbb=read.csv(r('bgbb.csv'))
mod=fm(data.bgbb,'bgbb',1);mod;rmse(mod);mean(mod)
myplot(mod)
cumtracking(mod)

#BG/NBD     Pareto/NBD
#r 0.243     0.553
#alpha 4.414 10.578
#a 0.793
#b 2.426
#s           0.606
#beta        11.669
#LL -9582.4 -9595.0
data.pnbd=read.csv(r('pnbd.csv'))
mod.pnbd=fm(data.pnbd,'pnbd',1);mod.pnbd;rmse(mod.pnbd);mean(mod.pnbd)
myplot(mod.pnbd)
cumtracking(mod.pnbd)

data.bgnbd=data.pnbd
mod=fm(data.bgnbd,'bgnbd',1);mod;rmse(mod);mean(mod)
myplot(mod)
cumtracking(mod)

data.pexp=data.pnbd
mod=fm(data.pexp,'pexp',2);mod;rmse(mod);mean(mod)
mean(mod)
myplot(mod)
growth(cumtracking(mod))

x=c(rnorm(800,15,2),rnorm(200,-10,4))
mod=fm(x,'norm',nseg=2,tries=50,b=2);mod

# TODO tracking plot for pnbd is still wrong. 
# TODO Maybe an exponential model
# TODO clean up dirichlet, and check it