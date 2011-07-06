#Construct package

#install instructions
#PATH C:\Program Files (x86)\GnuWin32\bin;%path%
#PATH C:\Program Files\R\R-2.12.2\bin;%path%
#cd C:\Users\Andrew\Desktop\Workspace\models\fader
#R CMD check fadermodels
#R CMD build fadermodels
#R CMD INSTALL --build fadermodels.tar.gz
rm(list=ls())
prompt(data.conint,file='data.conint.Rd')
prompt(data.disint,file='data.disint.Rd')
prompt(data.count,file='data.count.Rd')
prompt(data.contime,file='data.contime.Rd')
prompt(data.distime,file='data.distime.Rd')
prompt(data.choice2,file='data.choice2.Rd')
prompt(data.dir,file='data.dir.Rd')


setwd('C:\\Users\\Andrew\\Desktop\\Workspace\\models\\fader')
files=c(paste(sep='\\','C:\\Users\\Andrew\\Desktop\\Workspace\\models\\fader',c('fmspec.R','utilities.R','counting.R','contime.R','distime.R','integrated.R','choice.R','dirmethods.R')))
package.skeleton('fadermodels',code_files=files)

install.packages('fadermodels_1.0.zip',repos=NULL)
library('fadermodels')


library(gsl)

data.conint=read.csv('raw//data.conint.csv')
data.disint=read.csv('raw//data.disint.csv')
data.count=read.csv('raw//data.count.csv')
data.contime=read.csv('raw//data.contime.csv')
data.distime=read.csv('raw//data.distime.csv')
data.choice=read.csv('raw//data.choice.csv')
data.choice2=read.csv('raw//data.choice2.csv')
data.choice3=read.csv('raw//data.choice3.csv')
data.dir=read.csv('raw//data.dir.csv')

standardtest=function(data,model,nseg=1,...){
	print('Model')
	mod=fm(data,model,nseg,...)
	print(mod)
	print('Prediction/RMSE')
	print(rmse(mod))
	barplot(mod)
	print('Mean/Var')
	print(mean(mod))
	print(vcov(mod))
	print(chitest(mod))
	Sys.sleep(1)
	print(post(mod))
	paramplot(mod)
	return(mod)
}

runallTests=function(){
	standardtest(data.count,'pois',nseg=1)
	standardtest(data.count,'nbd',nseg=1,spike=T)
	standardtest(data.contime,'exp',pop=1499)
	standardtest(data.contime,'eg',pop=1499)
	standardtest(data.contime,'wg',pop=1499)
	standardtest(data.contime,'gg',pop=1499)
	standardtest(data.distime,'geom',pop=1000)
	standardtest(data.distime,'bg',pop=1000)
	standardtest(data.distime,'dw',pop=1000)
	standardtest(data.distime,'bdw',pop=1000)
	standardtest(data.choice,'bb')
	standardtest(data.dir,'dir')
	
	standardtest(data.conint,'pexp')
	standardtest(data.conint,'pnbd')
	standardtest(data.disint,'bgbb')
	standardtest(data.conint,'bgnbd')
	
	a=fm(data.conint,'pexp',1)
	b=fm(data.conint,'pexp',2)
	c=fm(data.conint,'bgnbd',1)
	d=fm(data.conint,'bgnbd',2)
	e=fm(data.conint,'pnbd',1)
	f=fm(data.conint,'pnbd',2)
	g=fm(data.disint,'bgbb',1)
	h=fm(data.disint,'bgbb',2,tries=10)
	
	par(mfrow=c(2,4))
	c1=cumtracking(a)
	c2=cumtracking(b)
	c3=cumtracking(c)
	c4=cumtracking(d)
	c5=cumtracking(e)
	c6=cumtracking(f)
	c7=cumtracking(g)
	c8=cumtracking(h)
	plot(c1)
	plot(c2)
	plot(c3)
	plot(c4)
	plot(c5)
	plot(c6)
	plot(c7)
	plot(c8)
	plot(diff(c1))
	plot(diff(c2))
	plot(diff(c3))
	plot(diff(c4))
	plot(diff(c5))
	plot(diff(c6))
	plot(diff(c7))
	plot(diff(c8))
	
}


dir=standardtest('dirich.csv','dir',2)
dir=update(dir,x1=.195,x2=.054,x3=0.06,x4=0.116,x5=.362,x6=.139,x7=.151,x8=.175)
predict(dir)
mean(dir)
pen.dir(dir)
loyal.dir(dir)
once.dir(dir)
mkt.dir(dir)
scr.dir(dir)
post(dir)




