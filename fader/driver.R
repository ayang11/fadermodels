source('fmspec.R')
source('utilities.R')
source('counting.R')
source('contime.R')
source('distime.R')
source('integrated.R')
source('other.R')
source('choice.R')
source('dirmethods.R')

data.conint=read.csv('raw//data.conint.csv')
data.disint=read.csv('raw//data.disint.csv')
data.count=read.csv('raw//data.count.csv')
data.contime=read.csv('raw//data.contime.csv')
data.distime=read.csv('raw//data.distime.csv')
data.choice=read.csv('raw//data.choice.csv')
data.choice2=read.csv('raw//data.choice2.csv')
data.choice3=read.csv('raw//data.choice3.csv')
data.dir=read.csv('raw//data.dir.csv')

# TODO tracking plot for pnbd is still different 
# TODO some means/variances are still off
# TODO write instructions
