#install instructions
#PATH C:\Program Files (x86)\GnuWin32\bin;%path%
#PATH C:\Program Files\R\R-2.12.2\bin;%path%
#R CMD check name
#R CMD build name
#R CMD INSTALL --build name.tar.gz

rm(list=ls())

setwd('C:\\Users\\Andrew\\Desktop\\Workspace\\models\\fader')
#files=c(paste(sep='\\','C:\\Users\\Andrew\\Desktop\\Workspace\\models\\fader',c('fmspec.R','utilities.R','counting.R','contime.R','distime.R','integrated.R','other.R','choice.R','dirmethods.R')))
package.skeleton('fadermodels',code_files='C:\\Users\\Andrew\\Desktop\\Workspace\\models\\fader\\aggregate.R')


source('fmspec.R')
source('utilities.R')
source('counting.R')
source('contime.R')
source('distime.R')
source('integrated.R')
source('other.R')
source('choice.R')
source('dirmethods.R')

# TODO tracking plot for pnbd is still different 
# TODO some means/variances are still off
# TODO write instructions
