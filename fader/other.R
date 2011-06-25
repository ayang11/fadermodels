# TODO: Add comment
# 
# Author: Andrew
###############################################################################

#The poisson model is a discrete one parameter counting model. It answers how many
control.norm=function(model,...){
	return(list(x=model$raw,y=model$raw,names=c('mean','sigma','p'),...))
}
ll.norm=function(model,param=NULL,x=model$control$x){
	return(log(dnorm(x,param[1],param[2])))
}
model.norm=function(model) model.fm(model,pos=2)
mean.norm=function(model) return(sum(model$param$mean*model$param$p))
var.norm=function(model) return(model$param$lambda)