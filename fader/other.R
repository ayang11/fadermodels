# TODO: Add comment
# 
# Author: Andrew
###############################################################################

#The poisson model is a discrete one parameter counting model. It answers how many
control.norm=function(model,...){
	return(list(x=model$raw,mult=1,...))
}
ll.norm=function(model,param=NULL,x=model$control$x){
	return(log(dnorm(x,param[1],param[2])))
}
predict.norm=function(model,...) standardpredict(model,...)
model.norm=function(model) standardmodel(model,c('mean','sigma','p'),pos=2)
mean.norm=function(model) return(sum(model$param$mean*model$param$p))
var.norm=function(model) return(model$param$lambda*t)
residuals.norm=function(model) predict(model)-model$control$x
print.norm=function(model) standardprint(model)
myplot.norm=function(model,...) standardplot(model,...)