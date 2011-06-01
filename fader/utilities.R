setwd('C:\\Users\\Andrew\\Desktop\\Workspace\\models\\fader')
SSE=function(res){
	act=res$act
	model=res$model
	return(sum((act-model)^2))
}
RMSE=function(res){
	return(sqrt(SSE(res)/nrow(res)))
}
plotModel=function(res){
	len=nrow(res)-1
	act=res$act
	model=res$model
	plot(0:len,act,type='l',col='blue')
	lines(0:len,model,col='red')
}
plotCountModel=function(res){
	len=nrow(res)-1
	act=res$num
	model=res$model
	plot(res$count,act,type='l',col='blue')
	lines(res$count,model,col='red')
}
cov=function(x,y) UseMethod("cov")
var=function(x) UseMethod('var')
ll=function(x) UseMethod('ll')
model=function(x) UseMethod('model')
#var,cov,mean,sd,resid,plot,predict,print

#bg,wg,geom go with time.csv
#pois,nbd go with count.csv
