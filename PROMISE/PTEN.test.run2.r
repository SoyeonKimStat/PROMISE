source("divide.data.r")
source("PROMISE2.r")
source("predict.PROMISE4.r")
name <- "PTEN"
fileName=paste("transEffectParamOptTest_4Indiv4", name, ".txt", sep="")

data=read.delim(fileName)[,-1]
set.seed(1234)
div.data <- divide.data(data, c(3/4))
x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test;


option="lasso"
object <- PROMISE(x.train, y.train, option=option)
promise.result <- predict.PROMISE(object, x.train, y.train, x.test, y.test)

loss <- promise.result$loss
R.squared <- cor(promise.result$response,y.test)**2
nonzero <- promise.result$nonzero

result=list(loss=loss,R.squared=R.squared,nonzero=nonzero)

set.seed(1234)
div.data <- divide.data(data, c(3/4))
x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test;

option="elastic"
object <- PROMISE(x.train, y.train, option=option)
promise.result <- predict.PROMISE(object, x.train, y.train, x.test, y.test)

loss <- promise.result$loss
R.squared <- cor(promise.result$response,y.test)**2
nonzero <- promise.result$nonzero

result=list(loss=loss,R.squared=R.squared,nonzero=nonzero)

