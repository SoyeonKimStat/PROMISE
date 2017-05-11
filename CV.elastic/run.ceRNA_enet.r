

repi=100;lambda.rule="lambda.min";family="gaussian";p.train=c(3/4);alpha=NA
source("predict.cv.1se.elastic_response_nfold_alpha_std.r")
data.name <- "transEffectParamOptTest_ver2"
data=read.delim(paste(data.name,".txt",sep=""))[,-1]

p=ncol(data)-1
loss <- rep(0, repi)
corr <- rep(0, repi)

nonzero = matrix(data = NA, nrow = repi, ncol = p)

set.seed(1234)
for(i in 1:repi) {
  print(paste("Running", i, "repeatition time"))  
  result.temp <- predict.cv.elastic(data, p.train, alpha = alpha,lambda.rule,family=family,nfolds=10,  standardize = TRUE)

    
  corr[i] <- cor(result.temp$response,result.temp$y.test)
  loss[i] <- result.temp$loss
  nonzero[i,] <- result.temp$nonzero
  
}
colnames(nonzero) <- colnames(data[,1:p])
result=list(mse = loss, nonzero=nonzero, corr=corr) 


save(result,file=paste(data.name, ".CV.alpha=", alpha, ".lambda.rule=",lambda.rule, ".std.RData", sep="")) #data.method.Rdata


