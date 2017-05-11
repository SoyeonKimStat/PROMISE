# 02/22/16
# stabilty path for PROMISE
# given lambda range
# given alpha, grid search (pi and max lambda) that maximize prediction accuracy
# any time of response

# Arguments (input options)
# data: input matrix. the last column is the response
# p.train: what percentage of data will be used for train data (default:0.75)
# option: "lasso"(default) or "elastic net"
# lambda.rule: "lambda.1se" or "lambda.min"
# family: "gussian"(default) or "binomial"


# stratified subsmapling for binary outcome
library(glmnet)



# repi =1000;p.train = c(3/4);family = "binomial"; lambda.rule = "lambda.1se"
# n=200
# nm=250;nnonzero=12;between.cor <- 0.2;beta.coef <- 0.5
# nsim =".four.trt.normal.mycoding.beta"
# source("stratified.divide.data.r")
# source(paste("sim",nsim,".data.r",sep="")) #source("sim1.r")
# within.cor <- 0.8
# data <- sim.data(n) 
# div.data <- s.divide.data(data, p.train)
# x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test;
# alpha <- 1

#x.train=matrix(rnorm(100*20),100,20)
#y.train=rnorm(100)
#family="gaussian";alpha=1;lambda.range = glmnet(x.train, y.train, alpha = alpha, family = family)$lambda


stability.cv.path <- function(x.train, y.train, family="binomial",alpha = alpha, lambda.range, standardize=FALSE){

  

 
  lambda.length = length(lambda.range)
  p <- ncol(x.train)
  usage <- matrix(0, lambda.length,p)
  total <- rep(0, lambda.length)
  n.train <- length(y.train)
  n.half.train <- floor(n.train/2)
  if(family == "binomial") {
    n.y1 <- sum(y.train)
    halfsize.y1 <- floor(n.y1/2)
  }
  for(i in 1:100){# fix it!
    # size n/2 stratified subsampling. 
    # Two sampling procedure
    # First from y=1 and second sampling is from y=0
    if(family == "binomial") {
      which <- c(sample(1:n.y1, halfsize.y1),sample(c(n.y1+1):n.train, c( n.half.train-halfsize.y1))) 
    } else {
      which <- sample(1:n.train,  n.half.train)
    }
    fit = glmnet(x.train[which,], y.train[which], alpha = alpha, lambda = lambda.range, family  =family, standardize = standardize) 
    
    cf <- coef(fit)[-1,]
    m.coef <- t(as(cf, "matrix"))
    lambda.act.length = length(fit$lambda)
    if(lambda.act.length !=lambda.length) {# when glmnet does not converge for a small lambda(lambda_e)
      # it does not gives result corresponds to smaller lambda than lambda_e
      # Assume that for smaller lambda than lambda_e,
      # the result is the same as lambda_e
      pre.use <- abs(sign(m.coef))
      nempty.row <- lambda.length - lambda.act.length 
      last.pre.use <- pre.use [lambda.act.length ,]
    
      fempty.row <- t(matrix(rep(last.pre.use, nempty.row ), ncol = nempty.row ))
      use <- rbind(pre.use, fempty.row)
      
     } else {
      use <- abs(sign(m.coef)) # 100 * p matrix with indicators whether a variable is selected
    }
       # total: total number of selected variables for each lambda over 100 times of run
    usage <- usage + use # how many times each variable are selected
  }
  phi <- usage/100 # selection probability of each variable
  
  phi

}
# 
# test
# repi =1000;p.train = c(3/4);family = "binomial"; lambda.rule = "lambda.1se"
# n=200
# nm=250;nnonzero=6;between.cor <- 0.2
# within.cor <- 0.8
# 
# nsim =".four.trt.normal.mycoding"
# 
# source(paste("sim",nsim,".data.r",sep=""))
# source("divide.data.r")
# data <- sim.data(n) 
# div.data <- divide.data(data,p.train)
# x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test;
# source("stab.cv.r")
# test0 <- stab.cv(x.train,y.train)
#  nfolds cross validation to find alpha, lambda_min, 
