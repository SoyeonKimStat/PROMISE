# 02/24/16
library(glmnet)

# apply stability selection after selecting parameters (alpha, lambda.min, cutoff) using stab.cv function
# family and standardize option
stability.path <- function(x.train, y.train, family="binomial", alpha = alpha, standardize=FALSE){

  p <- ncol(x.train)
  pre.fit <- glmnet(x.train, y.train, alpha = alpha, family = family, standardize = standardize)
  lambda.range <- pre.fit$lambda
  
  lambda.length = length(lambda.range)
  usage <- matrix(0, lambda.length, p)
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
    fit = glmnet(x.train[which,], y.train[which], alpha = alpha, lambda = lambda.range, family  =family, standardize = FALSE) 
    
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
  
  list(phi=phi,lambda.range=lambda.range)
}
