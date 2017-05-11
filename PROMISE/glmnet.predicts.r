# glmnet prediction using ridge regression without penalty 
source("auc.r")
glmnet.predicts <- function(nonzero.ind, x.sub, y.sub, x.val, y.val, family) {
  x.sel <- x.sub[, nonzero.ind]
  if(sum(nonzero.ind) >1) { # if two or more variable selected use glmnet
    #       if(is.null(nrow(x.sel))) {
    #         x.sel <- as.matrix(x.sel)
    #       } 
    test <- glmnet(x.sel, y.sub, family =family, alpha = 0, lambda = seq(0, 10, length = 100), standardize = FALSE)
    pred <- predict(test, x.val[,nonzero.ind], s=0, type ="response")
    if(family == "binomial") {
      loss <- auc(pred, y.val)
    } else if(family == "gaussian") {
      loss <- mean((pred-y.val)^2)
    } else {
      stop(paste("family=", family, " is not available option", sep=""))
    }
  } else if(sum(nonzero.ind) == 0) {
    if(family == "binomial") {
      loss <- 0.5 
    } else if(family == "gaussian") {
      loss <- mean((mean(y.val)-y.val)^2)
    } else {
      stop(paste("family=", family, " is not available option", sep=""))
    }
  
  } else {   # one variable selected: use glm
    train = data.frame(v=x.sub[,nonzero.ind],  res=y.sub)
    fit <- glm(res~., data=train, family=family)
    pred <- predict(fit, newdata=data.frame(v=x.val[,nonzero.ind]), type ="response")     
    if(family == "binomial") {
      loss <- auc(pred, y.val)
    } else if(family == "gaussian") {
      loss <- mean((pred-y.val)^2)
    } else {
      stop(paste("family=", family, " is not available option", sep=""))
    }
  }
  loss 
}