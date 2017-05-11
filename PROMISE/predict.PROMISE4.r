#02/24/16
# get final auc and nonzero list after subsampling procedure after doing CV
# to get alpha, lambda, cutoffs
# object is from PROMISE function
# binary outcome pred: mode of y.train
# continous outcome pred: mean(y.train)
library(glmnet)
library(ROCR)
source("auc.r")


predict.PROMISE <- function(object, x.train, y.train, x.test, y.test, family="gaussian", standardize=FALSE)  { 
  p <- ncol(x.train)

  phi <- object$phi
  final.phi <- phi[object$lambda.range >=object$lambda.cutoff,]
  
  
  if(sum(final.phi) != 0) {
    if(!is.null(nrow(final.phi))){
      result <- apply(final.phi,2,max)    
      # what percentage of lambda has been used 
    } else { 
      result <- final.phi
    }
    ind <- result > object$cutoff
  } else {
    ind <- rep(FALSE, p)
  }
  
 
  
  coefs <-as.numeric(ind)
  names(coefs) <- colnames(x.train)

  if(sum(ind) == 0) { # no variable is selected
    intercept <- 0
    if(family == "binomial") {
      pred <- rep(as.numeric(names(sort(table(y.train),decreasing=TRUE)[1])), length(y.test))
    } else if(family == "gaussian") {
      pred <- rep(mean(y.train), length(y.test))
    } else {
      stop(paste("family=", family, " is not available option", sep=""))
    }
  } else if(sum(ind) > 1) {
    x.sel <- x.train[, ind]
    if(is.null(nrow(x.sel))) {
      x.sel <- as.matrix(x.sel)
    }
    fit <- glmnet(x.sel, y.train, family =family, alpha = 0, lambda = seq(0, 10, length = 100), standardize = standardize)
    
    #coef(test2, s=0)
    pred <- as.numeric(predict(fit, x.test[,ind], s=0, type ="response"))
    
    pre_coef <- predict(fit, x.test[,ind], s=0, type ="coefficients")
    coefs[ind]  <- as.numeric(pre_coef)[-1]
    intercept = as.numeric(pre_coef)[1]
  } else {
    train = data.frame(v=x.train[,ind],  res=y.train)
    fit <- glm(res~., data=train, family=family)
    pred <- predict(fit, newdata=data.frame(v=x.test[,ind]), type ="response")     
    pre_coef<- fit$coefficients
    coefs[ind] <- pre_coef[-1]
    intercept = as.numeric(pre_coef)[1]
  }
  if(family == "binomial") {
    loss <- auc(pred, y.test)
  } else if(family == "gaussian") {
    loss <- mean((pred-y.test)^2)
  } else {
    stop(paste("family=", family, " is not available option", sep=""))
  }
  names(intercept) = "Intercept"
  
  list(loss=loss, nonzero=abs(sign(coefs)), coefficients=c(intercept, coefs), response=pred)

}

