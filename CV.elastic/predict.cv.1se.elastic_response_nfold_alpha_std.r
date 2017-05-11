# 07/17/13
# make prediction using cross validation and elastic-net/lasso
# using 5folds CV, lambda.min rule and glmnet package
# data: predictors + response (the last column)
# option: either "lasso" or "elastic"
# family: "gussian" or "binomial"
# value (result)
# loss:mean squared error for gaussian family and auc for binomial family
# varcount: indicator vector length p; 1 - selected and 0 - not selected
## nsim = 1; source(paste("sim",nsim,".data.r",sep="")) ;data <- sim.data(240) #
#p.train <- 1/6; family = "gaussian";option ="lasso"
# lasso var_sel is changed
# no stratified sampling for train vs test set
library(glmnet)
library(ROCR)
source("auc.r")
source("divide.data.r")
source("var.count.r")
source("cv.1se.elastic3.r")
########################################################


# Variable selction and mse using lambda and alpha gotten from C.V. using test data
# lambda.rule either "lambda.min" or "lambda.1se"
# lambda.rule = "lambda.min"
predict.cv.elastic <- function(data, p.train=c(3/4), alpha= NA ,lambda.rule="lambda.1se",family="binomial", nfolds=5, standardize = TRUE) {
  div.data <- divide.data(data, p.train)
  x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test;  
  if(family == "gaussian") {
    type1 <- "link"   
    type.measure = "mse"
  } else {
    type1 <- "response"
    type.measure = "auc"
  }
  if(is.na(alpha)) {
    min.four <- cv.elastic(x.train, y.train, family,nfolds = nfolds,standardize=standardize)
    if(lambda.rule == "lambda.min") {
      min.comb <- min.four[1:2]
    } else {
      min.comb <- min.four[3:4]
    }
    cvtest <- glmnet(x.train, y.train,  family=family,alpha = min.comb[1], standardize = standardize)
    pred <- predict(cvtest, newx=x.test, s = min.comb[2], type = type1)
    var_sel <- predict(cvtest, newx=x.test,s = c(min.comb[2], 0.01), type = "nonzero")[[1]]
  } else { # option = "lasso"
    cvtest = cv.glmnet(x.train, y.train, alpha=alpha, family=family,type.measure=type.measure,  nfolds = nfolds, standardize = standardize)
    final.lambda <- get(lambda.rule,envir=list2env(cvtest))
    if(!is.na(final.lambda)) { 
      pred <- predict(cvtest, newx=x.test, s = lambda.rule, type = "response")
      var_sel <- predict(cvtest, newx=x.test, s = c(final.lambda, 0.01),type = "nonzero")[[1]]  
      
       # fixed June 28/13
    } else {
      if(family == "gaussian"){
        pred <- rep(mean(y.train), lenth(y.test))
      } else {
        pred <- rep(as.numeric(names(sort(table(y.train),decreasing=TRUE)[1])), length(y.test))
      }
      var_sel <- NULL 
    }  
  } 
   #  loss function for precition
  if(family == "binomial") { 
    loss <- auc(pred, y.test) # auc
  } else {
    loss <- mean((pred-y.test)^2) # mse
  }
  # variable selection
  varcount <- var_count(var_sel, ncol(x.train))

  list(loss=loss, nonzero=varcount, response=pred, y.test=y.test)
}

  

