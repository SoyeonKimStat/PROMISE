# 02/25/16
# Combine two functions
source("stab.cv.glmnet.faster_fixed_continous2.r")
source("final.stability.path_any.response_entire.lambda.r")
#family="gaussian"; nfolds=5; option = "lasso"; nlambda =20; standardize=FALSE; lambda.rule = "lambda.1se"
PROMISE <- function(x.train, y.train, family="gaussian", nfolds=5, option = "lasso", nlambda =20, standardize=FALSE, lambda.rule = "lambda.1se") {
  best <- stab.cv(x.train, y.train, family=family, nfolds=nfolds, option = option, nlambda=nlambda)
  
  if(lambda.rule == "lambda.1se") {
    best.com <- best$combi.1se
  } else {
    best.com <- best$combi.min
  }
  lambda.cutoff <- best.com$lambda
  alpha <- best.com$alpha
  cutoff <- best.com$cutoff
  stab.path <- stability.path(x.train, y.train, family=family, alpha = alpha)
  list(alpha=alpha, lambda.cutoff=lambda.cutoff, cutoff=cutoff, phi=stab.path$phi, lambda.range=stab.path$lambda.range)
}

