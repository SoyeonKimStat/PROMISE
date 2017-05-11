# 09/12/16
# add standardized =TRUE option
# Based on http://stats.stackexchange.com/questions/17609/cross-validation-with-two-parameters-elastic-net-case
# changed step 2,3, 4
# Q. which one is more parcimoneous? larger alpha or larger lambda
# No one knows unless you apply it into the model.
# calculate lambda.min and alpha.min first
# find lambda.1se at alpha.min
# famly = "gaussian" or "binomial"
# procedure is described in simulation(a) document
# changed 1se standard

# source("stratified.divide.data.r")
# data <-  read.csv("sub100_array_interaction.csv", header=TRUE)
# 
# 
# repi =50;p.train = c(3/4);family = "binomial"; lambda.rule = "lambda.1se";nfolds=5
# option = "elastic"
# div.data <- s.divide.data(data, p.train)
# x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test;
# test <- cv.elastic(x.train, y.train, family)
# 

library(glmnet)
cv.elastic <- function(x.train, y.train, family, nfolds = 5, standardize=TRUE) {
  seq.alpha <- seq(0.1, 1, by=0.1)  # alpha = 0 is excluded for variable selection
  if (family == "gaussian") {
    type.measure = "mse"
  } else if(family == "binomial") {
    type.measure = "auc"
  } else {
    stop(paste("family=", family, " is not available option", sep="")) 
  }
  # step 1: do all crossvalidations for each alpha
  cvs <- lapply(seq.alpha, function(one.alpha){
    cv.glmnet(x.train, y.train, alpha=one.alpha, family=family, nfolds = nfolds, type.measure = type.measure, standardize = standardize)
  })
  # step 2: collect the optimum lambda for each alpha 
  best<-sapply(seq_along(seq.alpha), function(curi){
    curcvs<-cvs[[curi]]
    one.alpha<-seq.alpha[curi] 
    min.id<-match(curcvs$lambda.min, curcvs$lambda) # index number 
    names(curcvs$nzero)<- NULL
    c(lam=curcvs$lambda.min, alph=one.alpha, cvm=curcvs$cvm[min.id], cvsd=curcvs$cvsd[min.id], nonzero = curcvs$nzero[min.id]) # upper curve
  })
  #step 3: find the overall optimum
  if(family == "gaussian") {
    best.id<-which.min(best["cvm",]) 
    se.min <- best["cvm",best.id]  + best["cvsd",best.id] 
  } else if(family == "binomial") {
    best.id <- which.max(best["cvm",])
    se.min <- best["cvm",best.id]  - best["cvsd",best.id] 
  } else {
    stop(paste("family=", family, " is not available option", sep=""))
  }
  overall.lambda.min<-best["lam",best.id]
  overall.alpha.min<-best["alph",best.id]
  nzero.min <- best["nonzero",best.id] #nonzero of best model
  cvm.min <- best["cvm", best.id]
  names(overall.alpha.min) <- NULL
  names(overall.lambda.min) <- NULL
  names(nzero.min) <- NULL
  names(cvm.min) <- NULL
  best.combi.min <- data.frame(alpha = overall.alpha.min, lambda=overall.lambda.min, df=nzero.min, cvm=cvm.min )
  #step 4: For each alpha, which lambda is the best within 1SE difference from the best model
  set.1se<-lapply(seq_along(seq.alpha), function(curi){
    curcvs<-cvs[[curi]]
    if(family == "gaussian") {
      idmin = curcvs$cvm <= se.min & curcvs$nzero <= nzero.min
    } else {
      idmin = curcvs$cvm >= se.min & curcvs$nzero <= nzero.min 
    }
    if(any(idmin)) { 
      lam1se = max(curcvs$lambda[idmin], na.rm = TRUE)
      min.id <- match(lam1se, curcvs$lambda)
      nzero1se <- curcvs$nzero[min.id]
      cvm1se <- curcvs$cvm[min.id]
    } else {
      lam1se <- NA
      min.id <- NA # index number 
      nzero1se <- NA
      cvm1se <- NA
    }  
    c( alpha=seq.alpha[curi], lambda=lam1se,df = nzero1se,cvm=cvm1se)
  })
  dfr.1se <- as.data.frame(do.call(rbind,set.1se))
  colnames(dfr.1se) = c("alpha", "lambda", "df", "cvm")
  min.df.id <- which(dfr.1se$df == min(dfr.1se$df, na.rm = TRUE))
  sub.dfr.1se <- dfr.1se[min.df.id ,]
  if(!is.null(nrow(min.df.id))) {
    
    best.combi.1se <- sub.dfr.1se[which.max(sub.dfr.1se$lam),]
  } else {
    best.combi.1se <- sub.dfr.1se
  }
  
  #step 5: find the best (lowest) of these lambdas
  c(alpha.min = best.combi.min$alpha, lambda.min = best.combi.min$lambda, alpha.1se = best.combi.1se$alpha, lambda.1se= best.combi.1se$lambda)

}