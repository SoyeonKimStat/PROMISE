# 02/24/16
# Stabiilty selection with cross validation
# Choose alpha, lambda, cutoff using cross validation using 1se rule and minimum rule
# using stratified cross-validation
# return alpha, lambda, cutoff after CV
# change 1se rule: when df is the same, select the one with df.se (standard error of nnonzero)
# change 1se rule: 1)smallest df 2)largest cutoffs and lambda 3)smallest df.se

# glmnet fit after variable selection
# lambda = 100
# faster one

library(glmnet)
library(ROCR)
library(plyr)

source("stability.cv.path_any.response.r")
source("glmnet.predicts.r")


# source("stratified.divide.data2.r")
# 
# repi =1;p.train = c(3/4);family = "binomial"; lambda.rule = "lambda.1se"
# n=200
# nm=250;nnonzero=12;between.cor <- 0.2;beta.coef <- 0.5; within.cor <- 0.8;  option = "elastic"
# 
# 
# nsim =".four.trt.normal.mycoding.random.beta"
# 
# data <- sim.data(n) 
# div.data <- s.divide.data(data, p.train)
# x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test;
# nfolds=5;option="elastic"
# 
# test <- stab.cv(x.train,y.train)


# source("stratified.divide.data2.r")
# data <-  read.csv("sub100_array_interaction.csv", header=TRUE)
# repi =1000;p.train = c(3/4);family = "binomial"; lambda.rule = "lambda.1se"
# div.data <- s.divide.data(data, p.train)
# x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test;
# nfolds =5; option ="elastic"



stab.cv <- function(x.train, y.train, family="binomial", nfolds = 5, option = "elastic", nlambda=20, standardize=FALSE) {
  #ptm <- proc.time()
  if(option == "elastic") {
    alphas <- seq(0.1, 1, by=0.1)
  } else if(option == "lasso") {
    alphas <- 1
  } else {
    stop(paste("option=", option, " is not available option", sep="")) 
  }
  
  cutoffs <- seq(0.8, 0.3, by=-0.1); ncutoff <- length(cutoffs)
  nalpha = length(alphas); nlambda <- nlambda; ncanlambda <- nlambda-1
  p <- ncol(x.train);n.train<-length(y.train)
  
  lambda.names <- paste("s",1:ncanlambda,sep="")
  
  seq.lambda <- matrix(nrow = ncanlambda, ncol= nalpha)
  colnames(seq.lambda) <- alphas
  rownames(seq.lambda) <- lambda.names
  
  nonzero.mat <- matrix(NA, nrow=nfolds, ncol= ncutoff)
  rownames(nonzero.mat) <- 1:nfolds
  colnames(nonzero.mat) <- cutoffs
  
  kfold.df <- rep(list(nonzero.mat), ncanlambda)
  names(kfold.df) <-  lambda.names
  
  kfold.loss <- kfold.df
  
  lmean.loss <- vector("list", nalpha)
  lse.loss <- vector("list", nalpha)
  lmedian.df <- vector("list", nalpha)
  lse.df <- vector("list", nalpha)
  
  
  
  se <- function(x) sqrt(var(x)/length(x))
  
  # stratified CV (random sampling was already done when sepeperating train vs test set)
  if(family == "binomial") {
    ny1 <- sum(y.train)
    foldid <- c(rep(1:nfolds,length.out=ny1),rep(1:nfolds,length.out=c(n.train-ny1))) 
  } else {
    foldid = sample(rep(seq(nfolds), length = n.train))
  }
  
  
  
  one.cv.result <- function(alpha) {
    pre.fit <- glmnet(x.train, y.train, alpha = alpha, family = family, nlambda=nlambda, standardize = standardize) #nlambda
    lambda.range <- pre.fit$lambda
    n.real.lambda = length(lambda.range)
    if(n.real.lambda !=ncanlambda) {
      lambda.range = c(rep(100,c(nlambda -n.real.lambda)),lambda.range)
    }
    for(k in 1:nfolds) {
      
      which = foldid == k
      x.sub = x.train[!which,]
      y.sub = y.train[!which]
      x.val = x.train[which,]
      y.val = y.train[which]
      
      phi <- stability.cv.path(x.sub, y.sub, family=family,alpha = alpha, lambda.range) 
      #using grid of cutoffs and lambda, save names of selected variables in a list
      
      for(i in 1:ncanlambda) {
        limit.phi <- phi[1:c(i+1),]
        result <- apply(limit.phi,2,max)    
        # what percentage of lambda has been used 
        cut.values <-sapply(seq_along(cutoffs), function(curi){
          cf<-cutoffs[curi]
          nonzero.ind <- result > cf
          loss <- glmnet.predicts(nonzero.ind, x.sub, y.sub, x.val, y.val, family)
          df<- sum(nonzero.ind)
          c(loss=loss, df=df) # upper curve
        })
        kfold.loss[[i]][k,] <- cut.values["loss",]
        kfold.df[[i]][k,] <- cut.values["df",]
      }
    }
    
    mean.loss <- sapply(1:ncanlambda, function(curi) {
      apply(kfold.loss[[curi]], 2,mean) 
    })
    colnames(mean.loss) <- lambda.names
    
    se.loss <- sapply(1:ncanlambda, function(curi) {
      apply(kfold.loss[[curi]], 2,se)
    })
    colnames(se.loss) <- lambda.names
    
    median.df <- sapply(1:ncanlambda, function(curi) {
      apply(kfold.df[[curi]], 2,median)
    })
    colnames(median.df) <- lambda.names
    
    se.df <- sapply(1:ncanlambda, function(curi) {
      apply(kfold.df[[curi]], 2,se)
    })
    colnames(se.df) <- lambda.names
    
    list(mean.loss=mean.loss, se.loss=se.loss, median.df=median.df, se.df=se.df, lambdas = lambda.range[2:nlambda])
  }
  
  for(a in seq_along(alphas)) {
    print(paste("Running alpha=", alphas[a], sep=""))
    pre.result <- one.cv.result(alphas[a])
    lmean.loss[[a]] <- pre.result$mean.loss
    lse.loss[[a]] <- pre.result$se.loss
    lmedian.df[[a]] <- pre.result$median.df
    lse.df[[a]] <- pre.result$se.df
    seq.lambda[,a] <- pre.result$lambdas
  }
  
  
  
  
  # get maximum losss in each alpha
  seq.best <- lapply(seq_along(alphas), function(curi) {
    mat.mean.loss<-lmean.loss[[curi]]
    mat.se.loss <- lse.loss[[curi]]
    mat.md.nnonzero <- lmedian.df[[curi]]
    mat.se.nnonzero <- lse.df[[curi]]
    one.seq.lambda <- seq.lambda[,curi]
    names(one.seq.lambda) <- NULL
    if(family == "binomial") {
      best.loss <- max(mat.mean.loss)
    } else {
      best.loss <- min(mat.mean.loss)
    } 
    
    best.index <- which(mat.mean.loss == best.loss, arr.ind=T)
    if(nrow(best.index) !=1) { # if the combination is not one take only one
      # Just take the first one
      best.index <- best.index[1,]
    }
    
    # # variables that maximize loss under this alpha
    best.lambdas = one.seq.lambda[best.index[2]] 
    best.cutoffs = cutoffs[best.index[1]]
    best.mean.losss=best.loss
    best.se.losss <- mat.se.loss[best.index[1],best.index[2]]
    best.md.nnonzeros <- mat.md.nnonzero[best.index[1],best.index[2]]
    best.se.nnonzeros <- mat.se.nnonzero[best.index[1],best.index[2]]
    c(alpha =alphas[curi], lambda=best.lambdas, cutoff = best.cutoffs, mean.loss = best.mean.losss, se.loss = best.se.losss,
      md.df = best.md.nnonzeros, se.df=best.se.nnonzeros)
  })
  best <- as.data.frame(do.call(rbind,seq.best))
  if(family == "binomial") {
    best.combi <- best[which.max(best$mean.loss),]
    se.min <- best.combi$mean.loss - best.combi$se.loss
  } else {
    best.combi <- best[which.min(best$mean.loss),]
    se.min <- best.combi$mean.loss + best.combi$se.loss
  }
  
  
  
  
  
  #step 4: For each alpha, which lambda is the best within 1SE difference from the best model
  set.1se<-lapply(seq_along(alphas), function(curi){
    mat.mean.loss<-lmean.loss[[curi]]
    mat.md.nnonzero <- lmedian.df[[curi]]
    mat.se.nnonzero <- lse.df[[curi]]
    
    if(family == "binomial") {
      idmin = which(mat.mean.loss >= se.min & mat.md.nnonzero <= best.combi$md.df, arr.ind=T)
    } else {
      idmin = which(mat.mean.loss <= se.min & mat.md.nnonzero <= best.combi$md.df, arr.ind=T)
    }
    if(any(idmin)) {       # if there are multiple cases whose loss is larger than max.loss-SE 
      #and more parsimoninous than maximum model
      idmin<-data.frame(idmin)
      # find the largest lambda for each cutoff
      smaller.idmin<- ddply(idmin, .(row), summarise,
            min_lambda = min(col, na.rm = TRUE))
     
#       par.row <- which(idmin[,2] == min(idmin[,2])) # 1.For a largest lambda 
#       par.row2 <- which.min(idmin[par.row,1]) # Among them, largest cutoff
#       par.id <- idmin[par.row[par.row2],] # par.id: most parcimonioius id given the condition
#       if(any(idmin[,1] < par.id[1])) { # if there is  model with lower lambda but larger cutoff, this chooses model with the smallest median nonzero
#       
      s.idmin<-as.matrix(smaller.idmin)
      if(nrow(s.idmin) !=1) {
        id.df <-mat.md.nnonzero[s.idmin]
    
        #par.id <- idmin[which.min(id.df),]
        temp.par.id <- s.idmin[which(id.df == min(id.df)),]
        if(!is.null(nrow(temp.par.id))) { # if there are multiple rows whose med.nnonzero is minimum number
          # select with smallest se for nnonzero
          if(length(unique(temp.par.id[,2])) == 1 | length(unique(temp.par.id[,1])) == 1) {
            par.id <- temp.par.id[1,]
          } else {
            id.se.nnonzero <- mat.se.nnonzero[temp.par.id]
            par.id <- temp.par.id[which.min(id.se.nnonzero),]
          }
        } else {
          par.id <- temp.par.id
        }
      } else {
        par.id <- s.idmin
      }
      
      lamda.1se <- seq.lambda[par.id[2],curi]
      cutoff.1se <- cutoffs[par.id[1]]
      loss.1se <-  mat.mean.loss[par.id[1], par.id[2]]
      df.md.1se <- mat.md.nnonzero[par.id[1], par.id[2]]
      df.se.1se <- mat.se.nnonzero[par.id[1], par.id[2]]
    }  else {
      lamda.1se <- NA #fixed part
      cutoff.1se <- NA
      loss.1se <-  NA
      df.md.1se <- NA
      df.se.1se <- NA
    }  
    c(alpha=alphas[curi],lambda=lamda.1se, cutoff= cutoff.1se, mean.loss =  loss.1se, md.df=df.md.1se, se.df =df.se.1se)
  })
  dfr.1se <- as.data.frame(do.call(rbind,set.1se))
  min.df.id <- which(dfr.1se$md.df == min(dfr.1se$md.df, na.rm = TRUE))
  sub.dfr.1se <- dfr.1se[min.df.id ,]
  if(nrow(sub.dfr.1se) != 1) {
    
    best.combi.1se <- sub.dfr.1se[which.min(sub.dfr.1se$se.df),]
  } else {
    best.combi.1se <- sub.dfr.1se
  }
  list(combi.min = best.combi[c("alpha","lambda","cutoff", "mean.loss", "md.df", "se.df")], combi.1se = best.combi.1se)
}