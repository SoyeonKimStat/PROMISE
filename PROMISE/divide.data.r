# train.portion : what percent of data is used for train set
# response variable is in the last column
divide.data <- function(data, train.portion) {
  data <- apply(data, 2, as.double)
  n <- nrow(data)
  ntr <- floor(n*train.portion)
  trainid <- sample(seq(n), size = ntr)
  ntrainid <- !seq(n)%in%trainid
  x.train <- data[trainid, -ncol(data)]
  y.train <- data[trainid, ncol(data)]
  x.test <- data[ntrainid, -ncol(data)]
  y.test <- data[ntrainid, ncol(data)]
  p <- ncol(data) - 1
  list(x.train = x.train, y.train = y.train, x.test = x.test, y.test = y.test, p = p)
}
# x.train=div.data$x.train;y.train= div.data$y.train; x.test= div.data$x.test;  y.test= div.data$y.test; p= div.data$p;