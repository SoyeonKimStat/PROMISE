# find AUC scores
library(ROCR)
auc <- function(predic, y) {
  pred <- prediction(predic, y)
  auc <- as.numeric(performance(pred, measure = "auc", x.measure = "cutoff")@y.values)
  auc
}