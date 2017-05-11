var_count <- function(var_sel, p) {
  if(length(var_sel) == 0) {
    varc <- rep(0, p)
  } else {
    varc <- as.numeric(1:p %in% var_sel)
  }
  varc
}