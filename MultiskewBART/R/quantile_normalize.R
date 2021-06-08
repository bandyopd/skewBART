quantile_normalize <- function(X, test_X) {
  X_trans <- skewBART::quantile_normalize_bart(rbind(X, test_X))
  idx_train <- 1:nrow(X)
  X <- X_trans[idx_train,,drop=FALSE]
  test_X <- X_trans[-idx_train,,drop=FALSE]
  return(list(X = X, test_X = test_X))
}