# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rmvnorm <- function(mean, Precision) {
    .Call(`_LVBart_rmvnorm`, mean, Precision)
}

rlgam <- function(shape) {
    .Call(`_LVBart_rlgam`, shape)
}

rcpparma_hello_world <- function() {
    .Call(`_LVBart_rcpparma_hello_world`)
}

rcpparma_outerproduct <- function(x) {
    .Call(`_LVBart_rcpparma_outerproduct`, x)
}

rcpparma_innerproduct <- function(x) {
    .Call(`_LVBart_rcpparma_innerproduct`, x)
}

rcpparma_bothproducts <- function(x) {
    .Call(`_LVBart_rcpparma_bothproducts`, x)
}

update_b <- function(Z, b_new, mu_new, b_old, mu_old, cluster_idx, sigma) {
    .Call(`_LVBart_update_b`, Z, b_new, mu_new, b_old, mu_old, cluster_idx, sigma)
}

