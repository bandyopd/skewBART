MakeForest <- function(hypers, opts) {
  mf <- Module(module = "mod_forest", PACKAGE = "MVBart")
  return(new(mf$Forest, hypers, opts))
}
