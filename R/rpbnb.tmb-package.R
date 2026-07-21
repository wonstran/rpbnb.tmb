#' rpbnb.tmb: Random-Parameter Bivariate Negative Binomial with TMB
#'
#' Maximum-simulated-likelihood estimation of bivariate random-parameter
#' negative binomial models using Template Model Builder (TMB). Supports
#' Famoye/Sarmanov and discrete-copula (Frank, Gaussian, Clayton) dependence.
#'
#' @docType package
#' @name rpbnb.tmb-package
NULL

#' @importFrom stats nlminb coef vcov logLik AIC BIC predict
#' @importFrom TMB MakeADFun sdreport
NULL

# On-load hook: load the TMB DLL (package name has a dot, so useDynLib
# would look for the wrong R_init_ symbol — load via library.dynam instead).
.onLoad <- function(libname, pkgname) {
  library.dynam("rpbnb.tmb", pkgname, libname)
}
