#' rpbnb.tmb: Random-Parameter Bivariate Negative Binomial with TMB
#'
#' Maximum-simulated-likelihood estimation of bivariate random-parameter
#' negative binomial models using Template Model Builder (TMB). Supports
#' Famoye/Sarmanov and discrete-copula (Frank, Gaussian, Clayton) dependence.
#'
#' @keywords internal
#' @importFrom stats AIC BIC coef dnbinom logLik nlminb pnbinom pnorm predict
#' @importFrom stats qnbinom qnorm rnbinom rnorm runif symnum vcov
#' @importFrom numDeriv jacobian
#' @importFrom TMB MakeADFun
#' @useDynLib rpbnb.tmb, .registration = TRUE
"_PACKAGE"
