#' Estimates epsilon, knowing delta and  eta=X%*%beta, by least squares
#' @param a A Vector a of 0s and 1s of dimension n
#' @param eta A vector of linear predictors of dimension n
#' @param delta A real number between 0 and 1
#' @param phi A function that maps the real line onto the interval (0,1)
#' @return a real number between 0 and 1
#' @export
estiepsLS <- function(a, eta, delta, phi = expit) {
  phiev <- phi(eta)
  max(sum(a * (1 - phiev) - (1 - delta) * phiev * (1 - phiev)) / sum((1 -
                                                                        phiev) ^ 2), 0)
}
#' Estimates delta, knowing epsilon and  eta=X%*%beta, by least squares.
#' @param a A Vector a of 0s and 1s of dimension n
#' @param eta A vector of linear predictors of dimension n
#' @param epsilon A real number between 0 and 1
#' @param phi A function that maps the real line onto the interval (0,1)
#' @return a real number between 0 and 1
#' @export
estidelLS <- function(a, eta, epsilon, phi = expit) {
  phiev <- phi(eta)
  max(sum((epsilon) * phiev + (1 - epsilon) * phiev ^ 2 - a * phiev) / sum(phiev ^ 2), 0)
}
#' Estimates beta, knowing epsilon and delta, by least squares.
#' @param x A matrix of covariates of dimension n x p
#' @param a A vector a of 0s and 1s of dimension n
#' @param epsilon A real number between 0 and 1 such that \code{epsilon}+\code{delta}<1
#' @param delta A real number between 0 and 1
#' @param phi A function that maps the real line onto the interval (0,1)
#' @param dphi The derivative of \code{phi}
#' @param phiinv The inverse of \code{phi}
#' @param maxit The maximum number of iterations to be performed. 
#' @return a real number between 0 and 1
#' @export
estibetaLS <-
  function(x,
           a,
           epsilon,
           delta,
           beta0 = NULL,
           phi = expit,
           dphi = dexpit,
           phiinv = expitinv,
           maxit = 3) {
    control = nls.control(
      maxiter = maxit,
      tol = 1e-05,
      minFactor = 1 / 1024,
      printEval = FALSE,
      warnOnly = TRUE
    )
    p <- ncol(x) + 1
    if (is.null(beta0) == TRUE) {
      beta0 <- rep(0, p)
    }
    ajuste <-
      nls(
        a ~ pifunctionSPPSATEMatrixforbeta(b, epsilon, delta, x),
        start = list(b = beta0),
        control = control
      )
    ans <- list()
    ans$fitted <- fitted(ajuste)
    ans$coefficients <- coef(ajuste)
    ans$residuals <- residuals(ajuste)
    ans
  }
#' Estimates all the parameters in the SPPS model by least squares
#' @param x A matrix of covariates of dimension n x p
#' @param a A Vector a of 0s and 1s of dimension n
#' @param method The method to be used "ATE" (default) or "IPW"
#' @param phi A function that maps the real line onto the interval (0,1)
#' @param dphi The derivative of \code{phi}
#' @param phiinv The inverse of \code{phi}
#' @param maxit The maximum number of iterations to be performed. 
#' @return a list of model output values as in glm.fit
#' @export
estiLS <-
  function(x,
           a,
           method = "ATE",
           phi = expit,
           dphi = dexpit,
           phiinv = expitinv,
           maxit = 20) {
    n <- length(a)
    pi0ini <- (1 / 4) * (a == 0) + (3 / 4) * (a == 1)
    eta0ini <- phiinv(pi0ini)
    eps0ini <- estiepsLS( a, eta0ini, 1 / 4, phi = phi)
    del0ini <- estidelLS( a, eta0ini, 1 / 4, phi = phi)
    phispb <- function(t, epsilon, delta)
      (1 - delta - epsilon) * phi(t) + epsilon
    dphispb <- function(x, epsilon, delta)
      (1 - epsilon - delta) * dphi(x)
    phispbinv <-
      function(y, epsilon, delta)
        phiinv((y - epsilon) / (1 - epsilon - delta))
    p <- ncol(x)
    XI <- cbind(1, x)
    ajuste0 <-
      try(estibetaLS(
        x = x,
        a = a,
        epsilon = eps0ini,
        delta = del0ini,
        phi = expit,
        dphi = dexpit,
        phiinv = expitinv
      ),
      TRUE)
    beta0 <- ajuste0$coefficients
    pig0 <- ajuste0$fitted
    eps00 <- min(pig0)
    eps0 <- eps00
    if (method == "ATE") {
      del00 <- 1 - max(pig0)
    } else{
      del00 <- 0
    }
    del0 <- del00
    res0 <- ajuste0$residuals
    res0 <- residuals(ajuste0, type = "pearson")
    eps1 <- 0
    del1 <- 0
    res1 <- rep(0, nrow(x))
    pig1 <- rep(0, nrow(x))
    beta1 <- rep(0, p + 1)
    dise <- abs(eps0 - eps1)
    disd <- abs(del0 - del1)
    disres <- sum(abs(res0 - res1))
    dispig <- sum((pig0 - pig1) ^ 2)
    ran <- 1 - eps0 - del0
    contador = 0
    while ((dispig > 0.001) & (ran > 0.5) & contador < maxit) {
      contador <- contador + 1
      if (method == "ATE") {
        ajuste1 <-
          estibetaLS(
            x = x,
            a = a,
            epsilon = eps0,
            delta = del0,
            phi = expit,
            dphi = dexpit,
            phiinv = expitinv
          )
      }
      if (method == "IPW") {
        ajuste1 <-
          estibetaLS(
            x = x,
            a = a,
            epsilon = eps0,
            delta = 0,
            phi = expit,
            dphi = dexpit,
            phiinv = expitinv
          )
      }
      res1 <- residuals(ajuste1, type = "pearson")
      beta1 <- coef(ajuste1)
      pig1 <- fitted(ajuste1)
      #eps10<-min(eps0,min(pig1));del10<-min(min(1-pig1),del0)
      eta1 <- phispbinv(pig1, eps0, del0)
      eps1 <- estiepsLS( a, eta = eta1, delta = del0, phi = phi)
      if (method == "ATE") {
        del1 <- max(estidelLS( a, eta = eta1, epsilon = eps1, phi = phi), del00)
      }
      if (method == "IPW") {
        del1 <- 0
      }
      dise <- abs(eps0 - eps1)
      disd <- abs(del0 - del1)
      disres <- sum(abs(res0 - res1))
      dispig <- sum((pig0 - pig1) ^ 2)
      eps0 <- eps1
      del0 <- del1
      res0 <- res1
      beta0 <- beta1
      pig0 <- pig1
      ajuste0 <- ajuste1
    }
    if (contador == maxit) {
      print("no convergence")
    }
    ans <- list()
    ans$beta <- beta0
    ans$epsilon <- eps0
    
    if (method == "ATE")
      ans$delta <- del0
    ans$fitted <- fitted(ajuste0)
    if (ran <= 0.5)
      ans$beta <- NA
    ans
  }
