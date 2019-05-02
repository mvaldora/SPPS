
#estimates all the parameters in the SPPS model for missing data
#(method="IPW") and for treatment assignment (method="ATE")
#by least squares (estim.method="LS") or by maximum likelihood (estim.method="ML")
#' Estimates all the parameters in the SPPS model
#' @param x A matrix of covariates of dimension n x p
#' @param a A Vector a of 0s and 1s of dimension n
#' @param method The method to be used "ATE" (default) or "IPW"
#' @param estim.method The estimation method to be used for estimating the parameters "ML" (default) or "LS"
#' @param phi A function that maps the real line onto the interval (0,1)
#' @param dphi The derivative of \code{phi}
#' @param phiinv The inverse of \code{phi}
#' @param maxit The maximum number of iterations to be performed. 
#' @return a list of model output values as in glm.fit
#' @export
estimall <-
  function(x,
           a,
           method = "ATE",
           estim.method = "ML",
           phi = expit,
           dphi = dexpit,
           phiinv = expitinv,
           maxit = 20) {
    if (estim.method == "ML") {
      ans <-
        estiML(
          x,
          a,
          method = method,
          phi = phi,
          dphi = dphi,
          phiinv = phiinv,
          maxit = maxit
        )
    }
    if (estim.method == "LS") {
      ans <-
        estiLS(
          x,
          a,
          method = method,
          phi = phi,
          dphi = dphi,
          phiinv = phiinv,
          maxit = maxit
        )
    }
    
    ans
  }
#' Estimates beta, knowing epsilon and delta.
#' @param x A matrix of covariates of dimension n x p
#' @param a A vector a of 0s and 1s of dimension n
#' @param epsilon A real number between 0 and 1 such that \code{epsilon}+\code{delta}<1
#' @param delta A real number between 0 and 1
#' @param estim.method The estimation method to be used for estimating the parameters "ML" (default) or "LS"
#' @param phi A function that maps the real line onto the interval (0,1)
#' @param dphi The derivative of \code{phi}
#' @param phiinv The inverse of \code{phi}
#' @param maxit The maximum number of iterations to be performed. 
#' @return a real number between 0 and 1
#' @export
estibeta <-
  function(x,
           a,
           epsilon,
           delta,
           estim.method = "ML",
           phi = expit,
           dphi = dexpit,
           phiinv = expitinv,
           maxit = 3) {
    ini <-
      estibetaML(
        x,
        a,
        epsilon = epsilon,
        delta = delta,
        phi = phi,
        dphi = dphi,
        phiinv = phiinv,
        maxit = maxit
      )
    if (estim.method == "ML") {
      ans <- ini
    }
    if (estim.method == "LS") {
      ans <-
        estibetaLS(
          x,
          a,
          epsilon,
          delta,
          beta0 = ini$coefficients,
          phi = phi,
          dphi = dphi,
          phiinv = phiinv,
          maxit = maxit
        )
    }
    ans
  }
#' Estimates the avarage treatment effect (ATE)  using the SPPS model
#' @param x A matrix of covariates of dimension n x p
#' @param a A vector a of 0s and 1s of dimension n
#' @param y A vector of responses or outcomes
#' @param phi A function that maps the real line onto the interval (0,1)
#' @param dphi The derivative of \code{phi}
#' @param estim.method The estimation method to be used for estimating the parameters "ML" (default) or "LS"
#' @param phiinv The inverse of \code{phi}
#' @param maxit The maximum number of iterations to be performed. 
#' @param type The type of IPW estimator to be used: type 1 , type 2 and type 3  combine SPPS with, respectively, IPW1, IPW2 and IPW3 in Lunceford and Davidian(2004)
#' @return a list containing the ATE and the estimators of \code{beta}, \code{epsilon} and \code{delta}
#' @export
atemod <-
  function(x,
           a,
           y,
           phi = expit,
           dphi = dexpit,
           phiinv = expitinv,
           estim.method = "ML",
           maxit = 25,
           type = 1)
  {
    p <- ncol(x)
    ajuste <-
      try(estimall(
        x = x,
        a = a,
        phi = phi,
        dphi = dphi,
        phiinv = phiinv,
        maxit = maxit,
        estim.method = estim.method
      ),
      TRUE)
    pred <- ajuste$fitted
    if (type == 1) {
      dis1 <- mean(y * a / pred)
      dis2 <- mean(y * ((1 - a) / (1 - pred)))
      dis <- dis1 - dis2
    }
    if (type == 2) {
      dis1 <- sum(y * a / pred) / sum(a / pred)
      dis2 <-
        sum(y * ((1 - a) / (1 - pred))) / sum((1 - a) / (1 - pred))
      dis <- dis1 - dis2
    }
    if (type == 3) {
      c1 <- sum((a - pred) / pred) / sum((a - pred) ^ 2 / pred ^ 2)
      cn <- sum(a / pred * (1 - c1 / pred))
      dis1 <- 1 / (cn) * sum(y * (a / pred) * (1 - c1 / pred))
      c0 <-
        -sum((a - pred) / (1 - pred)) / sum((a - pred) ^ 2 / (1 - pred) ^
                                              2)
      cn0 <- sum((1 - a) / (1 - pred) * (1 - c0 / (1 - pred)))
      dis2 <-
        1 / cn0 * sum(y * (1 - a) / (1 - pred) * (1 - c0 / (1 - pred)))
      dis <- dis1 - dis2
    }
    ans <- list()
    ans$ate <- dis
    ans$epsilon <- ajuste$epsilon
    ans$delta <- ajuste$delta
    ans$beta <- ajuste$beta
    ans
  }
#' Estimates the avarage treatment effect (ATE)  using the  classical model (epsilon=delta=0)
#' @param x A matrix of covariates of dimension n x p
#' @param a A vector a of 0s and 1s of dimension n
#' @param y A vector of responses or outcomes
#' @param phi A function that maps the real line onto the interval (0,1)
#' @param dphi The derivative of \code{phi}
#' @param estim.method The estimation method to be used for estimating the parameters "ML" (default) or "LS"
#' @param phiinv The inverse of \code{phi}
#' @param maxit The maximum number of iterations to be performed. 
#' @param type The type of IPW estimator to be used: type 1 , type 2 and type 3  combine SPPS with, respectively, IPW1, IPW2 and IPW3 in Lunceford and Davidian(2004)
#' @return a list containing the ATE and the estimator of \code{beta}
#' @export
ateclas <-
  function(x,
           a,
           y,
           phi = expit,
           dphi = dexpit,
           phiinv = expitinv,
           estim.method = "ML",
           maxit = 25,
           type = 1)
  {
    p <- ncol(x)
    ajuste <-
      estibeta(
        x = x,
        a = a,
        phi = phi,
        epsilon = 0,
        delta = 0,
        estim.method = estim.method,
        dphi = dphi,
        phiinv = phiinv,
        maxit = maxit
      )
    pred <- ajuste$fitted
    if (type == 1) {
      dis1 <- mean(y * a / pred)
      dis2 <- mean(y * ((1 - a) / (1 - pred)))
      dis <- dis1 - dis2
    }
    if (type == 2) {
      dis1 <- sum(y * a / pred) / sum(a / pred)
      dis2 <-
        sum(y * ((1 - a) / (1 - pred))) / sum((1 - a) / (1 - pred))
      dis <- dis1 - dis2
    }
    if (type == 3) {
      c1 <- sum((a - pred) / pred) / sum((a - pred) ^ 2 / pred ^ 2)
      cn <- sum(a / pred * (1 - c1 / pred))
      dis1 <- 1 / (cn) * sum(y * (a / pred) * (1 - c1 / pred))
      c0 <-
        -sum((a - pred) / (1 - pred)) / sum((a - pred) ^ 2 / (1 - pred) ^
                                              2)
      cn0 <- sum((1 - a) / (1 - pred) * (1 - c0 / (1 - pred)))
      dis2 <-
        1 / cn0 * sum(y * (1 - a) / (1 - pred) * (1 - c0 / (1 - pred)))
      dis <- dis1 - dis2
    }
    ans <- list()
    ans$weights <- 1 / pred
    ans$ate <- dis
    ans$beta <- ajuste$coefficients
    ans
  }
#' Estimates the mean of an outcome with missing values using the SPPS model
#' @param x A matrix of covariates of dimension n x p
#' @param a A vector a of 0s and 1s of dimension n
#' @param y A vector of responses or outcomes
#' @param phi A function that maps the real line onto the interval (0,1)
#' @param dphi The derivative of \code{phi}
#' @param phiinv The inverse of \code{phi}
#' @param type The type of IPW estimator to be used: type 1 , type 2 and type 3  combine SPPS with, respectively, IPW1, IPW2 and IPW3 in Lunceford and Davidian(2004)
#' @param estim.method The estimation method to be used for estimating the parameters "ML" (default) or "LS"
#' @return a list containing the estimate of the mean and the estimators of \code{beta} and \code{epsilon}
#' @export
ipwmod <-
  function(x,
           a,
           y,
           phi = expit,
           dphi = dexpit,
           phiinv = expitinv,
           type = 1,
           estim.method = "ML")
  {
    y[which(is.na(y))] <- 0
    p <- ncol(x)
    ajuste <-
      estimall(
        x = x,
        a = a,
        method = "IPW",
        estim.method = estim.method,
        phi = phi,
        dphi = dphi,
        phiinv = phiinv,
        maxit = 25
      )
    pred <- ajuste$fitted
    if (type == 1) {
      ipw <- mean(y * a / pred)
    }
    if (type == 2) {
      ipw <- sum(y * a / pred) / sum(a / pred)
    }
    if (type == 3) {
      c1 <- sum((a - pred) / pred) / sum((a - pred) ^ 2 / pred ^ 2)
      cn <- sum(a / pred * (1 - c1 / pred))
      ipw <- 1 / (cn) * sum(y * (a / pred) * (1 - c1 / pred))
    }
    ans <- list()
    ans$ipw <- ipw
    ans$epsilon <- ajuste$epsilon
    ans$beta <- ajuste$beta
    ans
  }
#' Estimates the avarage treatment effect (ATE)  using the  classical model (epsilon=delta=0)
#' @param x A matrix of covariates of dimension n x p
#' @param a A vector a of 0s and 1s of dimension n
#' @param y A vector of responses or outcomes
#' @param phi A function that maps the real line onto the interval (0,1)
#' @param dphi The derivative of \code{phi}
#' @param estim.method The estimation method to be used for estimating the parameters "ML" (default) or "LS"
#' @param phiinv The inverse of \code{phi}
#' @param maxit The maximum number of iterations to be performed. 
#' @param type The type of IPW estimator to be used: type 1 , type 2 and type 3  combine SPPS with, respectively, IPW1, IPW2 and IPW3 in Lunceford and Davidian(2004)
#' @return a list containing the estimate of the mean and the estimator of \code{beta}
#' @export
ipwclas <-
  function(x,
           a,
           y,
           phi = expit,
           dphi = dexpit,
           phiinv = expitinv,
           estim.method = "ML",
           type = 1)
  {
    y[which(is.na(y))] <- 0
    p <- ncol(x)
    ajuste <-
      estibeta(
        x = x,
        a = a,
        phi = phi,
        epsilon = 0,
        delta = 0,
        estim.method = estim.method,
        dphi = dphi,
        phiinv = phiinv,
        maxit = 3
      )
    pred <- ajuste$fitted
    if (type == 1) {
      ipw <- mean(y * a / pred)
    }
    if (type == 2) {
      ipw <- sum(y * a / pred) / sum(a / pred)
    }
    if (type == 3) {
      c1 <- sum((a - pred) / pred) / sum((a - pred) ^ 2 / pred ^ 2)
      cn <- sum(a / pred * (1 - c1 / pred))
      ipw <- 1 / (cn) * sum(y * (a / pred) * (1 - c1 / pred))
    }
    ans <- list()
    ans$ipw <- ipw
    ans$beta <- ajuste$coefficients
    ans
  }
