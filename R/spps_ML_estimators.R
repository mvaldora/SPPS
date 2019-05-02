#' Estimates epsilon, knowing delta and  eta=X%*%beta, by maximum likelihood.
#' @param a A vector a of 0s and 1s of dimension n
#' @param eta A vector of linear predictors of dimension n
#' @param delta A real number between 0 and 1
#' @param phi A function that maps the real line onto the interval (0,1)
#' @return a real number between 0 and 1
#' @export
estiepsML <- function(a, eta, delta, phi = expit) {
  epsmin <- 0
  epsmax <- 1
  phispb <-
    function(t, epsilon, delta)
      (1 - delta - epsilon) * phi(t) + epsilon
  Leps <- function(epsilon) {
    grandes <- which(phispb(eta, epsilon, delta) > 0.9999)
    chicos  <- which(phispb(eta, epsilon, delta) < 0.0001)
    phiesspb <- phispb(eta, epsilon, delta)
    phiesspb[grandes] <- 0.9999
    phiesspb[chicos] <- 0.0001
    phispbev <- phispb(eta, epsilon, delta)
    ans <-
      sum((1 - phi(eta)) * (a / phispbev - (1 - a) / (1 - phiesspb)))
    ans
  }
  Leps <- Vectorize(Leps)
  if (Leps(epsmin) * Leps(epsmax) < 0) {
    eps1 <- uniroot(Leps, c(epsmin, epsmax))$root
  } else{
    if (abs(Leps(epsmin)) < abs(Leps(epsmax))) {
      eps1 <- epsmin
    } else{
      eps1 <- epsmax
    }
  }
  eps1
}
#' Estimates delta, knowing epsilon and  eta=X%*%beta, by maximum likelihood.
#' @param a A Vector a of 0s and 1s of dimension n
#' @param eta A vector of linear predictors of dimension n
#' @param epsilon A real number between 0 and 1
#' @param phi A function that maps the real line onto the interval (0,1)
#' @return a real number between 0 and 1
#' @export
estidelML <- function(a, eta, epsilon, phi = expit) {
  delmin <- 0
  delmax <- 0.99
  phispb <-
    function(t, epsilon, delta)
      (1 - delta - epsilon) * phi(t) + epsilon
  Ldel <- function(delta) {
    phispbev <- phispb(eta, epsilon, delta)
    sum((-phi(eta)) * (a / phispbev - (1 - a) / (1 - phispbev)))
  }
  Ldel <- Vectorize(Ldel)
  if (Ldel(delmin) * Ldel(delmax) < 0) {
    del1 <- uniroot(Ldel, c(delmin, delmax))$root
  } else{
    if (abs(Ldel(delmin)) < abs(Ldel(delmax))) {
      del1 <- delmin
    } else{
      del1 <- delmax
    }
  }
  del1
}
#' Estimates beta, knowing epsilon and delta, by maximum likelihood.
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
estibetaML <- function(x,
           a,
           epsilon,
           delta,
           phi = expit,
           dphi = dexpit,
           phiinv = expitinv,
           maxit = 3) {
    HH <- function(t)
      t + (1 - t) * (epsilon + delta) - delta
    HHinv <- function(s)
      (s - epsilon) / (1 - epsilon - delta)
    expitacot <- function(x)
      HH(expit(x))
    expitacotinv <- function(y)
      expitinv(HHinv(y))
    dexpitacot <- function(x)
      (1 - epsilon - delta) * dexpit(x)
    logitacot <- function() {
      linkfun <- function(mu)
        expitacotinv(mu)
      linkinv <- function(eta)
        expitacot(eta)
      mu.eta <- function(eta)
        dexpitacot(eta)
      valideta <- function(eta)
        TRUE
      link <- "logitacot"
      structure(
        list(
          linkfun = linkfun,
          linkinv = linkinv,
          mu.eta = mu.eta,
          valideta = valideta,
          name = link
        ),
        class = "link-glm"
      )
    }
    ans <- list()
    ran <- 1 - delta - epsilon
    pi0 <-
      (epsilon + ran / 4) * (a == 0) + (1 - delta - ran / 4) * (a == 1)
    control <-
      glm.control(epsilon = 0.0001,
                  maxit = maxit,
                  trace = FALSE)
    ajuste <-
      glm.fit(
        cbind(1, x),
        a,
        family = binomial(link = logitacot()),
        mustart = pi0,
        control = control
      )
    ans <- ajuste
    ans
  }
#' Estimates all the parameters in the SPPS model by maximum likelihood
#' @param x A matrix of covariates of dimension n x p
#' @param a A Vector a of 0s and 1s of dimension n
#' @param method The method to be used "ATE" (default) or "IPW"
#' @param phi A function that maps the real line onto the interval (0,1)
#' @param dphi The derivative of \code{phi}
#' @param phiinv The inverse of \code{phi}
#' @param maxit The maximum number of iterations to be performed. 
#' @return a list of model output values as in glm.fit
#' @export
estiML <- function(x,
                   a,
                   method = "ATE",
                   phi = expit,
                   dphi = dexpit,
                   phiinv = expitinv,
                   maxit = 20) {
  phispb <- function(t, epsilon, delta) {
    (1 - delta - epsilon) * phi(t) + epsilon
  }
  phispbinv <- function(y, epsilon, delta) {
    phiinv((y - epsilon) / (1 - epsilon - delta))
  }
  p <- ncol(x)
  ajuste0 <- try(estibetaML(
    x = x,
    a = a,
    epsilon = 0,
    delta = 0,
    phi = expit,
    dphi = dexpit,
    phiinv = expitinv
  ),
  TRUE)
  beta0 <- ajuste0$coefficients
  pig0 <- ajuste0$fitted.values
  eps0 <- min(pig0)
  if (method == "ATE") {
    del0 <- 1 - max(pig0)
  } else{
    del0 <- 0
  }
  eps1 <- 0
  del1 <- 0
  pig1 <- rep(0, nrow(x))
  beta1 <- rep(0, p + 1)
  dispig <- sum((pig0 - pig1) ^ 2)
  ran <- 1 - eps0 - del0
  contador <- 0
  while ((dispig > 0.001) & (ran > 0.5) & contador < maxit) {
    contador <- contador + 1
    if (method == "ATE") {
      ajuste1 <-
        estibetaML(
          x = x,
          a = a,
          epsilon = eps0,
          delta = del0,
          phi = phi,
          dphi = dphi,
          phiinv = phiinv,
          maxit = 3
        )
    }
    if (method == "IPW") {
      ajuste1 <-
        estibetaML(
          x = x,
          a = a,
          epsilon = eps0,
          delta = 0,
          phi = phi,
          dphi = dphi,
          phiinv = phiinv,
          maxit = 3
        )
    }
    beta1 <- ajuste1$coefficients
    pig1 <- ajuste1$fitted
    eta1 <- phispbinv(pig1, eps0, del0)
    eps1 <- estiepsML(a, eta = eta1, delta = del0, phi = phi)
    if (method == "ATE")
      del1 <- estidelML( a, eta = eta1, epsilon = eps1, phi = phi)
    if (method == "IPW")
      del1 <- 0
    dispig <- sum((pig0 - pig1) ^ 2)
    eps0 <- eps1
    del0 <- del1
    beta0 <- beta1
    pig0 <- pig1
    ajuste0 <- ajuste1
  }
  if (contador == maxit)
    print("no convergence")
  ans <- list()
  ans$lm <- ajuste0
  ans$beta <- beta0
  ans$epsilon <- eps0
  if (method == "ATE")
    ans$delta <- del0
  ans$fitted <- ajuste0$fitted
  ans
}



