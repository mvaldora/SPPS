#' This function computes exp(x) / (1 + exp(x)), the inverse of the logistic function log(p / (1 - p))
#'
#' @param x A real number
#' @return A real number between 0 and 1
#' @export
expit    <- function(x)
  exp(x) / (1 + exp(x))
#' This function computes the derivatibe of the expit function
#'
#' @param x A real number
#' @return A real number 
#' @export
dexpit<-function(x){exp(x)/(1+exp(x))^2}
#' This function computes  the logistic function log(p / (1 - p))
#'
#' @param p A real number between 0 and 1 
#' @return A real number
#' @export
expitinv <- function(p)
  log(p / (1 - p))
#' This function computes \pi(x) for the SPPS model
#' @param x A real number x
#' @param epsilon A real number between 0 and 1
#' @param delta A real number between 0 and 1
#' @return A real number
#' @export
expacotm <- function(x, epsilon, delta) {
  ans <- (1 - epsilon - delta) * expit(x) + epsilon
  grandes <- which(x > 100)
  ans[grandes] <- rep(1 - delta, length(grandes))
  ans
}
#' This function computes the derivative of pi(x) for the SPPS model
#' @param x A real number 
#' @param epsilon A real number between 0 and 1
#' @param delta A real number between 0 and 1
#' @return A real number
#' @export
dexpacotm <- function(x, epsilon, delta)
  (1 - epsilon - delta) * dexpit(x)
#' Computes the inverse of the function pi(x) for  the spps  model.
#' @param y A real number between 0 and 1
#' @param espilon A real number between 0 and 1 such that \code{epsilon} + \code{delta} < 1
#' @param delta A real number between 0 and 1 
#' @return A real number
#' @export
expacotminv <-
  function(y, epsilon, delta)
    expitinv((y - epsilon) / (1 - epsilon - delta))


