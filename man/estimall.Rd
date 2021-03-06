% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spps_final_estimators.R
\name{estimall}
\alias{estimall}
\title{Estimates all the parameters in the SPPS model}
\usage{
estimall(x, a, method = "ATE", estim.method = "ML", phi = expit,
  dphi = dexpit, phiinv = expitinv, maxit = 20)
}
\arguments{
\item{x}{A matrix of covariates of dimension n x p}

\item{a}{A Vector a of 0s and 1s of dimension n}

\item{method}{The method to be used "ATE" (default) or "IPW"}

\item{estim.method}{The estimation method to be used for estimating the parameters "ML" (default) or "LS"}

\item{phi}{A function that maps the real line onto the interval (0,1)}

\item{dphi}{The derivative of \code{phi}}

\item{phiinv}{The inverse of \code{phi}}

\item{maxit}{The maximum number of iterations to be performed.}
}
\value{
a list of model output values as in glm.fit
}
\description{
Estimates all the parameters in the SPPS model
}
