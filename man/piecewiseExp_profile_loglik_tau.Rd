\name{piecewiseExp_profile_loglik_tau}
\alias{piecewiseExp_profile_loglik_tau}
\title{Profile maximum log-likelihood function for change points in piecewise Exponential survival model}
\description{In a piecewise Exponential survival model, the estimates for the value of the hazard function can be given explicitly, see Fang and Su (2011), and can be computed using \code{\link{lambda_j_Exp}}. The values of the change points can then be computed using the profile log-likelihood function. The function \code{piecewiseExp_profile_loglik_tau} computes the value of this profile log-likelihood function and can be used together with \code{\link{optim}} to compute the change points \eqn{\tau_1, \ldots, \tau_{k}}.}
\usage{piecewiseExp_profile_loglik_tau(tau, time, event)}
\arguments{
  \item{tau}{Single number or vector with values of change points.}
  \item{time}{Event times, censored or observed, in months.}
  \item{event}{Censoring indicator, 1 for event, 0 for censored.}
}
\value{Value of the profile likelihood function.}
\author{
Kaspar Rufibach (maintainer) \cr \email{kaspar.rufibach@roche.com}
}
\note{This function is not intended to be invoked by the end user.}
\references{
Fang, L., Zheng, S. (2011).
A hybrid approach to predicting events in clinical trials with time-to-event outcomes.
\emph{Contemp. Clin. Trials}, \bold{32}, 755--759. 

Goodman, M.S., Li, Y., Tiwari, R.C. (2011).
Detecting multiple change points in piecewise constant hazard functions.
\emph{J. Appl. Stat}, \bold{38(11)}, 2523--2532.
}
\keyword{htest}
