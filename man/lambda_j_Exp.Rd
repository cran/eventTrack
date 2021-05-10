\name{lambda_j_Exp}
\alias{lambda_j_Exp}
\title{Compute lambda_j}
\description{For given change points \eqn{\tau_1, \ldots, \tau_k}, compute the profile maximum likelihood estimates \eqn{\lambda_1, \ldots, \lambda_{k + 1}} for the values of the piecewise constant hazard function in a piecewise Exponential survival model. Standard errors for these estimates are provided as well, based on standard maximum (profile) maximum likelihood theory.}
\usage{lambda_j_Exp(tau, time, event)}
\arguments{
  \item{tau}{Single number or vector with values of change points.}
  \item{time}{Event times, censored or observed, in months.}
  \item{event}{Censoring indicator, 1 for event, 0 for censored.}
}
\value{A list containing the following elements:
\item{aj}{The esimtates of the \eqn{\lambda_j}.}
\item{hess.diag}{The diagonal of the corresponding Hessian matrix. Note that the off-diagonal elements of the Hessian are all equal to 0.}
}
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
