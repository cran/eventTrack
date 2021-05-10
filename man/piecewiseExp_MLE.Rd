\name{piecewiseExp_MLE}
\alias{piecewiseExp_MLE}
\title{Estimate hazard function in piecewise Exponential survival model}
\description{This function estimates the values of the hazard function and the change points in a piecewise Exponential survival model. The number of change points \eqn{K} needs to be pre-specified.}
\usage{piecewiseExp_MLE(time, event, K)}
\arguments{
  \item{time}{Event times, censored or observed, in months.}
  \item{event}{Censoring indicator, 1 for event, 0 for censored.}
  \item{K}{Number of change points to be used in the model.}
}
\value{A list containing the following objects:
  \item{tau}{The maximum likelihood estimates of the change points.}
  \item{lambda}{The maximum likelihood estimates of the value of the hazard function.}
  \item{lambda.SE}{The standard errors of \code{lambda}.}
  \item{K}{Number of change points to be used in the model.}
}

\author{
Kaspar Rufibach (maintainer) \cr \email{kaspar.rufibach@roche.com}
}
\references{
Fang, L., Zheng, S. (2011).
A hybrid approach to predicting events in clinical trials with time-to-event outcomes.
\emph{Contemp. Clin. Trials}, \bold{32}, 755--759. 

Goodman, M.S., Li, Y., Tiwari, R.C. (2011).
Detecting multiple change points in piecewise constant hazard functions.
\emph{J. Appl. Stat}, \bold{38(11)}, 2523--2532.
}

\examples{# see vignette}

\keyword{htest}
