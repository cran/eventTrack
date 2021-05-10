\name{bootSurvivalSample}
\alias{bootSurvivalSample}
\title{Bootstrap survival data}
\description{Generate bootstrap samples from a survival object, for right-censored data. See Efron (1981) and Akritas (1986) for details.}
\usage{bootSurvivalSample(surv.obj, n, M = 1000)}
\arguments{
  \item{surv.obj}{\code{Surv} object containing survival times and censoring indicator.}
  \item{n}{Size of the bootstrap samples to be drawn from \code{surv.obj}.}
  \item{M}{Number of bootstrap samples to be drawn.}
}
\value{Matrix with bootstrap samples as columns.}
\author{
Kaspar Rufibach (maintainer) \cr \email{kaspar.rufibach@roche.com}
}
\references{
Akritas, M.G. (1986).
Bootstrapping the Kaplan-Meier Estimator.
\emph{JASA}, \bold{81}, 1032-1038

Efron, B. (1981).
Censored Data and the Bootstrap.
\emph{JASA}, \bold{76}, 312-319.
}

\examples{
# --------------------------------------------------
# simulate data for illustration
# --------------------------------------------------
set.seed(2021)
n <- 600
time <- rexp(n, rate = log(2) / 20)
event <- sample(round(runif(n, 0, 1)))

# --------------------------------------------------
# draw 20 bootstrap samples of size 10
# --------------------------------------------------
bootSurvivalSample(surv.obj = Surv(time, event), n = 10, M = 20)
}

\keyword{htest}
\keyword{survival}
