\name{bootstrapTimeToEvent}
\alias{bootstrapTimeToEvent}
\title{Bootstrap the predicted time when a given number of events is reached, for hybrid Exponential model}
\description{Bootstrap the predicted time when a given number of events is reached, based on the hybrid Exponential model.}
\usage{bootstrapTimeToEvent(time, event, interim.gates, future.units, n, M0 = 1000, 
                            K0 = 5, alpha = 0.05, accrual, seed = 2014)}
\arguments{
  \item{time}{Event times.}
  \item{event}{Censoring indicator, 0 = censored, 1 = event.}
  \item{interim.gates}{Number of events for which confidence intervals should be computed. May be a vector.}
  \item{future.units}{Number of months for which predictions are to be made.}
  \item{n}{Size of the bootstrap samples to be drawn from \code{surv.obj}.}
  \item{M0}{Number of bootstrap samples to be drawn.}
  \item{K0}{Number of changepoints for the piecewise constant hazard.}
  \item{alpha}{Familywise error rate for the sequential test.}
  \item{accrual}{Vector of dates when future patients enter the study. As for the specification, assume 50 pts are to be recruited in Dec 2013 and 20pts in Jan 2014, then use \code{accrual <- c(rep(as.Date("2013-12-01"), 50), rep(as.Date("2014-01-01"), 20))}. Leave as \code{NULL} if accrual for the study is completed.}
  \item{seed}{Seed for generation of bootstrap samples.}
}
\value{A list containing the following objects:
  \item{pe.MLEs}{Piecewise Exponential MLE object for each bootstrap sample, output of function \code{\link{piecewiseExp_MLE}}.}
  \item{pe.tabs}{Result of sequential test for each bootstrap sample, output of function \code{\link{piecewiseExp_test_changepoint}}.}
  \item{changepoints}{For each bootstrap sample, changepoint as resulting from sequential test.}
  \item{estS}{Estimated survival function for each bootstrap sample.}
  \item{event.dates}{Event dates for each bootstrap sample, where columns relate to the number of events in \code{interim.gates}.}
}
\author{
Kaspar Rufibach (maintainer) \cr \email{kaspar.rufibach@roche.com}
}
\references{
Rufibach, K. (2016).
\emph{Event projection: quantify uncertainty and manage expectations of broader teams.}
Slides for talk given in Basel Biometric Section Seminar on 28th April 2016.
\url{http://bbs.ceb-institute.org/wp-content/uploads/2016/06/Kaspar-event_tracking.pdf}.
}

\examples{

\dontrun{
# --------------------------------------------------
# simulate data for illustration
# --------------------------------------------------
set.seed(2021)
n <- 600
time <- rexp(n, rate = log(2) / 20)
event <- sample(round(runif(n, 0, 1)))
accrual_after_ccod <- 1:(n - length(time)) / 30

# --------------------------------------------------
# run bootstrap, for M0 = 3 only, for illustration
# tune parameters for your own example
# --------------------------------------------------
boot1 <- bootstrapTimeToEvent(time, event, 
          interim.gates = c(330, 350), future.units = 50, n = length(time), 
          M0 = 3, K0 = 5, alpha = 0.05, accrual = accrual_after_ccod, 
          seed = 2014)

# median of bootstrap samples:
apply(boot1$event.dates, 2, median)
}
}

\keyword{htest}
\keyword{survival}
