\name{predictEventsUncond}
\alias{predictEventsUncond}
\title{compute expected number of events based on a fixed survival function and with no recruited patients yet}
\description{Based on a specified survival function and future accrual, compute for each month in the future the number of expected events reached by then, based on a fixed pre-specified survival function.}
\usage{predictEventsUncond(St, accrual, future.units = 50)}
\arguments{
  \item{St}{Function that specifies the survival function to be used.}
  \item{accrual}{\code{NULL} if trial is fully accrued. For potential future accural, see the example below how to specify it.}
  \item{future.units}{Number of future months to compute prediction for.}
}
\value{A \code{data.frame} with the months and corresponding expected events.}

\author{
Kaspar Rufibach (maintainer) \cr \email{kaspar.rufibach@roche.com}
}
\examples{

# compute date when 380 events are reached:
nevent <- 380
n <- 800
accrual.month <- 50
accrual <- rep(seq(1, 16), each = accrual.month)
rate <- log(2) / 12

St1 <- function(t0){
      res <- 1 - pexp(t0, rate = rate)
      return(res)
}
pred1 <- predictEventsUncond(St = St1, accrual, future.units = 25)
pred1
exactDatesFromMonths(predicted = pred1, nevent)

}

\keyword{htest}
