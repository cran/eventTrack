\name{kaplanMeier_at_t0}
\alias{kaplanMeier_at_t0}
\title{Compute value of Kaplan-Meier estimate at a given time}
\description{Compute value of Kaplan-Meier estimate at a given time \eqn{t_0}.}
\usage{kaplanMeier_at_t0(time, event, t0)}
\arguments{
  \item{time}{Event times, censored or observed.}
  \item{event}{Censoring indicator, 1 for event, 0 for censored.}
  \item{t0}{Vector (or single number) of time points to compute confidence interval for.}
}
\value{Matrix with values of Kaplan-Meier estimate at \eqn{t_0}.}
\author{
Kaspar Rufibach (maintainer) \cr \email{kaspar.rufibach@roche.com} 
}
\examples{
# use Acute Myelogenous Leukemia survival data contained in package 'survival'
time <- leukemia[, 1]
status <- leukemia[, 2]

tmp <- Surv(time, status) ~ 1
plot(survfit(tmp, conf.type = "none"), mark = "/", col = 1:2)
kaplanMeier_at_t0(time, status, t0 = c(10, 25, 50))
}

\keyword{htest}
\keyword{survival}
