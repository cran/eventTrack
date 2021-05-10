\name{predictEvents}
\alias{predictEvents}
\title{Compute timepoint when a certain number of events in a time-to-event study is reached}
\description{Based on a specified survival function and potential future accrual, compute for each month in the future the number of expected events reached by then.}
\usage{predictEvents(time, event, St, accrual, future.units = 50)}
\arguments{
  \item{time}{Event times, censored or observed, in months.}
  \item{event}{Censoring indicator, 1 for event, 0 for censored.}
  \item{St}{Function that specifies the survival function to be used.}
  \item{accrual}{\code{NULL} if trial is fully accrued. For potential future accural, see example in \code{\link{hybrid_Exponential}} of how to specify it.}
  \item{future.units}{Number of future months to compute prediction for.}
}
\value{A \code{data.frame} with the months and corresponding expected events.}

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
