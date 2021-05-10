\name{exactDatesFromMonths}
\alias{exactDatesFromMonths}
\title{Compute exact timepoint when a certain number of events is reached, based on monthly number of events}
\description{Based on monhtly number of events, compute exact day when the pre-specified number of events happens, using simple linear interpolation.}
\usage{exactDatesFromMonths(predicted, nevent)}
\arguments{
  \item{predicted}{\code{data.frame}, with first column monthly dates and second column the number of events in each month.}
  \item{nevent}{Number of targeted events.}
}
\value{The exact date.}

\author{
Kaspar Rufibach (maintainer) \cr \email{kaspar.rufibach@roche.com} 
}

\keyword{htest}
