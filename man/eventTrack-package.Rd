\name{eventTrack-package}
\alias{eventTrack-package}
\alias{eventTrack}
\docType{package}
\title{Event Prediction for Time-to-Event Endpoints}
\description{Implements the hybrid framework for event prediction described in Fang & Zheng (2011). To estimate the survival function the event prediction is based on, a piecewise Exponential hazard function is fit to the time-to-event data to infer the potential change points. Prior to the last identified change point, the survival function is estimated using Kaplan-Meier, and the tail after the change point is fit using piecewise Exponential. The Weibull version described in Fang and Zheng (2011) is not implemented here.

An example has been presented in this talk: \url{https://baselbiometrics.github.io/home/docs/talks/20160428/1_Rufibach.pdf}.}
\details{
\tabular{ll}{
Package:  \tab eventTrack  \cr
Type:     \tab Package    \cr
Version:  \tab 1.0.3      \cr
Date:     \tab 2023-08-03 \cr
License:  \tab GPL (>=2)  \cr
LazyLoad: \tab yes        \cr
}
}

\section{Validation status}{The functions in this package have been written by Kaspar Rufibach and are generally to be considered experimental.} 

\author{
Kaspar Rufibach (maintainer) \cr 
\email{kaspar.rufibach@roche.com}
}

\examples{

# --------------------------------------------------
# simulate data 
# --------------------------------------------------
set.seed(2021)
n <- 600
time0 <- rexp(n, rate = log(2) / 20)
cens <- rexp(n, rate = log(2) / 50)
time <- pmin(time0, cens)
event <- as.numeric(time0 < cens)
accrual_after_ccod <- 1:(n - length(time)) / 30

# --------------------------------------------------
# compute hybrid estimate and predict timepoint
# --------------------------------------------------
plot(survfit(Surv(time, event) ~ 1), mark = "", xlim = c(0, 200), 
     ylim = c(0, 1), conf.int = FALSE, xaxs = "i", yaxs = "i",
     main = "estimated survival functions", xlab = "time", 
     ylab = "survival probability", col = grey(0.75), lwd = 5)

# how far out should we predict monthly number of events?
future.units <- 15
tab <- matrix(NA, ncol = 2, nrow = future.units)
tab[, 1] <- 1:nrow(tab)
ts <- seq(0, 100, by = 0.01)

# --------------------------------------------------
# starting from a piecewise Exponential hazard with 
# K = 5 change points, infer the last "significant"
# change point
# --------------------------------------------------
pe5 <- piecewiseExp_MLE(time = time, event = event, K = 5)
pe5.tab <- piecewiseExp_test_changepoint(peMLE = pe5, alpha = 0.05)
cp.select <- max(c(0, as.numeric(pe5.tab[, "established change point"])), na.rm = TRUE)

# the function predictEvents takes as an argument any survival function
# hybrid exponential with cp.select
St1 <- function(t0, time, event, cp){
     return(hybrid_Exponential(t0 = t0, time = time, event = event, 
     changepoint = cp))}

pe1 <- predictEvents(time = time, event = event, 
          St = function(t0){St1(t0, time = time, 
          event = event, cp.select)}, accrual_after_ccod, 
          future.units = future.units)
tab[, 2] <- pe1[, 2]
lines(ts, St1(ts, time, event, cp.select), col = 2, lwd = 2)

# --------------------------------------------------
# compute exact date when we see targeted number of events
# for hybrid Exponential model, through linear interpolation
# --------------------------------------------------
exactDatesFromMonths(predicted = tab, 450)
}
