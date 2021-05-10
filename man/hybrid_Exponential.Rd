\name{hybrid_Exponential}
\alias{hybrid_Exponential}
\title{Estimate survival function, as hybrid between Kaplan-Meier and Exponential tail}
\description{This function estimates the values of the survival function as a hybrid of Kaplan-Meier for times smaller than the specified change point and an Exponential fit to the tail of the survival function. The Exponential tail fit is computed assuming a piecewise constant hazard with one change point.}
\usage{hybrid_Exponential(t0, time = time, event = event, changepoint)}
\arguments{
  \item{t0}{Value at which to compute value of survival function. Can be a vector.}
  \item{time}{Event times, censored or observed, in months.}
  \item{event}{Censoring indicator, 1 for event, 0 for censored.}
  \item{changepoint}{Pre-specified change point.}
}
\value{A vector of the same dimension as \code{t0} containing the values of the estimated survival function at \code{t0}.}

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

Rufibach, K. (2016).
\emph{Event projection: quantify uncertainty and manage expectations of broader teams.}
Slides for talk given in Basel Biometric Section Seminar on 28th April 2016.
\url{http://bbs.ceb-institute.org/wp-content/uploads/2016/06/Kaspar-event_tracking.pdf}.
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

\keyword{htest}
