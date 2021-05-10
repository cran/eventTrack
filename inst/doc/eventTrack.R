## ----setup, include = FALSE---------------------------------------------------
now <- as.character(as.POSIXlt(Sys.time()))
today <- as.Date(substr(now, 1, 10))
now <- paste(today, " at ", substr(now, 12, 19), sep = "")

## set some knitr options
knitr::opts_chunk$set(echo = TRUE, cache = TRUE)

## ---- include=TRUE, echo=TRUE-------------------------------------------------
citation("eventTrack")

## ---- include=TRUE, echo=TRUE, results="hide", message=F----------------------

# load necessary packages
library(rpact)          # to simulate illustrative dataset
library(eventTrack)     # event tracking methodology
library(fitdistrplus)   # to estimate Weibull fit
library(knitr)          # to prettify some outputs

## ---- include=TRUE, echo=TRUE, fig.width = 8, fig.height = 6------------------

# ------------------------------------------------------------
# trial details 
# ------------------------------------------------------------
alpha <- 0.05
beta <- 0.2

# design
design <- getDesignGroupSequential(informationRates = 1, typeOfDesign = "asOF", 
                                   sided = 1, alpha = alpha, beta = beta)

# max number of patients
maxNumberOfSubjects1 <- 500
maxNumberOfSubjects2 <- 500
maxNumberOfSubjects <- maxNumberOfSubjects1 + maxNumberOfSubjects2

# survival hazard
piecewiseSurvivalTime <- c(0, 6, 9)
lambda2 <- c(0.025, 0.04, 0.02)

# effect size
hr <- 0.75

# quantities that determine timing of reporting events
dout1 <- 0
dout2 <- 0
douttime <- 12
accrualTime <- 0
accrualIntensity <- 42

# maximal sample size
samplesize <- getSampleSizeSurvival(design = design, piecewiseSurvivalTime = piecewiseSurvivalTime, 
                        lambda2 = lambda2, hazardRatio = hr, dropoutRate1 = dout1, 
                        dropoutRate2 = dout2, dropoutTime = douttime, accrualTime = accrualTime, 
                        accrualIntensity = accrualIntensity,
                        maxNumberOfSubjects = maxNumberOfSubjects)  
nevents <- as.vector(ceiling(samplesize$eventsPerStage))

# expected number of events overall and per group at 1-60 months after trial start
time <- 0:60 
probEvent <- getEventProbabilities(time = time, 
    piecewiseSurvivalTime = piecewiseSurvivalTime, lambda2 = lambda2, 
    hazardRatio = hr,
    dropoutRate1 = dout1, dropoutRate2 = dout2, dropoutTime = douttime,
    accrualTime = accrualTime, accrualIntensity = accrualIntensity, 
    maxNumberOfSubjects = maxNumberOfSubjects)

# extract expected number of events
expEvent <- probEvent$overallEventProbabilities * maxNumberOfSubjects   # overall
expEventInterv <- probEvent$eventProbabilities1 * maxNumberOfSubjects1  # intervention
expEventCtrl <- probEvent$eventProbabilities2 * maxNumberOfSubjects2   # control

# cumulative number of recruited patients over time
nSubjects <- getNumberOfSubjects(time = time, accrualTime = accrualTime, 
                                 accrualIntensity = accrualIntensity, 
                                 maxNumberOfSubjects = maxNumberOfSubjects)

# plot results
plot(time, nSubjects$numberOfSubjects, type = "l", lwd = 2,
     main = "Number of patients and expected number of events over time",
     xlab = "Time from trial start (months)",
     ylab = "Number of patients / events")
lines(time, expEvent, col = "blue", lwd = 2)
lines(time, expEventCtrl, col = "green")
lines(time, expEventInterv, col = "orange")
legend(x = 23, y = 850,
       legend = c("Recruitment", "Expected events - All subjects",
                "Expected events - Control", "Expected events - Intervention"),
       col = c("black", "blue", "green", "orange"), lty = 1, lwd = c(2, 2, 1, 1),
       bty = "n")

# add lines with event prediction pre-trial
fa_pred_pretrial <- exactDatesFromMonths(data.frame(1:length(expEvent) - 1, expEvent), nevents)
segments(fa_pred_pretrial, 0, fa_pred_pretrial, nevents, lwd = 2, lty = 2, col = "blue")
segments(0, nevents, fa_pred_pretrial, nevents, lwd = 2, lty = 2, col = "blue")

## ---- include=TRUE, echo=TRUE-------------------------------------------------
seed <- 2020

# ------------------------------------------------------------
# generate a dataset with 1000 patients
# administratively censor after 100 events
# ------------------------------------------------------------

# number of events after which to administratively censor
nevents_admincens <- 100

# generate a dataset after the interim analysis, administratively censor after 100 events
trial1 <- getSimulationSurvival(design,
    piecewiseSurvivalTime = piecewiseSurvivalTime, lambda2 = lambda2,
    hazardRatio = hr, dropoutRate1 = dout1, dropoutRate2 = dout2, dropoutTime = douttime,
    accrualTime = accrualTime, accrualIntensity = accrualIntensity,
    plannedEvents = nevents_admincens, directionUpper = FALSE,
    maxNumberOfSubjects = maxNumberOfSubjects, maxNumberOfIterations = 1,
    maxNumberOfRawDatasetsPerStage = 1, seed = seed)
trial2 <- getRawData(trial1)

# clinical cutoff date for this analysis
ccod <- trial2$observationTime[1]

## ---- include=TRUE, echo=TRUE-------------------------------------------------
kable(head(trial2, 10)[, 3:ncol(trial2)], row.names = FALSE)

## ---- include=TRUE, echo=TRUE-------------------------------------------------

# ------------------------------------------------------------
# generate a dataset with PFS information:
# - blinded to treatment assignment
# - remove patients recruited *after* clinical cutoff date
# ------------------------------------------------------------
pfsdata <- subset(trial2, subset = (timeUnderObservation > 0), 
                  select = c("subjectId", "timeUnderObservation", "event"))
colnames(pfsdata) <- c("pat", "pfs", "event")
pfsdata[, "event"] <- as.numeric(pfsdata[, "event"])

# define variables with PFS time and censoring indicator
pfstime <- pfsdata[, "pfs"]
pfsevent <- pfsdata[, "event"]

## ---- include=TRUE, echo=TRUE, fig.width = 8, fig.height = 6------------------

# --------------------------------------------------
# analyze hazard functions
# --------------------------------------------------
oldpar <- par(las = 0, mar = c(5.1, 4.1, 4.1, 2.1))
par(las = 1, mar = c(4.5, 5.5, 3, 1))

# kernel estimate of hazard function
h.kernel1 <- muhaz(pfstime, pfsevent)
plot(h.kernel1, xlab = "time (months)", ylab = "", col = grey(0.5), lwd = 2, 
     main = "hazard function", xlim = c(0, max(pfstime)), ylim = c(0, 0.06))
par(las = 3)
mtext("estimates of hazard function", side = 2, line = 4)
rug(pfstime[pfsevent == 1])

# MLE of exponential hazard, i.e. no changepoints
exp_lambda <- sum(pfsevent) / sum(pfstime)
segments(0, exp_lambda, 20, exp_lambda, lwd = 1, col = 2, lty = 1)

# add estimated piecewise exponential hazards
select <- c(1, 3, 5)
legend("topleft", legend = c("kernel estimate of hazard", "exponential hazard", 
                             paste("#changepoints = ", 
     select, sep = ""), "true hazard"), col = c(grey(0.5), 2, 3:(length(select) + 2), 1), 
     lty = c(1, 1, rep(2, length(select)), 1), lwd = c(2, 1, rep(2, length(select)), 2), 
     bty = "n")
for (i in 1:length(select)){
  pe1 <- piecewiseExp_MLE(pfstime, pfsevent, K = select[i])
  tau.i <- c(0, pe1$tau, max(pfstime))
  lambda.i <- c(pe1$lambda, tail(pe1$lambda, 1))  
  lines(tau.i, lambda.i, type = "s", col = i + 2, lwd = 2, lty = 2)
  segments(max(tau.i), tail(lambda.i, 1), 10 * max(pfstime), 
           tail(lambda.i, 1), col = i + 2, lwd = 2)
}

# since we simulated data we can add the truth
lines(stepfun(piecewiseSurvivalTime[-1], lambda2), lty = 1, col = 1, lwd = 2)

## ---- include=TRUE, echo=TRUE-------------------------------------------------

# --------------------------------------------------
# starting from a piecewise Exponential hazard with 
# K = 5 change points, infer the last "significant"
# change point
# --------------------------------------------------
pe5 <- piecewiseExp_MLE(time = pfstime, event = pfsevent, K = 5)
pe5.tab <- piecewiseExp_test_changepoint(peMLE = pe5, alpha = 0.05)
cp.select <- max(c(0, as.numeric(pe5.tab[, "established change point"])), na.rm = TRUE)
kable(subset(pe5.tab, select = c("k", "H_0", "alpha_{k-1}", "p-value", "reject", 
                                 "established change point")), row.names = FALSE)

## ---- include=TRUE, echo=TRUE, fig.width = 8, fig.height = 6------------------

# ------------------------------------------------------------
# projection for patients to be recruited after CCOD
# 42 patients / month, and we have maxNumberOfSubjects - length(pfstime)
# left to randomize. Assume they show up uniformly distributed
# this accrual is relative to the CCOD
# ------------------------------------------------------------
accrual_after_ccod <- 1:(maxNumberOfSubjects - length(pfstime)) / accrualIntensity 

# Kaplan-Meier and hybrid Exponential tail fits, for different changepoints
par(las = 1, mar = c(4.5, 4.5, 3, 1))
plot(survfit(Surv(pfstime, pfsevent) ~ 1), mark = "", xlim = c(0, 30), 
     ylim = c(0, 1), conf.int = FALSE, xaxs = "i", yaxs = "i",
     main = "estimated survival functions", xlab = "time", 
     ylab = "survival probability", col = grey(0.75), lwd = 5)

# how far out should we predict monthly number of events?
future.units <- 15
tab <- matrix(NA, ncol = 2, nrow = future.units)
ts <- seq(0, 100, by = 0.01)

# the function predictEvents takes as an argument any survival function
# hybrid exponential with cp.select
St1 <- function(t0, time, event, cp){
     return(hybrid_Exponential(t0 = t0, time = time, event = event, 
     changepoint = cp))}

pe1 <- predictEvents(time = pfstime, event = pfsevent, 
          St = function(t0){St1(t0, time = pfstime, 
          event = pfsevent, cp.select)}, accrual_after_ccod, 
          future.units = future.units)
tab[, 1] <- pe1[, 2]
lines(ts, St1(ts, pfstime, pfsevent, cp.select), col = 2, lwd = 2)

# for illustration add piecewise exponential with five changepoints
St2 <- function(t0, time, event, tau, lambda){
      return(1 - getPiecewiseExponentialDistribution(time = t0, 
                        piecewiseSurvivalTime = tau, 
                        piecewiseLambda = lambda))}

pe2 <- predictEvents(time = pfstime, event = pfsevent, 
          St = function(t0){St2(t0, time = pfstime, 
          event = pfsevent, tau = c(0, pe5$tau), 
          lambda = pe5$lambda)}, accrual_after_ccod, 
          future.units = future.units)
tab[, 2] <- pe2[, 2]
lines(ts, St2(ts, pfstime, pfsevent, c(0, pe5$tau), pe5$lambda), col = 3, lwd = 2)

# since we have simulated data we can add the truth
tp <- seq(0, 100, by = 0.01) 
lines(tp, 1 - getPiecewiseExponentialDistribution(time = tp, 
                        piecewiseSurvivalTime = piecewiseSurvivalTime, 
                        piecewiseLambda = lambda2), col = 4, lwd = 2, lty = 1)

# legend
legS <- c("Kaplan-Meier", "hybrid exponential", 
          "piecewise exponential K = 5", "true S")
legend("bottomleft", legS, col = c(grey(0.75), 2:4), 
       lwd = c(5, 2, 2, 2), bty = "n", lty = 1)

# prettify the output
tab <- data.frame(pe1[, 1], round(tab, 1))
colnames(tab) <- c("month", legS[2:3])

## ---- include=TRUE, echo=TRUE-------------------------------------------------
kable(tab, row.names = FALSE)

## ---- include=TRUE, echo=TRUE-------------------------------------------------
fa_pred_ccod <- exactDatesFromMonths(tab[, c("month", "hybrid exponential")], nevents)
fa_pred_ccod

## ---- include=TRUE, echo=TRUE-------------------------------------------------
fa_pred_ccod + ccod

## ---- include=TRUE, echo=TRUE, fig.width = 8, fig.height = 6------------------

# ------------------------------------------------------------
# plot event predictions pre-trial and after CCOD
# ------------------------------------------------------------
par(las = 1, mar = c(4.5, 4.5, 3, 1))
plot(time, nSubjects$numberOfSubjects, type = "l", lwd = 5, col = grey(0.75), 
     main = "Number of patients and expected number of events over time",
     xlab = "Time from trial start (months)",
     ylab = "Number of patients / events")

# expected events at design stage
lines(time, expEvent, col = 4, lwd = 3, lty = 2)

# patients already recruited
tot_accrual1 <- subset(trial2, select = "accrualTime", subset = timeUnderObservation > 0)[, 1]

# add those that will only be recruited after CCOD
tot_accrual <- c(tot_accrual1, accrual_after_ccod + ccod)
tot_accrual_unit <- cumsum(table(cut(tot_accrual, 0:59, labels = FALSE)))

# actually not much to see, b/c simulated data perfectly uniform
lines(1:length(tot_accrual_unit), tot_accrual_unit, col = 1) 

# now add event prediction based on hybrid exponential and 100 events
lines(tab$month + ccod, tab[, "hybrid exponential"], col = 2, lwd = 3)

# add lines with prediction based on CCOD data
segments(fa_pred_ccod + ccod, 0, fa_pred_ccod + ccod, nevents, lwd = 3, col = 2, lty = 3)
segments(0, nevents, fa_pred_ccod + ccod, nevents, lwd = 3, col = 2, lty = 3)

# finally, we add observed events over time
# note the x-axis of the plot (time from trial start), so we need to add accrual time
events_trial_time <- sort((pfstime + tot_accrual1)[pfsevent == 1])
lines(stepfun(events_trial_time, 0:sum(pfsevent)), pch = NA, col = "orange", lwd = 2)

# legend
legend(x = 20, y = 930, c("Planned recruitment at design stage", 
                          "Recruitment (up to CCOD + prediction)", 
                          "Expected events at design stage", 
                          "Observed events", 
                          "Expected events (after CCOD prediction)"),
       col = c(grey(0.75), 1, 4, "orange", 2), lty = c(1, 1, 2, 1, 1), 
       lwd = c(3, 1, 3, 2, 3), bty = "n")

# revert back to initial plot parameters
par(oldpar)

## ---- include=TRUE, echo=TRUE-------------------------------------------------

# ------------------------------------------------------------
# compute a confidence interval for our event prediction using bootstrap
# ------------------------------------------------------------

# do not run as takes a few minutes, code only here for illustration
if (1 == 0){
boot1 <- bootstrapTimeToEvent(time = pfstime, event = pfsevent, 
          interim.gates = nevents, future.units = 20, n = length(pfstime), 
          M0 = 1000, K0 = 5, alpha = 0.05, accrual = accrual_after_ccod, seed = 2020)
ci <- quantile(boot1$event.dates[, 1], c(0.025, 0.975))
}

# hard-code result and add ccod
ci <- c(10.26989, 12.50578) + ccod
ci

## ---- include=TRUE, echo=TRUE-------------------------------------------------

# ------------------------------------------------------------
# event prediction based on pure Weibull assumption
# ------------------------------------------------------------

# prepare dataset to fit function "fitdistcens"
pfsdata2 <- data.frame(cbind("left" = pfstime, "right" = pfstime))
pfsdata2[pfsevent == 0, "right"] <- NA

# estimate parameters based on Weibull distribution
fit_weibull <- suppressWarnings(fitdistcens(pfsdata2, "weibull"))

# event prediction based on this estimated survival function
# simply specify survival function based on above fit
pred_weibull <- predictEvents(time = pfstime, event = pfsevent, 
                      St = function(t0){1 - pweibull(t0, 
                                                     shape = fit_weibull$estimate["shape"], 
                                                     scale = fit_weibull$estimate["scale"])}, 
                      accrual = accrual_after_ccod, future.units = future.units)
ccod_weibull <- exactDatesFromMonths(pred_weibull, nevents)
ccod_weibull
ccod_weibull + ccod

## ---- include=TRUE, echo=TRUE-------------------------------------------------

# ------------------------------------------------------------
# input data (e.g. OS or PFS)
# we use the data simulated above
# alternatively, simply use
# pfsdata <- read.csv("C:/my path/pfs.csv")
# --------------------------------------------------

# accrual after ccod
accrual_after_ccod <- 1:(maxNumberOfSubjects - length(pfstime)) / accrualIntensity 

# infer changepoint
pe5 <- piecewiseExp_MLE(time = pfstime, event = pfsevent, K = 5)
pe5.tab <- piecewiseExp_test_changepoint(peMLE = pe5, alpha = 0.05)
cp.selected <- suppressWarnings(max(as.numeric(pe5.tab[, "established change point"]), 
                                             na.rm = TRUE))
cp.selected <- max(cp.selected, 0)

# compute predictions for selected hybrid exponential model, i.e. chosen changepoint
future.units <- 15
St1 <- function(t0, time, event, cp){
     return(hybrid_Exponential(t0 = t0, time = time, event = event, 
     changepoint = cp))}

pe1 <- predictEvents(time = pfstime, event = pfsevent, 
          St = function(t0){St1(t0, time = pfstime, 
          event = pfsevent, cp.selected)}, accrual_after_ccod, 
          future.units = future.units)

# find exact timepoint when targeted number of events are reached through linear interpolation
# add computation of confidence interval using code from just above
fa_pred_ccod <- exactDatesFromMonths(tab[, c("month", "hybrid exponential")], nevents)
fa_pred_ccod + ccod

## ---- include=TRUE, echo=TRUE-------------------------------------------------
library(gestate)

# accrual: observed + predicted (observed accrual not needed when using eventTrack!)
tot_accrual1 <- subset(trial2, select = "accrualTime", subset = timeUnderObservation > 0)[, 1]
tot_accrual <- c(tot_accrual1, accrual_after_ccod + ccod)
tot_accrual_unit <- cumsum(table(cut(tot_accrual, 0:59, labels = FALSE)))

# Create an RCurve incorporating both observed and predicted recruitment 
recruit <- PieceR(matrix(c(rep(1, length(tot_accrual_unit)), 
                           c(tot_accrual_unit[1], diff(tot_accrual_unit))), ncol = 2), 1)

# do the prediction: need to supply CCOD (ccod), number of events at CCOD (nevents_admincens), 
# and numbers at risk then (length(tot_accrual1))
predictions <- event_prediction(data = pfsdata, Time = "pfs", Event = "event", 
                                censoringOne = FALSE, rcurve = recruit, max_time = 60, 
                                cond_Events = nevents_admincens, 
                                cond_NatRisk = length(tot_accrual1), 
                                cond_Time = ceiling(ccod), units = "Months")

# exact date (not exactly as Weibull fit above, as cond_Time argument needs to be an integer)
pred_gestate <- exactDatesFromMonths(predictions$Summary[, c("Time", "Predicted_Events")], nevents)
pred_gestate

