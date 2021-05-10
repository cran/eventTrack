bootstrapTimeToEvent <- function(time, event, interim.gates, future.units, n, M0 = 1000, K0 = 5, alpha = 0.05, 
                                 accrual, seed = 2014){
  
  
  # --------------------------------------------------------------
  # bootstrap estimated time to event
  # same sample size of bootstrapped samples as original sample size
  # --------------------------------------------------------------

     set.seed(seed)

     # generate bootstrap samples
     bootstrappedSamples <- bootSurvivalSample(surv.obj = Surv(time, event), n = n, M = M0)

     # collect piecewise Exponential estimates
     pe.MLEs <- NULL

     # collect corresponding changepoint tables
     pe.tabs <- NULL

     # collect chosen change points
     changepoints <- rep(NA, M0)

     # collect estimated survival functions
     t0s <- seq(0, 100, by = 0.1)
     estS <- matrix(NA, ncol = M0, nrow = length(t0s))
     rownames(estS) <- t0s

     # collect dates when we see the exact number of events
     event.dates <- data.frame(matrix(NA, nrow = M0, ncol = length(interim.gates)))

     for (i in 1:M0){
          boot.i <- as.matrix(bootstrappedSamples[, i])
          pe.i <- piecewiseExp_MLE(time = boot.i[, 1], event = boot.i[, 2], K = K0)
          pe.MLEs[[i]] <- pe.i
     
          # only do inference in case MLE of lambda has no entry of 0
          pe.i.tab <- piecewiseExp_test_changepoint(peMLE = pe.i, alpha = alpha)
          pe.tabs[[i]] <- pe.i.tab
          
          established.changepoint <- 0
          if (all(is.na(as.numeric(pe.i.tab[, "established change point"]))) == FALSE){
               established.changepoint <- max(as.numeric(pe.i.tab[, "established change point"]), na.rm = TRUE)
          }
          changepoints[i] <- established.changepoint
     
          St1 <- function(t0, time, event, cp){
               return(hybrid_Exponential(t0 = t0, time = time, event = event, 
                               changepoint = cp))}
          estS[, i] <- St1(t0s, time = boot.i[, 1], event = boot.i[, 2], cp = established.changepoint)
     
          pred.i <- predictEvents(time = boot.i[, 1], event = boot.i[, 2], 
                      St = function(t0){St1(t0, time, event, established.changepoint)}, accrual,
                      future.units = future.units)
     
          # day when we have enough events for the planned analysis
          for (j in 1:length(interim.gates)){event.dates[i, j] <- exactDatesFromMonths(predicted = pred.i, interim.gates[j])}
     
          print(paste("date ", i, " of ", M0, " computed", sep = ""))
     }

  # generate output
  res <- list("pe.MLEs" = pe.MLEs, "pe.tabs" = pe.tabs, "changepoints" = changepoints, 
              "estS" = estS, "event.dates" = event.dates)
  return(res)
}




