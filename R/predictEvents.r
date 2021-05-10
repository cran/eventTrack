predictEvents <- function(time, event, St, accrual, future.units = 50){

     # --------------------------------------------------
     # compute expected number of events according to Fang & Su (2011)
     # Hybrid approach: Kaplan-Meier until last change-point, 
     # from there Exponential fit to tail
     # --------------------------------------------------
     t.predict <- seq(0, length = future.units, by = 1) + 1

     # events that hapened so far:
     n <- sum(event)

     # indicator for patients still censored:
     cens.ind <- (event == 0)
     predict.events <- rep(NA, length(t.predict))
     for (i in 1:length(predict.events)){
     
          # events for patients already accrued
          t0 <- t.predict[i]
     
          Sti <- St(t0 = time[cens.ind])
          StiT <- St(t0 = time[cens.ind] + t0)
          e1 <- sum((Sti - StiT) / Sti)
     
          # events for those patients still to be accrued. Rescale accrual to match 
          # timescale, in months.
          accrual.i <- t.predict[i] - accrual
          accrual.i <- accrual.i[accrual.i >= 0]
          e2 <- 0
          if (length(accrual.i) > 0){e2 <- sum(1 - St(t0 = accrual.i))}
          m <- e1 + e2
          predict.events[i] <- n + m
     }

     predev <- data.frame(t.predict, predict.events)
     return(predev)    
}
