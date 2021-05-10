predictEventsUncond <- function(St, accrual, future.units = 50){

     # --------------------------------------------------
     # compute expected number of events according to Fang & Su (2011)
     # unconditional method
     # --------------------------------------------------
     t.predict <- seq(0, length = future.units, by = 1) + 1
     predict.events <- rep(NA, length(t.predict))
     for (i in 1:length(predict.events)){
     
          # rescale accrual to match timescale, in months.
          accrual.i <- t.predict[i] - accrual
          predict.events[i] <- sum(1 - St(t0 = accrual.i[accrual.i >= 0]))
     }

     predev <- data.frame(t.predict, predict.events)
     return(predev)    
}
