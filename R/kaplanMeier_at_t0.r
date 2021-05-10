kaplanMeier_at_t0 <- function(time, event, t0){
          
     # compute ci according to the formulas in Huesler and Zimmermann, Chapter 21
     obj <- survfit(Surv(time, event) ~ 1)
     Surv_t <- summary(obj)$surv
     t <- summary(obj)$time
     n <- summary(obj)$n.risk
     res <- rep(NA, length(t0))
     
     for (i in 1:length(t0)){
          ti <- t0[i]
          if (min(t) > ti){res[i] <- 1}
          if (min(t) <= ti){    
               if (ti %in% t){res[i] <- Surv_t[t == ti]} else {
                    Surv_ti <- min(Surv_t[t < ti])
                    res[i] <- Surv_ti
               }
          }
     } # end for
     
     res <- cbind(t0, res)
     dimnames(res)[[2]] <- c("t0", "S at t0")
     return(res)
}
