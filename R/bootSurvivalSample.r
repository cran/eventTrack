bootSurvivalSample <- function(surv.obj, n, M = 1000){
     
     # bootstrap right-censored survival data according to Efron (1981)
     # see Akritas (1986) for details
     # surv.obj: Surv object to sample from
     # n: sample size for samples to be drawn
     # M: number of samples
     
     n.surv <- nrow(surv.obj)
     res.mat <- data.frame(matrix(NA, ncol = M, nrow = n))
     
     for (i in 1:M){
          ind <- sample(1:n.surv, size = n, replace = TRUE)
          res.mat[, i] <- surv.obj[ind]
          print(paste("sample ", i, " of ", M, " done", sep = ""))
     }
     
     return(res.mat)
}

