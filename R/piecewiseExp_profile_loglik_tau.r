piecewiseExp_profile_loglik_tau <- function(tau, time, event){
  
  # Maximum likelihood function to compute lambdas in piecewise Exponential
  # time-to-event model. see Fang & Su (2011) for details
  
  k <- length(tau)
  tau <- c(0, tau, Inf)
  loglik <- 0
  
  for (j in 2:(k + 2)){
    Xtj <- sum(event[time <= tau[j]])
    Xtj1 <- sum(event[time <= tau[j - 1]])
    diff <- Xtj - Xtj1
    nenner <- sum((pmin(time, tau[j]) - tau[j - 1])[time > tau[j - 1]])
    aj <- diff / nenner
    if (diff != 0){loglik <- loglik + diff * log(aj)}
  }
  
  return(loglik) 
}
