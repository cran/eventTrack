lambda_j_Exp <- function(tau, time, event){

  k <- length(tau)
  aj <- rep(NA, k + 1)
  tau <- c(0, tau, Inf)
  hess.diag <- aj
  
  # compute lambdas and Hessian matrix accompanying lambdas 
  # (used to compute SEs)
  for (j in 2:(k + 2)){  
    
    # lambdas
    Xtj <- sum(event[time <= tau[j]])
    Xtj1 <- sum(event[time <= tau[j - 1]])
    diff <- Xtj - Xtj1
    nenner <- sum((pmin(time, tau[j]) - tau[j - 1]) * (time > tau[j - 1]))
    aj[j - 1] <- diff / nenner
    
    # diagonal for Hessian matrix for lambdas
    hess.diag[j - 1] <- - diff / aj[j - 1] ^ 2
  }
  
  res <- list("aj" = aj, "hess.diag" = hess.diag)
  return(res)
}
