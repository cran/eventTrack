hybrid_Exponential <- function(t0, time = time, event = event, changepoint){
  
  # Kaplan-Meier with Exponential tail beyond changepoint
  tau1 <- changepoint
  tau2 <- Inf
  
  # Kaplan-Meier prior to changepoint:
  km <- as.numeric(kaplanMeier_at_t0(time = time, event = event, t0 = t0)[, 2])
  res <- rep(NA, length(t0))
  res[t0 <= tau1] <- km[t0 <= tau1]
  
  # estimate of Exponential parameter beyond changepoint
  lambda2 <- lambda_j_Exp(tau = changepoint, time = time, event = event)$aj[2]
  
  # add tail
  Stau <- as.numeric(kaplanMeier_at_t0(time = time, event = event, t0 = tau1)[, 2])
  Sexpo <- exp(-lambda2 * (t0[t0 > tau1] - tau1))
  res[t0 > tau1] <- Stau * Sexpo
  
  return(res)    
}
