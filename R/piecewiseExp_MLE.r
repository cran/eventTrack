piecewiseExp_MLE <- function(time, event, K){
  
  # estimate piecewise Exponential hazard function
  # with implicit estimation of break points, for number of break points K
  tau0 <- seq(max(time) / 10, max(time) * 0.7, length.out = K)

  suppressWarnings({
    if (K == 1){tau <- optimize(f = piecewiseExp_profile_loglik_tau, interval = range(time), 
                              time = time, event = event, maximum = TRUE)$maximum}
    if (K > 1){tau <- optim(par = tau0, fn = piecewiseExp_profile_loglik_tau, time = time, 
                          event = event, control = list(fnscale = -1))$par}}
  )

  # compute Hessian matrix and standard errors for lambdas
  results <- lambda_j_Exp(tau, time, event)
  lambda.SE <- sqrt(-1 / results$hess.diag)
  
  res <- list("tau" = tau, "lambda" = results$aj, "lambda.SE" = lambda.SE, "K" = K)
  return(res)
}
