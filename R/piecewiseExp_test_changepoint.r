piecewiseExp_test_changepoint <- function(peMLE, alpha = 0.05){
  
  # test on number of changepoints, as described in Fang & Zheng (2011)
  K <- peMLE$K
  lambda <- peMLE$lambda
  lambda.SE <- peMLE$lambda.SE

  test.mat <- data.frame(matrix(NA, ncol = 8, nrow = K))
  colnames(test.mat) <- c("k", "H_0", "alpha_{k-1}", "chi2quantile", "Xw", "p-value", "reject", "established change point")  
  
  for (k in 2:(K + 1)){
    test.mat[k, "k"]                        <- k
    test.mat[k, "H_0"]                      <- paste("lambda", k - 1, " - lambda", k, sep = "")
    test.mat[k, "alpha_{k-1}"]                  <- alpha * 2 ^ (1 - k)
    test.mat[k, "chi2quantile"]             <- qchisq(p = test.mat[k, "alpha_{k-1}"], df = 1, lower.tail = FALSE) 
    test.mat[k, "Xw"]                       <- (lambda[[k - 1]] - lambda[[k]]) ^ 2 / (lambda.SE[[k - 1]] ^ 2 + lambda.SE[[k]] ^ 2)
    test.mat[k, "p-value"]                  <- pchisq(test.mat[k, "Xw"], df = 1, lower.tail = FALSE)
    test.mat[k, "reject"]                   <- as.numeric(test.mat[k, "p-value"] <= test.mat[k, "alpha_{k-1}"])
    test.mat[k, "established change point"] <- peMLE$tau[k - 1]  
  }
  
  test.mat <- test.mat[-1, ]  
  test.mat[is.na(test.mat[, "reject"]), "reject"] <- 0
  
  if (all(test.mat[, "reject"] == 1) == FALSE){
       first.non.rej <- min((1:nrow(test.mat))[test.mat[, "reject"] == 0])
       test.mat[first.non.rej:nrow(test.mat), "established change point"] <- ""
  }
  
  return(test.mat) 
}

