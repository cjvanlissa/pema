if(FALSE){
  B <- c(.5, .4)
  n <- 100
  p <- length(B)
  residual_sd <- sqrt(sum(B^2)) # Use L2 norm of B, so the R^2 is always .5

  X <- matrix(rnorm(p*n), ncol = p)

  Y <- X%*%B + rnorm(n, sd = residual_sd)

  var(Y)


  regr <- lm(Y~X)

  var(regr$residuals)

  Y2 <- X%*%c(.3, .2)

  resid <- Y-Y2
  var(resid)


  summary(lm(Y~X))



}
