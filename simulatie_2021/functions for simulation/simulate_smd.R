simulate_smd <- function(k_train = 20, k_test = 100, mean_n = 40, es = .5,
                         tau2 = 0.04, alpha_tau = 0, alpha_mod = 0, moderators = 5, distribution = "normal",
                         model = "es * x[, 1]")
{
  # Make label for training and test datasets
  training <- c(rep(1, k_train), rep(0, k_test))
  
  # Randomly determine n from a normal distribution with mean
  # mean_n and sd mean_n/3
  n <- rnorm(length(training), mean = mean_n, sd = mean_n/3)
  # Round to even number
  n <- ceiling(n) - ceiling(n)%%2
  # Truncate n so the lowest value can be 8; cases where n=2 will result
  # in errors
  n[n < 8] <- 8
  
  # Generate moderator matrix x:
  if(distribution == "normal") x <- matrix(rnorm(length(n) * moderators), ncol = moderators)
  if(distribution == "sn") x <- matrix(rsn(n = length(n) * moderators, xi = 0, omega = 1, alpha = alpha_mod), ncol = moderators)
  if(distribution == "bernoulli") x <- matrix(rbinom((length(n) * moderators), 1, .5), ncol = moderators)
  if(!(distribution %in% c("normal", "bernoulli", "sn"))) stop(paste0(distribution, "is not a valid distribution for SimulateSMD"))
  
  # Sample true effect sizes theta.i from a normal distribution with mean
  # mu, and variance tau2, where mu is the average
  # population effect size. The value of mu depends on the values of the
  # moderators and the true model mu <- eval(model)
  model <- parse(text = model)
  mu <- eval(model)
  
  # theta.i: true effect size of study i
  theta.i <- mu + rsn(n = length(n), xi = 0, omega = sqrt(tau2), alpha = alpha_tau)
  
  # Then the observed effect size yi is sampled from a non-central
  # t-distribution under the assumption that the treatment group and
  # control group are both the same size
  p_ntk <- 0.5  #Percentage of cases in the treatment group
  ntk <- p_ntk * n  #n in the treatment group for study i
  nck <- (1 - p_ntk) * n  #n in the control group for study i
  df <- n - 2  #degrees of freedom
  j <- 1 - 3/(4 * df - 1)  #correction for bias
  nk <- (ntk * nck)/(ntk + nck)
  ncp <- theta.i * sqrt(nk)  #Non-centrality parameter
  
  # Standardized mean difference drawn from a non-central t-distribution
  SMD <- mapply(FUN = rt, n = 1, df = df, ncp = ncp)
  
  # yi is Hedges' g for study i
  yi <- SMD/((j^-1) * (nk^0.5))
  
  # Calculate the variance of the effect size
  vi <- j^2 * (((ntk + nck)/(ntk * nck)) + ((yi/j)^2/(2 * (ntk + nck))))
  
  # Dersimonian and Laird estimate of tau2
  Wi <- 1/vi[1:k_train]
  tau2_est <- max(0, (sum(Wi * (yi[1:k_train] - (sum(Wi * yi[1:k_train])/sum(Wi)))^2) -
                        (k_train - 1))/(sum(Wi) - (sum(Wi^2)/sum(Wi))))
  
  data <- data.frame(training, vi, yi, x)
  
  list(training = subset(data, training == 1, -1), testing = subset(data,
                                                                    training == 0, -c(1, 2)), housekeeping = data.frame(n = n, mu_i = mu, theta_i = theta.i),
       tau2_est = tau2_est)
}

n_to_vi <- function(n){
  3.84417/n^.93728
}