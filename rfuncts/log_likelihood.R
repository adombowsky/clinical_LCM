log_lik <- function(y, mu, sigma2, eta) {
  # function to compute log-likelihood
  # use waic functions
  source("rfuncts/waic.R")
  l_l <- sum(apply(y, 1, cont_LCM_log_density, eta = eta, mu = mu, sigma2 = sigma2, K = length(eta)))
  return(l_l)
}