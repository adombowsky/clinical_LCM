# in which I tune a prior for CD4%
require(ggplot2)
require(dplyr)
# want CD4% to be clustered around 40, try not to go above 50

# first, write sample NI-G function
sample_noig_prior <- function(n, zeta, lambda, psi, nu){
  # sample variance
  sigma2 <- 1/rgamma(n = n, shape = psi, rate = nu)
  mu <- rnorm(n = n, mean = zeta, sd = sqrt(sigma2/lambda))
  return(list(mu = mu, sigma2 = sigma2))
}

# and write sample prior predict function
sample_noig_predict <- function(n, zeta, lambda, psi, nu){
  # sample variance
  sigma2 <- 1/rgamma(n = 1, shape = psi, rate = nu)
  mu <- rnorm(n = 1, mean = zeta, sd = sqrt(sigma2/lambda))
  x <- rnorm(n = n, mean = mu, sd = sqrt(sigma2))
  return(x)
}

### sampling from prior predictive
pps <- sample_noig_prior(n = 10000, zeta = 32, lambda = 10, psi = 2, nu = 50)
x <- sample_noig_predict(n = 10000, zeta = 32, lambda = 10, psi = 2, nu = 50)

### based off stats in here: https://journals.sagepub.com/doi/10.1177/2325957414530472