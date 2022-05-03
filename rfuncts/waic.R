### functions for calculating WAIC in the grade of Membership Model ###
require(abind)
# continuous

## notes
# pars is an R x K x 2*p + n array
# each matrix in pars has the first 2*p columns be parameters, everything else is weights
# the first 2*p columns have every odd thing  being the mean, every even thing being a variance

cont_LCM_log_density <- function(x, eta, mu, sigma2, K){
  d <- c()
  for (h in 1:K){
    d[h] <- eta[h] * prod(dnorm(x, mean = mu[h, ], sd = sqrt(sigma2[h,])))
  }
  d <- log(sum(d))
  return(d)
}

cont_LCM_density <- function(x, eta, mu, sigma2, K){
  d <- c()
  for (h in 1:K){
    d[h] <- eta[h] * prod(dnorm(x, mean = mu[h, ], sd = sqrt(sigma2[h,])))
  }
  d <- sum(d)
  return(d)
}


# logs of small numbers
# remember to ask about what to do for missing values that are imputed

cont_LCM_WAIC <- function(Beta=1, fit, burnin, K) {
  # fit should have theta and eta parameters
  theta <- fit$theta
  eta <- fit$eta[-(1:(burnin+1)), ]
  y_miss <- fit$y_array[-(1:(burnin+1)), , ]
  p <- length(theta)
  n <- dim(y_miss)[2]
  K <- ncol(eta)
  R <- nrow(eta)
  theta <- lapply(theta, aperm, perm = c(3, 1, 2))
  # convert to one big set of arrays to make easiest
  theta <- abind(theta[[1]][-(1:(burnin+1)), , ], theta[[2]][-(1:(burnin+1)), , ], theta[[3]][-(1:(burnin+1)), , ], theta[[4]][-(1:(burnin+1)), , ],
                 theta[[5]][-(1:(burnin+1)), , ], theta[[6]][-(1:(burnin+1)), , ], theta[[7]][-(1:(burnin+1)), , ], theta[[8]][-(1:(burnin+1)), , ],
                 theta[[9]][-(1:(burnin+1)), , ], theta[[9]][-(1:(burnin+1)), , ], theta[[10]][-(1:(burnin+1)), , ],
                 along = 3) # note: might be error here due to p
  # storage for the MCMC samples
  l_p_samples <- matrix(0, nrow = R, ncol = n)
  p_samples <- matrix(0, nrow = R, ncol = n)
  
  # computing MCMC estimates
  for (r in 1:R){
    theta_r <- theta[r, , ]
    y_r <- y_miss[r, , ]
    q <- ncol(theta_r)
    mu_r <- theta_r[, -seq(2, q, by = 2)]
    sigma2_r <- theta_r[, seq(2, q, by = 2)]
    l_p_samples[r, ] <- apply(y_r, 1, cont_LCM_log_density, eta = eta[r, ], mu = mu_r, sigma2 = sigma2_r, K = K)
    p_samples[r, ] <-  apply(y_r, 1, cont_LCM_density, eta = eta[r, ], mu = mu_r, sigma2 = sigma2_r,K = K)
  }
  l_p_estimate <- colMeans(l_p_samples)
  l_p_sq_estimate <- colMeans(l_p_samples^2)
  p_star_estimate <- colMeans(p_samples)
  BtLn <- -mean(log(p_star_estimate))
  V_n <- sum(l_p_sq_estimate - l_p_estimate^2)
  
  # return WAIC
  WAIC <- BtLn + (Beta/n)*V_n
  return(WAIC)
}

