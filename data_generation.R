# This is the file storing functions for data generation
# and relevant plots.

gen_binorm_data <- function(n, mu, Sigma, Th) {
  # n: number of samples
  # mu_{1} for tox, mu_2 for eff
  # Th: threshold to decide toxicity probability
  dat <- rmvnorm(n, mean = mu, sigma = Sigma)
  dat[, 1] <- (dat[, 1] <= Th)
  return(list(
    tox = dat[, 1],
    eff = dat[, 2]
  ))
}

logit <- function(x) {
  log(x / (1 - x))
}

logistic <- function(x) {
  1 / (1 + exp(-x))
}

# according to SE_k and ST_k, calculating mu_E'
par_for_trt <- function(muc, sek, stk, sigmaT, Th) {
  # muc = (mu_{1ck}, mu_{2ck})
  # sek, stk are CRDs for sub-trial k
  # sigmaT = (sigmaTc, sigmaTe) are std of ctrl and trt of tox
  # Th: threshold for tox
  
  ## for toxicity
  pck <- pnorm(Th, mean = muc[1], sd = sigmaT[1])
  pek <- logistic(logit(pck) + stk)
  mu1ek <- uniroot(f = function(x) {
    pnorm(Th, mean = x, sd = sigmaT[2]) - pek
  }, interval = c(-1e3, 1e3))$root
  
  ## for efficacy
  mu2ek <- muc[2] + sek
  
  return(list(
    mu = c(mu1ek, mu2ek),
    pt = c(pck, pek)
  ))
}

st_detect <- function(mu1ck, Th, sigma1ck, pup = 0.45, plow = 0.05) {
  # mu1ck: the mean value when generating original normal toxicity data before transformation
  # sigma1ck: the std of normal tox data (presumed to be fixed)
  # pup and plow: the upper and lower bound of toxicity probabilities
  
  ## get tox prob of the ctrl group
  tpck <- pnorm(Th, mean = mu1ck, sd = sigma1ck)
  
  STmax <- logit(pup) - logit(tpck)
  STmin <- logit(plow) - logit(tpck)
  
  return(c(STmin, STmax))
}

# determine the range of mu1ck given the p_up and p_low
decide_mu1ck <- function(pup, plow, Th, sigma1ck) {
  int_low <- uniroot(f = function(x) pnorm(Th, mean = x, sd = sigma1ck) - plow,
                     interval = c(-1e3, 1e3))$root
  int_up <- uniroot(f = function(x) pnorm(Th, mean = x, sd = sigma1ck) - pup,
                    interval = c(-1e3, 1e3))$root
  return(c(int_low, int_up))
}

gen_data <- function(ntrial, trial_size, muct, muce, st_all, se_all, Th, fcov) {
  # ntrial: number of sub-trials
  # trial_size: sample sizes in each ctrl and trt group (assuming to be the same)
  # muct, muce: all mean values for ctrl tox (eff) data (w.r.t each sub-trial)
  # st_all, se_all: all treatment effects
  # Th: threshold for transforming tox data
  # fcov: fixed covariance matrix of the bi-normal distr
  
  ## a final [max(trial_size) + 1, 2, n_mod] array for all sub-trials
  ## with [max(trial_size) + 1, , ] for toxicity observations
  ## and [, 1, ] for ctrl, [, 2, ] for trt groups
  fr <- array(, dim = c(max(trial_size) + 1, 2, ntrial))
  
  ## for control groups
  ctrl_dat_tox <- NULL; ctrl_dat_eff <- NULL
  for (i in 1:ntrial) {
    ## temporary observations
    temp_dat <- gen_binorm_data(trial_size[i], mu = c(muct[i], muce[i]),
                                Sigma = fcov, Th = Th)
    ctrl_dat_tox <- qpcR:::cbind.na(ctrl_dat_tox, temp_dat$tox)
    ctrl_dat_eff <- qpcR:::cbind.na(ctrl_dat_eff, temp_dat$eff)
  }
  ## processing
  ctrl_dat_tox <- colSums(ctrl_dat_tox[, -1], na.rm = T) # leaving out NA values
  ctrl_dat_eff <- ctrl_dat_eff[, -1]
  
  ## for treatment groups
  trt_dat_tox <- NULL
  trt_dat_eff <- NULL
  for (i in 1:ntrial) {
    ## obtain mean values of trt group data (bi-norm)
    muek <- par_for_trt(muc = c(muct[i], muce[i]), sek = se_all[i], stk = st_all[i],
                        sigmaT = sqrt(diag(fcov)), Th = Th)$mu
    
    temp_dat <- gen_binorm_data(trial_size[i], mu = muek, Sigma = fcov, Th = Th)
    trt_dat_tox <- qpcR:::cbind.na(trt_dat_tox, temp_dat$tox)
    trt_dat_eff <- qpcR:::cbind.na(trt_dat_eff, temp_dat$eff)
  }
  ## processing
  trt_dat_tox <- colSums(trt_dat_tox[, -1], na.rm = T)
  trt_dat_eff <- trt_dat_eff[, -1]
  
  ## combine all data into the final array
  fr[max(trial_size) + 1, 1, ] <- ctrl_dat_tox
  fr[max(trial_size) + 1, 2, ] <- trt_dat_tox
  fr[1:max(trial_size), 1, ] <- ctrl_dat_eff
  fr[1:max(trial_size), 2, ] <- trt_dat_eff
  
  return(fr)
}

out_corr <- function(mu, Sigma, Th) {
  # create a function to calculate the theoretical correlation
  # between patients' toxicity and efficacy outcomes
  # mu: mean of BiNormal distr
  # Sigma: covariance matrix
  # Th: threshold to make Bernoulli r.m.
  bi_den <- function(x1, x2) {
    sapply(x2, function(t) mvtnorm::dmvnorm(x = c(x1, t), mean = mu, sigma = Sigma))
  }
  marg_expectation <- function(x1) {
    res <- sapply(x1, FUN = function(lt) {
      integrate(f = function(y) y * bi_den(lt, y), lower = -Inf, upper = Inf)$value
    })
    return(res)
  }
  # integration part
  r1 <- integrate(f = marg_expectation, lower = -Inf, upper = Th)
  #print(paste0("Outer abs err: ", r1$abs.error))
  # std
  std <- sqrt(diag(Sigma))
  # pr1
  pr1 <- pnorm(Th, mu[1], std[1])
  
  (r1$value - mu[2] * pr1) / (std[2] * sqrt(pr1 - pr1 ** 2))
}

# Same as out_corr, but much faster and more stable
out_corr2 <- function(mu, Sigma, Th) {
  # create a function to calculate the theoretical correlation
  # between patients' toxicity and efficacy outcomes
  # mu: mean of BiNormal distr
  # Sigma: covariance matrix
  # Th: threshold to make Bernoulli r.m.
  marg_expect <- function(mu, Sigma, Th) {
    r1 <- integrate(f = function(x) x * dnorm(x, mean = mu[1], sd = sqrt(Sigma[1, 1])),
                    lower = -Inf, upper = Th)$value
    
    std <- sqrt(diag(Sigma))
    rho <- Sigma[1, 2] / (std[1] * std[2])
    
    (rho * std[2] * r1 / std[1]) + (mu[2] - std[2] * mu[1] * rho / std[1]) * pnorm(Th, mu[1], std[1])
  }
  # marginal expectation
  r1 <- marg_expect(mu, Sigma, Th)
  # std
  std <- sqrt(diag(Sigma))
  # pr1
  pr1 <- pnorm(Th, mu[1], std[1])
  
  (r1 - mu[2] * pr1) / (std[2] * sqrt(pr1 - pr1 ** 2))
}

# first derivative of the correlation (Pearson)
# w.r.t the threshold T
dcor1 <- function(t, mu, std, rho) {
  # t: threshold
  # mu and std: mean and std value vectors
  # rho: BiNorm correlation
  intT <- integrate(f = function(x) x * dnorm(x, mean = mu[1], sd = std[1]),
                    lower = -Inf, upper = t)$value
  pr1 <- pnorm(t, mu[1], std[1]) # Phi(t; mu1, sd1)
  deno <- std[2] ** 2 * (pr1 - pr1 ** 2) ** 1.5 # denominator
  
  r1 <- rho * std[2] ** 2 * (t - mu[1]) * dnorm(t, mu[1], std[1]) * (pr1 - pr1 ** 2) / std[1]
  r2 <- rho * std[2] ** 2 * dnorm(t, mu[1], std[1]) * (1 - 2 * pr1) * (intT - mu[1] * pr1) / (2 * std[1])
  
  (r1 - r2) / deno
}