# Example functions for data generation #

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


par_for_trt <- function(muc, sek, stk, sigmaT, Th) {
  # muc = (mu_{1ck}, mu_{2ck})
  # sek, stk are tox and eff treatment effects for sub-trial k
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

gen_data <- function(ntrial, trial_size, muct, muce, st_all, se_all, Th, fcov) {
  # ntrial: number of sub-trials
  # trial_size: sample sizes in each ctrl and trt group (equal allocation rate)
  # muct, muce: mean of ctrl tox (eff) data (w.r.t each sub-trial)
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
  
  r1 <- rho * std[2] ** 2 * (t - mu[1]) * dnorm(t, mu[1], std[1]) * (pr1 - pr1 ** 2) / std[1]
  r2 <- rho * std[2] ** 2 * dnorm(t, mu[1], std[1]) * (1 - 2 * pr1) * (intT - mu[1] * pr1) / (2 * std[1])
  
  (r1 - r2) / deno
}
