# This file is to store simulation and analytical functions for different models

crd_BiEXNEX <- function(ntrial, trial_size, dat_temp, w) {
  # create two plots to see the estimations of CRDs in a single simulation
  # only suitable for the combined-weight model
  # dat_temp: observations for model fitting
  # w: prior exchangeability weights
  
  # all data needed for fitting
  dat_joint <- list(
    y = dat_temp,
    n_mod = ntrial,
    n = trial_size,
    n_max = max(trial_size),
    m0 = c(0, 0),
    s0 = c(5, 5),
    m = c(0, 0),
    s = c(5, 5),
    z = c(0.5, 0.5),
    paT = c(0, 10),
    paE = c(0, 10),
    std_var_eff = c(1, 1),
    pw = c(w, 1 - w)
  )
  
  # fit the model once
  fit_joint <- jags.parallel(data = dat_joint,
                             parameters.to.save = c("theta"),
                             model.file = "../model/EXNEX_COMB.txt",
                             n.chains = 4,
                             n.iter = 1e4,
                             n.burnin = 2e3,
                             n.thin = 2,
                             n.cluster = 4)
  
  # plot 90% credible intervals
  mcmc_fit <- as.mcmc(fit_joint)
  p1 <- mcmc_intervals(mcmc_fit, pars = paste0("theta[", 1, ",", 1:ntrial, "]"),
                       prob = 0.9, prob_outer = 0.95) + ggtitle(
                         "Treatment effects for toxicity outcomes",
                         "with 90% posterior credible intervals"
                       )
  
  p2 <- mcmc_intervals(mcmc_fit, pars = paste0("theta[", 2, ",", 1:ntrial, "]"),
                       prob = 0.9, prob_outer = 0.95) + ggtitle(
                         "Treatment effects for efficacy outcomes",
                         "with 90% posterior credible intervals"
                       )
  grid.arrange(p1, p2, ncol = 2)
}

crd_Stdalone <- function(ntrial, trial_size, dat_temp) {
  # create two plots to see the estimations of CRDs in a single simulation
  # only suitable for the stand-alone model
  
  # parameter definitions are the same as those of gen_data
  
  # all data needed for fitting
  dat_joint <- list(
    y = dat_temp,
    n_mod = ntrial,
    n = trial_size,
    n_max = max(trial_size),
    paT = c(0, 10),
    paE = c(0, 10),
    std_var_eff = c(1, 1),
    pte_tox = matrix(c(rep(0, n_mod), rep(1, n_mod)), ncol = 2),
    pte_eff = matrix(c(rep(0, n_mod), rep(1, n_mod)), ncol = 2)
  )
  
  # fit the model once
  fit_joint <- jags.parallel(data = dat_joint,
                             parameters.to.save = c("theta"),
                             model.file = "../model/Stand-alone.txt",
                             n.chains = 4,
                             n.iter = 1e4,
                             n.burnin = 2e3,
                             n.thin = 2,
                             n.cluster = 4)
  
  # plot 90% credible intervals
  mcmc_fit <- as.mcmc(fit_joint)
  p1 <- mcmc_intervals(mcmc_fit, pars = paste0("theta[", 1, ",", 1:ntrial, "]"),
                       prob = 0.9, prob_outer = 0.95) + ggtitle(
                         "Treatment effects for toxicity outcomes",
                         "with 90% posterior credible intervals"
                       )
  
  p2 <- mcmc_intervals(mcmc_fit, pars = paste0("theta[", 2, ",", 1:ntrial, "]"),
                       prob = 0.9, prob_outer = 0.95) + ggtitle(
                         "Treatment effects for efficacy outcomes",
                         "with 90% posterior credible intervals"
                       )
  grid.arrange(p1, p2, ncol = 2)
}

crd_IndEXNEX <- function(ntrial, trial_size, dat_temp, w1, w2) {
  # create two plots to see the estimations of CRDs in a single simulation
  # only suitable for the combined-weight model
  # dat_temp: observations for model fitting
  # w1, w1: prior exchangeability weights for tox and eff
  
  # all data needed for fitting
  dat_tox <- list(
    y = dat_temp[dim(dat_temp)[1], , ],
    n_mod = ntrial,
    n = trial_size,
    pw = c(w1, 1 - w1),
    mu_ex_prior = c(0, 5),
    std_var_ex = 0.5,
    mu_nex = 0,
    std_nex = 5,
    pa = 0,
    std_alpha = 10
  )
  
  dat_eff <- list(
    y = dat_temp[-dim(dat_temp)[1], , ],
    n_mod = ntrial,
    n = trial_size,
    pw = c(w2, 1 - w2),
    mu_ex_prior = c(0, 5),
    std_var_ex = 0.5,
    mu_nex = 0,
    std_nex = 5,
    pa = 0,
    std_alpha = 10,
    std_var_eff = c(1, 1)
  )
  # fit the model once
  fit_tox <- jags.parallel(data = dat_tox,
                             parameters.to.save = c("theta"),
                             model.file = "../model/EXNEX_binom.txt",
                             n.chains = 4,
                             n.iter = 1e4,
                             n.burnin = 2e3,
                             n.thin = 2,
                             n.cluster = 4)
  fit_eff <- jags.parallel(data = dat_eff,
                           parameters.to.save = c("theta"),
                           model.file = "../model/EXNEX_norm.txt",
                           n.chains = 4,
                           n.iter = 1e4,
                           n.burnin = 2e3,
                           n.thin = 2,
                           n.cluster = 4)
  
  # plot 90% credible intervals
  mcmc_tox <- as.mcmc(fit_tox)
  mcmc_eff <- as.mcmc(fit_eff)
  
  p1 <- mcmc_intervals(mcmc_tox, pars = paste0("theta[", 1:ntrial, "]"),
                       prob = 0.9, prob_outer = 0.95) + ggtitle(
                         "Treatment effects for toxicity outcomes",
                         "with 90% posterior credible intervals"
                       )
  
  p2 <- mcmc_intervals(mcmc_eff, pars = paste0("theta[", 1:ntrial, "]"),
                       prob = 0.9, prob_outer = 0.95) + ggtitle(
                         "Treatment effects for efficacy outcomes",
                         "with 90% posterior credible intervals"
                       )
  grid.arrange(p1, p2, ncol = 2)
}

par_summary <- function(fit, model_name, params = "theta", alpha = 0.05) {
  # This function is used for summarising posterior samples
  # from given parameters and the rjags object
  # model_name: a string of the name of `fit`
  # alpha: posterior (alpha, 1 - alpha) interval
  
  # make data.frame
  mcmc_obj <- coda::as.mcmc(fit)
  df_obj <- as.data.frame(MCMCvis::MCMCchains(mcmc_obj, params = params))
  # summary
  df_sum <- t(do.call(cbind.data.frame,
                      args = lapply(as.list(df_obj),
                                    FUN = function(x) quantile(x, c(alpha, 0.5, 1 - alpha)))))
                                    
  colnames(df_sum) <- c("low", "med", "up")
  df_sum <- as.data.frame(df_sum)
  
  # add identity variable
  df_sum[, "par"] <- rownames(df_sum)
  df_sum[, "Model"] <- rep(model_name, nrow(df_sum))
  
  return(df_sum)
}

theta_compare <- function(fit_list, fit_name, tr_val = NULL, em_val) {
  # fit_list is the list of fitted models with "theta" stored
  # fit_name is a vector including corresponding model names
  # tr_val is a list(tox = ..., eff = ...) of true par values
  # em_val are the observed values of the true parameter
  
  #processing and combining
  lsum <- mapply(FUN = function(x, y) par_summary(x, y),
                 fit_list, fit_name,
                 SIMPLIFY = FALSE)
  # combine them into a summary data.frame
  df_sum <- do.call(rbind.data.frame, lsum)
  
  # add true parameter values
  df_trval <- data.frame(par = unique(df_sum$par),
                         #trval = c(tr_val$tox, tr_val$eff),
                         emval = c(em_val$tox, em_val$eff))
  df_sum <- merge(df_sum, df_trval, by = "par")
  
  # split it according to Tox and Eff
  df_sum_tox <- df_sum[which(df_sum$par %in% paste0("theta[1,", 1:6, "]")), ]
  df_sum_eff <- df_sum[which(df_sum$par %in% paste0("theta[2,", 1:6, "]")), ]
  
  # plotting
  p1 <- ggplot(data = df_sum_tox, mapping = aes(x = par, color = Model)) +
    geom_errorbar(width = .5, mapping = aes(ymin = low, ymax = up), position = position_dodge(0.5)) +
    geom_point(mapping = aes(y = med), position = position_dodge(0.5)) +
    #geom_point(mapping = aes(y = trval), shape = 4, size = 2.5, color = "black", position = position_dodge(0.5)) +
    geom_point(mapping = aes(y = emval), shape = 4, size = 2.5, color = "red", position = position_dodge(0.5)) +
    coord_flip() + labs(x = "", y = "") + guides(color = "none") +
    scale_x_discrete(labels = expression(theta[1]^t, theta[2]^t, theta[3]^t,
                                         theta[4]^t, theta[5]^t, theta[6]^t)) +
    ggtitle(paste0("Comparison of treatment differences between ", length(fit_list), " models"),
            "with 90% posterior intervals")
  
  p2 <- ggplot(data = df_sum_eff, mapping = aes(x = par, color = Model)) +
    geom_errorbar(width = .5, mapping = aes(ymin = low, ymax = up), position = position_dodge(0.5)) +
    geom_point(mapping = aes(y = med), position = position_dodge(0.5)) +
    #geom_point(mapping = aes(y = trval), shape = 4, size = 2.5, color = "black", position = position_dodge(0.5)) +
    geom_point(mapping = aes(y = emval), shape = 4, size = 2.5, color = "red", position = position_dodge(0.5)) +
    scale_x_discrete(labels = expression(theta[1]^e, theta[2]^e, theta[3]^e,
                                         theta[4]^e, theta[5]^e, theta[6]^e)) +
    coord_flip() + labs(x = "", y = "") +
    theme(legend.position = "top")
  grid.arrange(p1, p2, ncol = 2)
}

# explore the correlation between crd in the combined weight model
corr_BiEXNEX <- function(ntrial, trial_size, dat_temp, w) {
  # only suitable for the combined-weight model
  # dat_temp: observations for model fitting
  # w: prior exchangeability weights
  
  # all data needed for fitting
  dat_joint <- list(
    y = dat_temp,
    n_mod = ntrial,
    n = trial_size,
    n_max = max(trial_size),
    m0 = c(0, 0),
    s0 = c(5, 5),
    m = c(0, 0),
    s = c(5, 5),
    z = c(0.5, 0.5),
    paT = c(0, 10),
    paE = c(0, 10),
    std_var_eff = c(1, 1),
    pw = c(w, 1 - w)
  )
  
  # fit the model once
  fit_joint <- jags.parallel(data = dat_joint,
                             parameters.to.save = c("rho1", "rho2"),
                             model.file = "../model/EXNEX_COMB.txt",
                             n.chains = 4,
                             n.iter = 1e4,
                             n.burnin = 2e3,
                             n.thin = 2,
                             n.cluster = 4)
  
  # plot 90% credible intervals
  mcmc_fit <- as.mcmc(fit_joint)
  p1 <- mcmc_areas(mcmc_fit, pars = "rho1", prob = 0.9) + ggtitle(
                         "Correlation for EX",
                         "with 90% posterior credible intervals"
                       )
  
  p2 <- mcmc_areas(mcmc_fit, pars = "rho2", prob = 0.9) + ggtitle(
                         "Correlation for NEX",
                         "with 90% posterior credible intervals"
                       )
  grid.arrange(p1, p2, ncol = 2)
}

corr_PE <- function(ntrial, trial_size, dat_temp, wt = 0.5, we = 0.5) {
  # dat_temp: observations for model fitting
  # wt, we: prior exchangeability weights
  
  # all data needed for fitting
  dat_joint <- list( y = dat_temp,
                  n_mod = ntrial,
                  n = trial_size,
                  n_max = max(trial_size),
                  m0 = c(0, 0),
                  s0 = c(5, 5),
                  m = c(0, 0),
                  s = c(5, 5),
                  z = c(1, 1),
                  paT = c(0, 10),
                  paE = c(0, 10),
                  std_var_eff = c(5, 5),
                  pwT = c(wt, 1 - wt),
                  pwE = c(we, 1 - we)
  )
  fit_joint <- jags.parallel(data = dat_joint,
                          parameters.to.save = c("rho1", "rho2"),
                          model.file = "../model/Panel-EXNEX.txt",
                          n.chains = 4,
                          n.iter = 1e4,
                          n.burnin = 2e3,
                          n.thin = 2,
                          n.cluster = 4)
  
  # plot 90% credible intervals
  mcmc_fit <- as.mcmc(fit_joint)
  p1 <- mcmc_areas(mcmc_fit, pars = "rho1", prob = 0.9) + ggtitle(
    "Correlation for EX",
    "with 90% posterior credible intervals"
  )
  
  p2 <- mcmc_areas(mcmc_fit, pars = "rho2", prob = 0.9) + ggtitle(
    "Correlation for NEX",
    "with 90% posterior credible intervals"
  )
  grid.arrange(p1, p2, ncol = 2)
}

# explore posterior exchangeability weights
postw_BiEXNEX <- function(ntrial, trial_size, dat_temp, w) {
  # all data needed for fitting
  dat_joint <- list(
    y = dat_temp,
    n_mod = ntrial,
    n = trial_size,
    n_max = max(trial_size),
    m0 = c(0, 0),
    s0 = c(5, 5),
    m = c(0, 0),
    s = c(5, 5),
    z = c(0.5, 0.5),
    paT = c(0, 10),
    paE = c(0, 10),
    std_var_eff = c(5, 5),
    pw = c(w, 1 - w)
  )
  
  # fit the model once
  fit_joint <- jags.parallel(data = dat_joint,
                             parameters.to.save = c("cat"),
                             model.file = "../model/EXNEX_COMB.txt",
                             n.chains = 4,
                             n.iter = 8e3,
                             n.burnin = 2e3,
                             n.thin = 2,
                             n.cluster = 4)
  
  # mcmc object
  mcmc_fit <- as.mcmc(fit_joint)
  extract_w <- codatools::coda_df(mcmc_fit,
               parameters = grep("cat", coda::varnames(mcmc_fit), value = T))[, -c(1, 2)]
  
  fr <- colMeans(2 - extract_w)
  names(fr) <- paste0("Sub-trial", 1:6)
  return(fr)
}

postw_IndEXNEX <- function(ntrial, trial_size, dat_temp, w1, w2) {
  # all data needed for fitting
  dat_tox <- list(
    y = dat_temp[dim(dat_temp)[1], , ],
    n_mod = ntrial,
    n = trial_size,
    pw = c(w1, 1 - w1),
    mu_ex_prior = c(0, 5),
    std_var_ex = 1,
    mu_nex = 0,
    std_nex = 5,
    pa = 0,
    std_alpha = 10
  )
  
  dat_eff <- list(
    y = dat_temp[-dim(dat_temp)[1], , ],
    n_mod = ntrial,
    n = trial_size,
    pw = c(w2, 1 - w2),
    mu_ex_prior = c(0, 5),
    std_var_ex = 1,
    mu_nex = 0,
    std_nex = 5,
    pa = 0,
    std_alpha = 10,
    std_var_eff = c(5, 5)
  )
  # fit the model once
  fit_tox <- jags.parallel(data = dat_tox,
                           parameters.to.save = c("catT"),
                           model.file = "../model/EXNEX_binom.txt",
                           n.chains = 4,
                           n.iter = 8e3,
                           n.burnin = 2e3,
                           n.thin = 2,
                           n.cluster = 4)
  fit_eff <- jags.parallel(data = dat_eff,
                           parameters.to.save = c("catE"),
                           model.file = "../model/EXNEX_norm.txt",
                           n.chains = 4,
                           n.iter = 8e3,
                           n.burnin = 2e3,
                           n.thin = 2,
                           n.cluster = 4)
  
  # mcmc calculation
  mcmc_tox <- as.mcmc(fit_tox)
  mcmc_eff <- as.mcmc(fit_eff)
  extract_w1 <- codatools::coda_df(mcmc_tox,
                                  parameters = grep("catT", coda::varnames(mcmc_tox), value = T))[, -c(1, 2)]
  extract_w2 <- codatools::coda_df(mcmc_eff,
                                   parameters = grep("catE", coda::varnames(mcmc_eff), value = T))[, -c(1, 2)]
  
  res <- as.data.frame(
    rbind(colMeans(2 - extract_w1), colMeans(2 - extract_w2))
  )
  colnames(res) <- paste0("Sub-trial", 1:6)
  rownames(res) <- c("Tox", "Eff")
  
  return(res)
}

postw_PE <- function(dat, n_mod, n, wt = 0.5, we = 0.5) {
  # posterior exchangeability weights for PE model
  # default prior w = 0.5
  # for Panel EXNEX
  dat_pi <- list(
    y = dat,
    n_mod = n_mod,
    n = n,
    n_max = max(n),
    m0 = c(0, 0),
    s0 = c(5, 5),
    m = c(0, 0),
    s = c(5, 5),
    z = c(1, 1),
    paT = c(0, 10),
    paE = c(0, 10),
    std_var_eff = c(5, 5),
    pwT = c(wt, 1 - wt),
    pwE = c(we, 1 - we)
  )
  
  fit_pi <- jags.parallel(data = dat_pi,
                          parameters.to.save = c("catT", "catE"),
                          model.file = "../model/Panel-EXNEX.txt",
                          n.chains = 4,
                          n.iter = 2e4,
                          n.burnin = 2e3,
                          n.thin = 2,
                          n.cluster = 4)
  
  mcmc_pi <- as.mcmc(fit_pi)
  # weight
  df_pi <- coda_df(mcmc_pi)
  
  # result list
  res <- colMeans(2 - df_pi[, 3:(2 + 2 * n_mod)])
  
  # trun into a data.frame
  res <- rbind(res[-(1:n_mod)], res[1:n_mod])
  colnames(res) <- paste0("Sub-trial", 1:6)
  rownames(res) <- c("Tox", "Eff")
  
  return(res)
}

postw_Panel_4C <- function(dat, n_mod, n) {
  # given the observation array input
  # compare the posterior weight for each sub-trial
  # between Panel and 4C EXNEX mixture models
  dat_4c <- list(
    y = dat,
    n_mod = n_mod,
    n = n,
    n_max = max(n),
    m0 = c(0, 0),
    s0 = c(5, 5),
    m = c(0, 0),
    s = c(5, 5),
    z = c(1, 1),
    paT = c(0, 10),
    paE = c(0, 10),
    std_var_eff = c(5, 5),
    pw = rep(0.25, 4)
  )
  
  fit_4c <- jags.parallel(data = dat_4c,
                          parameters.to.save = c("cat"),
                          model.file = "../model/4C-EXNEX.txt",
                          n.chains = 4,
                          n.iter = 2e4,
                          n.burnin = 2e3,
                          n.thin = 2,
                          n.cluster = 4)
  mcmc_4c <- as.mcmc(fit_4c)
  # weight
  df_4c <- coda_df(mcmc_4c)
  # for Panel EXNEX
  dat_pi <- list(
    y = dat,
    n_mod = n_mod,
    n = n,
    n_max = max(n),
    m0 = c(0, 0),
    s0 = c(5, 5),
    m = c(0, 0),
    s = c(5, 5),
    z = c(1, 1),
    paT = c(0, 10),
    paE = c(0, 10),
    std_var_eff = c(5, 5),
    pwT = c(0.5, 1 - 0.5),
    pwE = c(0.5, 1 - 0.5)
  )
  
  fit_pi <- jags.parallel(data = dat_pi,
                          parameters.to.save = c("catT", "catE"),
                          model.file = "../model/Panel-EXNEX.txt",
                          n.chains = 4,
                          n.iter = 2e4,
                          n.burnin = 2e3,
                          n.thin = 2,
                          n.cluster = 4)
  
  mcmc_pi <- as.mcmc(fit_pi)
  # weight
  df_pi <- coda_df(mcmc_pi)
  
  # result list
  res <- list()
  res$grid <- colMeans(2 - df_pi[, 3:(2 + 2 * n_mod)])
  # store 4c results
  res$fc <- matrix(, ncol = 4, nrow = n_mod)
  for (k in 1:n_mod) {
    for (j in 1:4) {
      res$fc[k, j] <- sum(df_4c[, 2 + k] == j) / nrow(df_4c)
    }
  }
  # grid to 4c results - post. w.
  res$grid2fc <- matrix(, ncol = 4, nrow = n_mod)
  res$grid2fc[, 1] <- res$grid[1:n_mod] * res$grid[-(1:n_mod)]
  res$grid2fc[, 2] <- (1 - res$grid[1:n_mod]) * res$grid[-(1:n_mod)]
  res$grid2fc[, 3] <- res$grid[1:n_mod] * (1 - res$grid[-(1:n_mod)])
  res$grid2fc[, 4] <- (1 - res$grid[1:n_mod]) * (1 - res$grid[-(1:n_mod)])
  
  return(res)
}

# explore frequentist properties of each model
sim_BiEXNEX <- function(ntrial, trial_size, muct, muce, st_all, se_all, Th, fcov,
                num_cl = 2, n_sim = 500, delta = 0, w = 0.5, eta = c(0.95, 0.95)) {
  # num_cl: number of threads you want to use
  # n_sim: number of simulations
  # mp: the location where the model exists
  # other parameters are the same as those in "gen_data"
  
  # set the number of clusters used in the process; more != better
  cl <- makeCluster(num_cl)
  clusterExport(cl, varlist = c("gen_data", "gen_binorm_data", "par_for_trt",
                                "logit", "logistic"))
  registerDoParallel(cl)
  
  # main body
  sim_comb <- foreach(i = 1:n_sim, .combine = rbind,
                      .packages = c("R2jags", "qpcR", "mvtnorm")) %dopar% {
                        
    # generate the outcomes
    sim_obs <- gen_data(ntrial, trial_size, muct, muce, st_all, se_all, Th, fcov)
    # combine all relevant data for fitting the combined model
    sim_dat <- list(
      y = sim_obs,
      n_mod = ntrial,
      n = trial_size,
      n_max = max(trial_size),
      m0 = c(0, 0),
      s0 = c(5, 5),
      m = c(0, 0),
      s = c(5, 5),
      z = c(0.5, 0.5),
      paT = c(0, 10),
      paE = c(0, 10),
      std_var_eff = c(5, 5),
      pw = c(w, 1 - w)
    )
                        
    # fitting
    sim_fit_joint <- jags(data = sim_dat,
                          parameters.to.save = c("theta"),
                          model.file = "~/PhaseII_basket_EXNEX/model/EXNEX_COMB.txt",
                          n.chains = 4,
                          n.iter = 1e4,
                          n.burnin = 2e3,
                          n.thin = 2)
    # transform the fitting object into a mcmc object
    sim_mcmc_joint <- coda::as.mcmc(sim_fit_joint)
    df_sim_comb <- codatools::coda_df(sim_mcmc_joint,
                   parameters = grep("theta", coda::varnames(sim_mcmc_joint), value = T))[, -c(1, 2)]
    
    tsim_comb <- c()                    
    tsim_comb[1:ntrial] <- (colMeans(df_sim_comb[, 1:ntrial] < 0) > eta[1])
    tsim_comb[(ntrial + 1):(2 * ntrial)] <- (colMeans(df_sim_comb[, (ntrial + 1):(2 * ntrial)] > delta) > eta[2])
    
    tsim_comb <- tsim_comb[1:ntrial] * tsim_comb[(ntrial + 1):(2 * ntrial)]
    tsim_comb
  }
  stopCluster(cl)
  return(sim_comb)
}

sim_IndEXNEX <- function(ntrial, trial_size, muct, muce, st_all, se_all, Th, fcov,
                        num_cl = 2, n_sim = 500, delta = 0, eta = c(0.95, 0.95), w1 = 0.5, w2 = 0.5) {
  # num_cl: number of threads you want to use
  # n_sim: number of simulations
  # mp: the location where the model exists
  # other parameters are the same as those in "gen_data"
  
  # set the number of clusters used in the process; more != better
  cl <- makeCluster(num_cl)
  clusterExport(cl, varlist = c("gen_data", "gen_binorm_data", "par_for_trt",
                                "logit", "logistic"))
  registerDoParallel(cl)
  
  # main body
  sim_comb <- foreach(i = 1:n_sim, .combine = rbind,
                      .packages = c("R2jags", "qpcR", "mvtnorm")) %dopar% {
                        
                        # generate the outcomes
                        sim_obs <- gen_data(ntrial, trial_size, muct, muce, st_all, se_all, Th, fcov)
                        # combine all relevant data for fitting the combined model
                        sim_tox <- list(
                          y = sim_obs[dim(sim_obs)[1], , ],
                          n_mod = ntrial,
                          n = trial_size,
                          pw = c(w1, 1 - w1),
                          mu_ex_prior = c(0, 5),
                          std_var_ex = 0.5,
                          mu_nex = 0,
                          std_nex = 5,
                          pa = 0,
                          std_alpha = 10
                        )
                        
                        sim_eff <- list(
                          y = sim_obs[-dim(sim_obs)[1], , ],
                          n_mod = ntrial,
                          n = trial_size,
                          pw = c(w2, 1 - w2),
                          mu_ex_prior = c(0, 5),
                          std_var_ex = 0.5,
                          mu_nex = 0,
                          std_nex = 5,
                          pa = 0,
                          std_alpha = 10,
                          std_var_eff = c(5, 5)
                        )
                        
                        # fitting
                        sim_fit_tox <- jags(data = sim_tox,
                                                 parameters.to.save = c("theta"),
                                                 model.file = "~/PhaseII_basket_EXNEX/model/EXNEX_binom.txt",
                                                 n.chains = 4,
                                                 n.iter = 1e4,
                                                 n.burnin = 2e3,
                                                 n.thin = 2)
                        sim_fit_eff <- jags(data = sim_eff,
                                                 parameters.to.save = c("theta"),
                                                 model.file = "~/PhaseII_basket_EXNEX/model/EXNEX_norm.txt",
                                                 n.chains = 4,
                                                 n.iter = 1e4,
                                                 n.burnin = 2e3,
                                                 n.thin = 2)
                        # transform the fitting object into a mcmc object
                        sim_mcmc_tox <- coda::as.mcmc(sim_fit_tox)
                        sim_mcmc_eff <- coda::as.mcmc(sim_fit_eff)
                        
                        # extract fitting results
                        df_sim_tox <- codatools::coda_df(sim_mcmc_tox,
                                       parameters = grep("theta", coda::varnames(sim_mcmc_tox), value = T))[, -c(1, 2)]
                        df_sim_eff <- codatools::coda_df(sim_mcmc_eff,
                                                         parameters = grep("theta", coda::varnames(sim_mcmc_eff), value = T))[, -c(1, 2)]
                        
                        df_sim_tox <- (colMeans(df_sim_tox < 0) > eta[1])
                        df_sim_eff <- (colMeans(df_sim_eff > delta) > eta[2])
                        
                        df_results <- df_sim_tox * df_sim_eff
                        df_results
                        
                      }
  stopCluster(cl)
  return(sim_comb)
}

sim_Stdalone <- function(ntrial, trial_size, muct, muce, st_all, se_all, Th, fcov,
                        num_cl = 2, n_sim = 500, delta = 0, eta = c(0.95, 0.95)) {
  # num_cl: number of threads you want to use
  # n_sim: number of simulations
  # mp: the location where the model exists
  # other parameters are the same as those in "gen_data"
  
  # set the number of clusters used in the process; more != better
  cl <- makeCluster(num_cl)
  clusterExport(cl, varlist = c("gen_data", "gen_binorm_data", "par_for_trt",
                                "logit", "logistic"))
  registerDoParallel(cl)
  
  # main body
  sim_comb <- foreach(i = 1:n_sim, .combine = rbind,
                      .packages = c("R2jags", "qpcR", "mvtnorm")) %dopar% {
                        
                        # generate the outcomes
                        sim_obs <- gen_data(ntrial, trial_size, muct, muce, st_all, se_all, Th, fcov)
                        # combine all relevant data for fitting the combined model
                        sim_dat <- list(
                          y = sim_obs,
                          n_mod = ntrial,
                          n = trial_size,
                          n_max = max(trial_size),
                          paT = c(0, 10),
                          paE = c(0, 10),
                          std_var_eff = c(5, 5),
                          pte_tox = matrix(c(rep(0, ntrial), rep(5, ntrial)), ncol = 2),
                          pte_eff = matrix(c(rep(0, ntrial), rep(5, ntrial)), ncol = 2)
                        )
                        
                        # fitting
                        sim_fit_joint <- jags(data = sim_dat,
                                              parameters.to.save = c("theta"),
                                              model.file = "~/PhaseII_basket_EXNEX/model/Stand-alone.txt",
                                              n.chains = 4,
                                              n.iter = 1e4,
                                              n.burnin = 2e3,
                                              n.thin = 2)
                        # transform the fitting object into a mcmc object
                        sim_mcmc_joint <- coda::as.mcmc(sim_fit_joint)
                        df_sim_comb <- codatools::coda_df(sim_mcmc_joint,
                                       parameters = grep("theta", coda::varnames(sim_mcmc_joint), value = T))[, -c(1, 2)]
                        
                        tsim_comb <- c()                    
                        tsim_comb[1:ntrial] <- (colMeans(df_sim_comb[, 1:ntrial] < 0) > eta[1])
                        tsim_comb[(ntrial + 1):(2 * ntrial)] <- (colMeans(df_sim_comb[, (ntrial + 1):(2 * ntrial)] > delta) > eta[2])
                        
                        tsim_comb <- tsim_comb[1:ntrial] * tsim_comb[(ntrial + 1):(2 * ntrial)]
                        tsim_comb
                        
                      }
  stopCluster(cl)
  return(sim_comb)
}

sim_PanelEXNEX <- function(ntrial, trial_size, muct, muce, st_all, se_all, Th, fcov,
                        num_cl = 2, n_sim = 500, delta = 0, wt = 0.5, we = 0.5, eta1 = 0.95, eta2 = 0.95) {
  # num_cl: number of threads you want to use
  # n_sim: number of simulations
  # mp: the location where the model exists
  # other parameters are the same as those in "gen_data"
  
  # set the number of clusters used in the process; more != better
  cl <- makeCluster(num_cl)
  clusterExport(cl, varlist = c("gen_data", "gen_binorm_data", "par_for_trt",
                                "logit", "logistic"))
  registerDoParallel(cl)
  
  # main body
  sim_comb <- foreach(i = 1:n_sim, .combine = rbind,
                      .packages = c("R2jags", "qpcR", "mvtnorm")) %dopar% {
                        
                        # generate the outcomes
                        sim_obs <- gen_data(ntrial, trial_size, muct, muce, st_all, se_all, Th, fcov)
                        # combine all relevant data for fitting the combined model
                        sim_dat <- list(
                          y = sim_obs,
                          n_mod = ntrial,
                          n = trial_size,
                          n_max = max(trial_size),
                          m0 = c(0, 0),
                          s0 = c(5, 5),
                          m = c(0, 0),
                          s = c(5, 5),
                          z = c(0.5, 0.5),
                          paT = c(0, 10),
                          paE = c(0, 10),
                          std_var_eff = c(5, 5),
                          pwT = c(wt, 1 - wt),
                          pwE = c(we, 1 - we)
                        )
                        
                        # fitting
                        sim_fit_joint <- jags(data = sim_dat,
                                              parameters.to.save = c("theta"),
                                              model.file = "~/PhaseII_basket_EXNEX/model/Panel-EXNEX.txt",
                                              n.chains = 4,
                                              n.iter = 1e4,
                                              n.burnin = 2e3,
                                              n.thin = 2)
                        # transform the fitting object into a mcmc object
                        sim_mcmc_joint <- coda::as.mcmc(sim_fit_joint)
                        df_sim_comb <- codatools::coda_df(sim_mcmc_joint,
                                                          parameters = grep("theta", coda::varnames(sim_mcmc_joint), value = T))[, -c(1, 2)]
                        
                        tsim_comb <- c()                    
                        tsim_comb[1:ntrial] <- (colMeans(df_sim_comb[, 1:ntrial] < 0) > eta1)
                        tsim_comb[(ntrial + 1):(2 * ntrial)] <- (colMeans(df_sim_comb[, (ntrial + 1):(2 * ntrial)] > delta) > eta2)
                        
                        tsim_comb <- tsim_comb[1:ntrial] * tsim_comb[(ntrial + 1):(2 * ntrial)]
                        tsim_comb
                      }
  stopCluster(cl)
  return(sim_comb)
}

sim_EBiEXNEX <- function(ntrial, trial_size, muct, muce, st_all, se_all, Th, fcov,
                         w1 = 0.25, w2 = 0.25, w3 = 0.25,
                         num_cl = 2, n_sim = 500, delta = 0, eta = c(0.95, 0.95),
                         model_file = "~/PhaseII_basket_EXNEX/model/E-BiEXNEX.txt") {
  # num_cl: number of threads you want to use
  # n_sim: number of simulations
  # mp: the location where the model exists
  # other parameters are the same as those in "gen_data"
  
  # set the number of clusters used in the process; more != better
  cl <- makeCluster(num_cl)
  clusterExport(cl, varlist = c("gen_data", "gen_binorm_data", "par_for_trt",
                                "logit", "logistic"))
  registerDoParallel(cl)
  
  # main body
  sim_comb <- foreach(i = 1:n_sim, .combine = rbind,
                      .packages = c("R2jags", "qpcR", "mvtnorm")) %dopar% {
                        
                        # generate the outcomes
                        sim_obs <- gen_data(ntrial, trial_size, muct, muce, st_all, se_all, Th, fcov)
                        # combine all relevant data for fitting the combined model
                        sim_dat <- list(
                          y = sim_obs,
                          n_mod = ntrial,
                          n = trial_size,
                          n_max = max(trial_size),
                          m0 = c(0, 0),
                          s0 = c(5, 5),
                          m = c(0, 0),
                          s = c(5, 5),
                          z = c(0.5, 0.5),
                          paT = c(0, 10),
                          paE = c(0, 10),
                          std_var_eff = c(5, 5),
                          pw = c(w1, w2, w3, 1 - w1 - w2 - w3)
                        )
                        
                        # fitting
                        sim_fit_joint <- jags(data = sim_dat,
                                              parameters.to.save = c("theta"),
                                              model.file = model_file,
                                              n.chains = 4,
                                              n.iter = 1e4,
                                              n.burnin = 2e3,
                                              n.thin = 2)
                        # transform the fitting object into a mcmc object
                        sim_mcmc_joint <- coda::as.mcmc(sim_fit_joint)
                        df_sim_comb <- codatools::coda_df(sim_mcmc_joint,
                                                          parameters = grep("theta", coda::varnames(sim_mcmc_joint), value = T))[, -c(1, 2)]
                        
                        tsim_comb <- c()                    
                        tsim_comb[1:ntrial] <- (colMeans(df_sim_comb[, 1:ntrial] < 0) > eta[1])
                        tsim_comb[(ntrial + 1):(2 * ntrial)] <- (colMeans(df_sim_comb[, (ntrial + 1):(2 * ntrial)] > delta) > eta[2])
                        
                        tsim_comb <- tsim_comb[1:ntrial] * tsim_comb[(ntrial + 1):(2 * ntrial)]
                        tsim_comb
                      }
  stopCluster(cl)
  return(sim_comb)
}

sim_BMA2 <- function(ntrial, trial_size, muct, muce, st_all, se_all, Th, fcov,
                    num_cl = 2, n_sim = 500, delta = 0, wt = 0.5, we = 0.5, eta = c(0.95, 0.95),
                    model_file = "~/PhaseII_basket_EXNEX/model/bma_4c.txt") {
  # num_cl: number of threads you want to use
  # n_sim: number of simulations
  # mp: the location where the model exists
  # other parameters are the same as those in "gen_data"
  
  # set the number of clusters used in the process; more != better
  cl <- makeCluster(num_cl)
  clusterExport(cl, varlist = c("gen_data", "gen_binorm_data", "par_for_trt",
                                "logit", "logistic"))
  registerDoParallel(cl)
  
  # main body
  sim_comb <- foreach(i = 1:n_sim, .combine = rbind,
                      .packages = c("R2jags", "qpcR", "mvtnorm"),
                      .errorhandling = "remove") %dopar% {
                        
                        # generate the outcomes
                        sim_obs <- gen_data(ntrial, trial_size, muct, muce, st_all, se_all, Th, fcov)
                        # combine all relevant data for fitting the combined model
                        sim_dat <- list(
                          # observations
                          y = array(rep(sim_obs, 5), dim = c(dim(sim_obs), 5)),
                          n = trial_size,
                          n_mod = ntrial,
                          n_max = max(trial_size),
                          # model 1
                          paT1 = c(0, 10),
                          paE1 = c(0, 10),
                          std_var_eff1 = c(5, 5),
                          m1 = c(0, 0),
                          s1 = c(5, 5),
                          z1 = c(0.5, 0.5),
                          # model 2
                          paT2 = c(0, 10),
                          paE2 = c(0, 10),
                          std_var_eff2 = c(5, 5),
                          m2 = c(0, 0),
                          s2 = c(5, 5),
                          z2 = 0.5,
                          # model 3
                          paT3 = c(0, 10),
                          paE3 = c(0, 10),
                          std_var_eff3 = c(5, 5),
                          m3 = c(0, 0),
                          s3 = c(5, 5),
                          z3 = 0.5,
                          # model 4
                          paT4 = c(0, 10),
                          paE4 = c(0, 10),
                          std_var_eff4 = c(5, 5),
                          m4 = c(0, 0),
                          s4 = c(5, 5),
                          # model selection
                          pwm = rep(1, 4)
                        )
                        
                        # fitting
                        sim_fit_joint <- jags(data = sim_dat,
                                              parameters.to.save = c("theta"),
                                              model.file = model_file,
                                              n.chains = 4,
                                              n.iter = 1e4,
                                              n.burnin = 2e3,
                                              n.thin = 2)
                        # transform the fitting object into a mcmc object
                        sim_mcmc_joint <- coda::as.mcmc(sim_fit_joint)
                        df_sim_comb <- codatools::coda_df(sim_mcmc_joint,
                                                          parameters = grep("theta", coda::varnames(sim_mcmc_joint), value = T))[, -c(1, 2)]
                        
                        tsim_comb <- c()                    
                        tsim_comb[1:ntrial] <- (colMeans(df_sim_comb[, 1:ntrial] < 0) > eta[1])
                        tsim_comb[(ntrial + 1):(2 * ntrial)] <- (colMeans(df_sim_comb[, (ntrial + 1):(2 * ntrial)] > delta) > eta[2])
                        
                        tsim_comb <- tsim_comb[1:ntrial] * tsim_comb[(ntrial + 1):(2 * ntrial)]
                        tsim_comb
                      }
  stopCluster(cl)
  return(sim_comb)
}