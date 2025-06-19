# An example simulation function for E-BiEXNEX model #
# Required R package: R2jags, qpcR, mvtnorm

sim_EBiEXNEX <- function(ntrial, trial_size, muct, muce, st_all, se_all, Th, fcov,
                         w1 = 0.25, w2 = 0.25, w3 = 0.25,
                         num_cl = 2, n_sim = 500, delta = 0, eta = c(0.95, 0.95),
                         model_file = "./E-BiEXNEX.txt") {
  # num_cl: number of threads you want to use
  # n_sim: number of simulations
  # model_file: the location model texts
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
