## Purpose: Run stan model
## Project: High_Severity_Birds
## Author: Zack Steel

runstan <- function(data_path = "Data/msom_modeldata.RData",
                    samples = 500, #n sample iterations
                    warmup = 500, #n warmup samples
                    nc = 4, #n chains
                    pc = 4, #number of parallel chains (cores)
                    cmdstanr = F,
                    outfile #path for saving the file

) {
  library(tidyverse)
  library(rstan)
  
  ## read in stan_data
  load(data_path)
  
  source('Scripts/stan_utility.R') #http://mc-stan.org/users/documentation/case-studies/rstan_workflow.html
  
  ## Parameters to track
  pars <- c('beta', 'alpha', #occ and detect means
            'b_fire', 'b_patch', 'b_loc', #varying intercepts, mu = 0
            ## community mean effects
            'mu_b_dist', 'mu_b_ysf', 'mu_b_ysfsq', 'mu_b_dy',
            'mu_b_parea', 'mu_b_elev', 'mu_b_elevsq', 
            'mu_a_snagba', 'mu_a_shrub', 'mu_a_jday', 'mu_a_jdaysq',
            ## species-level effects
            'b_dist', 'b_ysf', 'b_ysfsq', 'b_dy', 
            'b_parea', 'b_elev', 'b_elevsq', 'b_lat',
            'a_snagba', 'a_shrub', 'a_jday', 'a_jdaysq',
            'logit_psi'
  )
  
  ## Run using RStan or CmdStanR?
  if(cmdstanr){
    library(cmdstanr)
    
    ## trying cmdstanr https://mc-stan.org/cmdstanr/reference/model-method-sample.html
    m <- cmdstan_model("Code/stan_msom.stan")
    fit <- m$sample(
      data = stan_data,
      seed = 123,
      init = 0,
      chains = nc,
      parallel_chains = pc,
      iter_sampling = samples,
      iter_warmup = warmup
    )
    
    # fit$cmdstan_diagnose()
    
    ## create an rstan object from cmdstanr fit
    out <- rstan::read_stan_csv(fit$output_files())
  } else
  {
    rstan_options(auto_write = TRUE)
    # options(mc.cores = parallel::detectCores())
    # options(mc.cores = pc)
    
    ni <- samples + warmup
    
    ## Call Stan from R
    out <- stan("Code/stan_msom.stan",
                data = stan_data,
                chains = nc, cores = pc,
                iter = ni, warmup = warmup, 
                init = 0, pars = pars, 
                seed = 1, 
                open_progress = FALSE)
  }
  
  ## Save big model
  save(out, stan_data, spp.tab, 
       file = outfile)
  
  # =========================
  # = Timing and Efficiency =
  # =========================
  (timing <- cbind(get_elapsed_time(out), rowSums(get_elapsed_time(out))))
  (max_time <- max(timing))/60
  (neff_md <- median(summary(out)$summary[,"n_eff"], na.rm = T))
  (rhat_md <- median(summary(out)$summary[,"Rhat"], na.rm = T))
  (rhat_max <- max(summary(out)$summary[,"Rhat"], na.rm = T))
  perf <- max_time/neff_md
  
  list(neff_median = neff_md, 
       rhat_median = rhat_md, rhat_max = rhat_max, 
       performance = perf)
}
