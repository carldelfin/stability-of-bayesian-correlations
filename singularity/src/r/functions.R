# ==================================================================================================
# Create data with specific correlation
# ==================================================================================================

create_data <- function(rho) {
    set.seed(outer_seed)
    data <- as.data.table(mvrnorm(1000000,
                                  mu = c(0, 0),
                                  Sigma = matrix(c(1, rho, rho, 1), ncol = 2),
                                  empirical = TRUE)) %>%
            setnames(., c("V1", "V2"), c("x", "y"))
    data[, id := seq(1, 1000000, 1)]
    return(data)
}

# ==================================================================================================
# Main simulation function
# ==================================================================================================

correlation_sim <- function(mod_name, 
                            mod, 
                            outer_seed, 
                            inner_seed,
                            sample_out, 
                            initial_size,
                            step_size,
                            n_sim,
                            stan_chains,
                            stan_cores,
                            stan_warmup,
                            stan_iter,
                            mapply_cores) {
    # start timer
    step_start <- Sys.time()

    # create storage matrix
    # allocating with NA_real_ should be faster than using NAs:
    # https://www.r-bloggers.com/2014/06/pitfall-did-you-really-mean-to-use-matrixnrow-ncol/
    matrix_length <- 1 + ((sample_out - initial_size) / step_size)
    M = matrix(NA_real_, nrow = matrix_length, ncol = 27)
   
    # prepare data
    set.seed(inner_seed)
    initial_sample <- data[sample(.N, initial_size)]
    remaining_sample <- data[!data$id %in% initial_sample$id, ]
    data_for_stan <- list(x = initial_sample[[1]],
                          y = initial_sample[[2]],
                          N = nrow(initial_sample))
    
    # begin sampling
    tmp <- sampling(mod,
                    data = data_for_stan,
                    chains = stan_chains,
                    iter = stan_iter,
                    warmup = stan_warmup,
                    cores = stan_cores,
                    seed = as.integer(inner_seed + 10000),
                    init_r = 1,
                    #control = list(adapt_delta = 0.99),
                    verbose = FALSE,
                    show_messages = FALSE,
                    refresh = 0)
    
    # MCMC diagnostics
    # see stan_utility.R for details, taken from:
    # https://github.com/betanalpha/knitr_case_studies/blob/master/rstan_workflow/stan_utility.R
    num_div <- get_num_divergent(tmp)
    num_treedepth <- get_num_max_treedepth(tmp)
    num_bfmi <- get_low_bfmi_chains(tmp)
    num_neff <- as.numeric(1 - as.integer(check_n_eff(tmp, quiet = TRUE)))
    num_rhat <- as.numeric(1 - as.integer(check_rhat(tmp, quiet = TRUE)))
    
    # sum will be > 0 if there are any issues 
    sum_diagnostics <- sum(num_div + num_treedepth + num_bfmi + num_neff + num_rhat)

    # ESS values for beta (i.e., estimated rho)
    # see: https://avehtari.github.io/rhat_ess/rhat_ess.html
    tmp_diag <- monitor(tmp, print = FALSE)
    tail_ess <- as.numeric(tmp_diag[1, "Tail_ESS"])
    bulk_ess <- as.numeric(tmp_diag[1, "Bulk_ESS"])

    # extract posterior samples
    post <- rstan::extract(tmp,
                           permuted = TRUE,
                           inc_warmup = FALSE,
                           include = TRUE) %>%
        as.data.frame()

    # freq model
    tmp_freq <- cor.test(initial_sample[[1]], initial_sample[[2]])

    # point estimates
    M[1, 1] <- mean_hdi(post[, "beta"], .width = 0.95)[, "y"]
    M[1, 2] <- median_hdi(post[, "beta"], .width = 0.95)[, "y"]
    M[1, 3] <- mode_hdi(post[, "beta"], .width = 0.95)[, "y"]
    M[1, 4] <- cor.test(initial_sample[[1]], initial_sample[[2]])$estimate[[1]]
   
    # intervals, hdi
    M[1, 5] <- mean_hdi(post[, "beta"], .width = 0.95)[, "ymin"]
    M[1, 6] <- mean_hdi(post[, "beta"], .width = 0.95)[, "ymax"]
    M[1, 7] <- mean_hdi(post[, "beta"], .width = 0.90)[, "ymin"]
    M[1, 8] <- mean_hdi(post[, "beta"], .width = 0.90)[, "ymax"]
    M[1, 9] <- mean_hdi(post[, "beta"], .width = 0.66)[, "ymin"]
    M[1, 10] <- mean_hdi(post[, "beta"], .width = 0.66)[, "ymax"]
    
    # intervals, qi
    M[1, 11] <- mean_qi(post[, "beta"], .width = 0.95)[, "ymin"]
    M[1, 12] <- mean_qi(post[, "beta"], .width = 0.95)[, "ymax"]
    M[1, 13] <- mean_qi(post[, "beta"], .width = 0.90)[, "ymin"]
    M[1, 14] <- mean_qi(post[, "beta"], .width = 0.90)[, "ymax"]
    M[1, 15] <- mean_qi(post[, "beta"], .width = 0.66)[, "ymin"]
    M[1, 16] <- mean_qi(post[, "beta"], .width = 0.66)[, "ymax"]

    # intervals, frequentist
    M[1, 17] <- cor.test(initial_sample[[1]], initial_sample[[2]], conf.level = 0.95)$conf.int[[1]]
    M[1, 18] <- cor.test(initial_sample[[1]], initial_sample[[2]], conf.level = 0.95)$conf.int[[2]]
    M[1, 19] <- cor.test(initial_sample[[1]], initial_sample[[2]], conf.level = 0.90)$conf.int[[1]]
    M[1, 20] <- cor.test(initial_sample[[1]], initial_sample[[2]], conf.level = 0.90)$conf.int[[2]]
    M[1, 21] <- cor.test(initial_sample[[1]], initial_sample[[2]], conf.level = 0.66)$conf.int[[1]]
    M[1, 22] <- cor.test(initial_sample[[1]], initial_sample[[2]], conf.level = 0.66)$conf.int[[2]]

    # sample size and diagnostics
    M[1, 23] <- nrow(initial_sample)
    M[1, 24] <- sum_diagnostics 
    M[1, 25] <- bulk_ess 
    M[1, 26] <- tail_ess 
    
    # being looping through rest of sample
    for (i in 2:matrix_length) {
        set.seed(inner_seed)
        initial_sample <- rbind(initial_sample, 
                                dplyr::sample_n(remaining_sample, step_size))
        
        data_for_stan <- list(x = initial_sample[[1]],
                              y = initial_sample[[2]],
                              N = nrow(initial_sample))
    
        tmp <- sampling(mod,
                        data = data_for_stan,
                        chains = stan_chains,
                        iter = stan_iter,
                        warmup = stan_warmup,
                        cores = stan_cores,
                        seed = as.integer(inner_seed + 10000),
                        init_r = 1,
                        verbose = FALSE,
                        show_messages = FALSE,
                        refresh = 0)

        # diagnostics
        num_div <- get_num_divergent(tmp)
        num_treedepth <- get_num_max_treedepth(tmp)
        num_bfmi <- get_low_bfmi_chains(tmp)
        num_neff <- as.numeric(1 - as.integer(check_n_eff(tmp, quiet = TRUE)))
        num_rhat <- as.numeric(1 - as.integer(check_rhat(tmp, quiet = TRUE)))

        sum_diagnostics <- sum(num_div + num_treedepth + num_bfmi + num_neff + num_rhat)
        
        tmp_diag <- monitor(tmp, print = FALSE)
        tail_ess <- as.numeric(tmp_diag[1, "Tail_ESS"])
        bulk_ess <- as.numeric(tmp_diag[1, "Bulk_ESS"])
        
        post <- rstan::extract(tmp,
                               permuted = TRUE, 
                               inc_warmup = FALSE, 
                               include = TRUE) %>%
            as.data.frame()
        
        # point estimates
        M[i, 1] <- mean_hdi(post[, "beta"], .width = 0.95)[, "y"]
        M[i, 2] <- median_hdi(post[, "beta"], .width = 0.95)[, "y"]
        M[i, 3] <- mode_hdi(post[, "beta"], .width = 0.95)[, "y"]
        M[i, 4] <- cor.test(initial_sample[[1]], initial_sample[[2]])$estimate[[1]]
       
        # intervals, hdi
        M[i, 5] <- mean_hdi(post[, "beta"], .width = 0.95)[, "ymin"]
        M[i, 6] <- mean_hdi(post[, "beta"], .width = 0.95)[, "ymax"]
        M[i, 7] <- mean_hdi(post[, "beta"], .width = 0.90)[, "ymin"]
        M[i, 8] <- mean_hdi(post[, "beta"], .width = 0.90)[, "ymax"]
        M[i, 9] <- mean_hdi(post[, "beta"], .width = 0.66)[, "ymin"]
        M[i, 10] <- mean_hdi(post[, "beta"], .width = 0.66)[, "ymax"]
        
        # intervals, qi
        M[i, 11] <- mean_qi(post[, "beta"], .width = 0.95)[, "ymin"]
        M[i, 12] <- mean_qi(post[, "beta"], .width = 0.95)[, "ymax"]
        M[i, 13] <- mean_qi(post[, "beta"], .width = 0.90)[, "ymin"]
        M[i, 14] <- mean_qi(post[, "beta"], .width = 0.90)[, "ymax"]
        M[i, 15] <- mean_qi(post[, "beta"], .width = 0.66)[, "ymin"]
        M[i, 16] <- mean_qi(post[, "beta"], .width = 0.66)[, "ymax"]

        # intervals, frequentist
        M[i, 17] <- cor.test(initial_sample[[1]], initial_sample[[2]], conf.level = 0.95)$conf.int[[1]]
        M[i, 18] <- cor.test(initial_sample[[1]], initial_sample[[2]], conf.level = 0.95)$conf.int[[2]]
        M[i, 19] <- cor.test(initial_sample[[1]], initial_sample[[2]], conf.level = 0.90)$conf.int[[1]]
        M[i, 20] <- cor.test(initial_sample[[1]], initial_sample[[2]], conf.level = 0.90)$conf.int[[2]]
        M[i, 21] <- cor.test(initial_sample[[1]], initial_sample[[2]], conf.level = 0.66)$conf.int[[1]]
        M[i, 22] <- cor.test(initial_sample[[1]], initial_sample[[2]], conf.level = 0.66)$conf.int[[2]]

        # sample size and diagnostics
        M[i, 23] <- nrow(initial_sample)
        M[i, 24] <- sum_diagnostics 
        M[i, 25] <- bulk_ess 
        M[i, 26] <- tail_ess 

        # get remaining sample and add 1 to index
        remaining_sample <- data[!data$id %in% initial_sample$id, ]
    }

    # note which seed this was
    M[, 27] <- inner_seed

    # some rudimentary timing information 
    # it's not super accurate but will give a rough idea of elapsed/remaining time
    end_time <- Sys.time()
    step_time <- difftime(end_time, step_start, units = "secs")
    elapsed_time <- difftime(end_time, start_time, units = "secs") 
    est_total <- (step_time[[1]] * n_sim) / mapply_cores
    est_remaining <- est_total - elapsed_time[[1]]

    cat("\n")
    cat("Current model:", mod_name, "\n")
    cat("Completed seed:", inner_seed, "\n")
    cat("This step took:", round(step_time / 60, 2), "minutes", "\n")
    cat("Time since start:", round(elapsed_time / 60, 2), "minutes", "\n")
    cat("Estimated total time:", round(est_total / 60, 2), "minutes", "\n")
    cat("Estimated time remaining:", round(est_remaining / 60, 2), "minutes", "\n")
    
    rm(end_time, step_time, elapsed_time, est_total, est_remaining, tmp, post)
    return(M)

}

# ==================================================================================================
# Simulation wrapper
# ==================================================================================================

sim_wrapper <- function(mod_name, 
                        mod, 
                        outer_seed, 
                        sample_out, 
                        initial_size,
                        step_size,
                        n_sim, 
                        stan_chains,
                        stan_cores,
                        stan_warmup,
                        stan_iter,
                        mapply_cores) {
    out <- mcmapply(FUN = correlation_sim,
             inner_seed = 1:n_sim,
             MoreArgs = list(mod_name = mod_name,
                             mod = mod,
                             outer_seed = outer_seed,
                             sample_out = sample_out, 
                             initial_size = initial_size, 
                             step_size = step_size, 
                             n_sim = n_sim,
                             stan_chains = stan_chains,
                             stan_cores = stan_cores,
                             stan_warmup = stan_warmup,
                             stan_iter = stan_iter,
                             mapply_cores = mapply_cores),
             mc.cores = mapply_cores,
             SIMPLIFY = FALSE)

    # clean up 
    out <- as.data.frame(do.call(rbind, out)) %>%
         setNames(., c("mean_rho",
                       "median_rho",
                       "mode_rho",
                       "freq_rho",

                       "mean_hdi_lower_95",
                       "mean_hdi_upper_95",
                       "mean_hdi_lower_90",
                       "mean_hdi_upper_90",
                       "mean_hdi_lower_66",
                       "mean_hdi_upper_66",
                       
                       "mean_qi_lower_95",
                       "mean_qi_upper_95",
                       "mean_qi_lower_90",
                       "mean_qi_upper_90",
                       "mean_qi_lower_66",
                       "mean_qi_upper_66",
                       
                       "freq_lower_95",
                       "freq_upper_95",
                       "freq_lower_90",
                       "freq_upper_90",
                       "freq_lower_66",
                       "freq_upper_66",

                       "sample_size",
                       "diagnostic_error",
                       "bulk_ess",
                       "tail_ess",
                       "seed"))

    # save
    saveRDS(out, file = paste0("/output/", mod_name, ".rds"))
    rm(out)
}
