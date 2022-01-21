library(MASS)
library(rstan)
library(dplyr)
library(tidybayes)
library(parallel)
library(data.table)
library(openblasctl)

# better RNG for parallel processing 
# https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/RngStream.html
# https://www.jottr.org/2020/09/22/push-for-statistical-sound-rng/
RNGkind("L'Ecuyer-CMRG")

# make sure DT and BLAS only use a single thread
setDTthreads(threads = 1)
openblas_set_num_threads(1)

# get Bash args
args <- commandArgs(trailingOnly = TRUE)

# set variables from Bash args
outer_seed      <- as.integer(args[[1]])
sample_out      <- as.numeric(args[[2]])
initial_size    <- as.numeric(args[[3]])
step_size       <- as.numeric(args[[4]])
n_sim           <- as.integer(args[[5]])
stan_chains     <- as.integer(args[[6]])
stan_cores      <- as.integer(args[[7]])
stan_warmup     <- as.integer(args[[8]])
stan_iter       <- as.integer(args[[9]])
mapply_cores    <- as.integer(args[[10]])
pop_rho         <- as.numeric(args[[11]])
model           <- as.character(args[[12]])

# source functions and compile stan files
source("/src/r/functions.R")
source("/src/r/stan_utility.R")

if (model == "weakly_inf") {
    stan_corr_mod <- stan_model(file = paste0("/src/stan/", 
                                              model, 
                                              ".stan"))
} else {
    stan_corr_mod <- stan_model(file = paste0("/src/stan/", 
                                              model, 
                                              "_", 
                                              gsub("\\.", "", as.character(pop_rho)),
                                              ".stan"))
}

# generate data with specified population rho
data <- create_data(pop_rho)

# start simulating
start_time <- Sys.time()
sim_wrapper(paste0(model, "_", gsub("\\.", "", as.character(pop_rho))),
            stan_corr_mod, 
            outer_seed = outer_seed, 
            sample_out = sample_out, 
            initial_size = initial_size, 
            step_size = step_size, 
            n_sim = n_sim,
            stan_chains = stan_chains,
            stan_cores = stan_cores,
            stan_warmup = stan_warmup,
            stan_iter = stan_iter,
            mapply_cores = mapply_cores)
