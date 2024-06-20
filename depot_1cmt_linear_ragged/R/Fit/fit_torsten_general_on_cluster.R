rm(list = ls())
cat("\014")

library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)
library(parallel)
library(furrr)
library(future.apply)
library(future)
library(future.batchtools)
library(here)

set_cmdstan_path("/data/Random/mpi_example1/Torsten/v0.91.0/cmdstan/")
mpi_cmd <- "/usr/bin/mpiexec" ## Default for Metworx

## The approach used for parallel computation on a grid in this example 
## requires at least one active compute node. If no compute node is
## running then wake up a compute node using the qtouch function.
## Wait for the compute node to come up before running the furrr::future_pmap
## step in the script.
qtouch <- function(name = 'touch'){
  system(sprintf('echo "sleep 5" | qsub -N %s', name), intern = TRUE) 
}

nonmem_data <- read_csv("depot_1cmt_linear/Data/depot_1cmt_prop.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

n_subjects <- nonmem_data %>%  # number of individuals
  distinct(ID) %>%
  count() %>%
  deframe()

n_total <- nrow(nonmem_data)   # total number of records

i_obs <- nonmem_data %>%
  mutate(row_num = 1:n()) %>%
  filter(evid == 0) %>%
  select(row_num) %>%
  deframe()

n_obs <- length(i_obs)

subj_start <- nonmem_data %>%
  mutate(row_num = 1:n()) %>%
  group_by(ID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(row_num) %>%
  deframe()

subj_end <- c(subj_start[-1] - 1, n_total)

stan_data <- list(n_subjects = n_subjects,
                  n_total = n_total,
                  n_obs = n_obs,
                  i_obs = i_obs,
                  ID = nonmem_data$ID,
                  amt = nonmem_data$amt,
                  cmt = nonmem_data$cmt,
                  evid = nonmem_data$evid,
                  rate = nonmem_data$rate,
                  ii = nonmem_data$ii,
                  addl = nonmem_data$addl,
                  ss = nonmem_data$ss,
                  time = nonmem_data$time,
                  dv = nonmem_data$DV,
                  subj_start = subj_start,
                  subj_end = subj_end,
                  lloq = nonmem_data$lloq,
                  bloq = nonmem_data$bloq,
                  location_tvcl = 0.6,
                  location_tvvc = 15,
                  location_tvka = 1,
                  scale_tvcl = 1,
                  scale_tvvc = 1,
                  scale_tvka = 1,
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  scale_omega_ka = 0.4,
                  lkj_df_omega = 2,
                  scale_sigma_p = 0.5,
                  prior_only = 0)

model <- cmdstan_model(
  "depot_1cmt_linear/Stan/Fit/depot_1cmt_prop_torsten_general.stan")

## number of cores per chain
## Each Metworx compute node should have a number of vCPUs >= (2 * n_jobs)
n_chains <- 4

RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()

model_name <- "depot_1cmt_linear_torsten_general"

sample_and_save_all_mpi <- function(chain,
                                    seed = 112356,
                                    run_number,
                                    n_jobs,
                                    stan_data,
                                    model,
                                    iter_warmup = 500,
                                    iter_sampling = 1000,
                                    output_basename,
                                    output_dir,
                                    max_treedepth = 10,
                                    adapt_delta = 0.80,
                                    refresh = 100,
                                    save_warmup = FALSE){
  
  set_cmdstan_path("/data/Random/mpi_example1/Torsten/v0.91.0/cmdstan/")
  
  ## Setup use of future with Grid Engine
  sge <- future::tweak(
    future.batchtools::batchtools_sge,
    label = model_name,
    template = "/data/Random/mpi_example1/batchtools.sge-mrg.tmpl",
    workers = n_chains, # This is really the number of "jobs" that we want to create; no clear relationship to number of worker nodes
    resources = list(slots = 2 * n_jobs) # this is the number of cores that each job will require
  )
  
  future::plan(sge)
  
  init_files <- str_c("depot_1cmt_linear/Data/Inits/inits_", run_number, "_",
                      chain, ".json")
  
  fit <- model$sample_mpi(data = stan_data,
                          mpi_cmd = mpi_cmd,
                          mpi_args = list("bind-to" = "core",
                                          "n" = n_jobs),
                          seed = seed,
                          chains = 1,
                          chain_ids = chain,
                          iter_warmup = iter_warmup,
                          iter_sampling = iter_sampling,
                          init = init_files,
                          output_dir = output_dir, 
                          max_treedepth = max_treedepth,
                          adapt_delta = adapt_delta,
                          refresh = refresh,
                          save_warmup = save_warmup)
  
  fit$save_output_files(dir = output_dir, 
                        basename = str_c(output_basename, run_number, "_", 
                                         chain),
                        timestamp = FALSE, random = FALSE)
  
}

# qtouch()

go_through_them_all <- function(n_jobs, run_number){
  
  furrr::future_pwalk(list(chain = 1:n_chains,
                           seed = 10^c(0:3)), 
                      .f = sample_and_save_all_mpi,
                      run_number = run_number,
                      n_jobs = n_jobs,
                      stan_data = stan_data,
                      model = model,
                      iter_warmup = 500,
                      iter_sampling = 1000,
                      output_dir = "depot_1cmt_linear/Stan/Fits/MPI/",
                      output_basename = str_c("torsten_general_",
                                              n_jobs, "_jobs_run_"),
                      max_treedepth = 10,
                      adapt_delta = 0.80,
                      refresh = 100,
                      save_warmup = FALSE,
                      .options = furrr_options(seed = TRUE))
  
  fit <- as_cmdstan_fit(
    str_c("depot_1cmt_linear/Stan/Fits/MPI/torsten_general_", n_jobs, 
          "_jobs_run_", run_number, "_", 1:4, "-1.csv"))
  
  fit$save_object(
    str_c("depot_1cmt_linear/Stan/Fits/MPI/torsten_general_", n_jobs, 
          "_jobs_run_", run_number, ".rds"))
}

expand_grid(n_jobs = c(1, 2, 4, 8, 12, 24, 48),
            run_number = 1:10) %>%
  pwalk(.f = go_through_them_all)
